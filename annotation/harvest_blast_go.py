#!/usr/bin/env python3
"""
Assign GO terms to Deinandra genes from the NCBI Gene annotation of each
transcript's best BLASTp hit, for one eggNOG bucket selected by --mode:
  no_go  (default) transcripts with an eggNOG hit but no GO term
  has_go           transcripts eggNOG already gave GO terms; output is a
                   supplement to be unioned with the eggNOG terms

Steps: pick target transcripts from the eggNOG annotations; take each one's
best (lowest e-value) BLASTp hit; resolve accession -> NCBI GeneID via elink;
look up GeneIDs in NCBI's bulk gene2go file; map transcripts to genes via the
GFF and write gene_id/go_term.

elink results are appended to ACC2GENE_CACHE after each batch, so a rerun
resumes from the cache; only failed requests are retried.

Usage: python3 harvest_blast_go.py [--mode no_go|has_go] [--limit N] [output.tsv]
  --mode M  : which eggNOG GOs bucket to target (default: no_go)
  --limit N : only process the first N target transcripts (for testing)
  output.tsv defaults to annotated_sequences_blastharvest_geneid_go.tsv for
  --mode no_go, and annotated_sequences_blastharvest_supplement_geneid_go.tsv
  for --mode has_go.
"""
import sys
import os
import re
import gzip
import time
import urllib.request
import urllib.parse
import urllib.error
from collections import defaultdict

BASE_DIR = "/private/groups/meyerlab/esteiger"

GFF = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
EGGNOG_ANN = f"{BASE_DIR}/genome/annotation/gvtp_eggnog_Di2_viridiplantae.emapper.annotations"
BLASTP = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.blastp.out"
GENE2GO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
GENE2GO_LOCAL = f"{BASE_DIR}/genome/annotation/gene2go.gz"
DEFAULT_OUT = {
    "no_go": f"{BASE_DIR}/genome/annotated_sequences_blastharvest_geneid_go.tsv",
    "has_go": f"{BASE_DIR}/genome/annotated_sequences_blastharvest_supplement_geneid_go.tsv",
}

# Accession -> GeneID cache, shared by both modes. Holds confirmed results only
# (a GeneID or confirmed "no linked gene"); failed requests are never cached.
ACC2GENE_CACHE = f"{BASE_DIR}/genome/annotation/harvest_blast_go_acc2gene_cache.tsv"

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "gvtp_annotation_pipeline"
EMAIL = "emma.c.steigerwald@gmail.com"
# With separate id= params, elink returns an empty response above ~30 ids
BATCH_SIZE = 25
SLEEP_BETWEEN_BATCHES = 0.4  # stay under the 3 req/sec unauthenticated limit


def load_tx2gene():
    tx2gene = {}
    with open(GFF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "mRNA":
                continue
            parent = re.search(r"Parent=([^;]+)", fields[8])
            mid = re.search(r"ID=([^;]+)", fields[8])
            if parent and mid:
                tx2gene[mid.group(1)] = parent.group(1)
    return tx2gene


def find_target_transcripts(mode):
    """eggNOG query rows matching mode's GOs filter:
    no_go  -> GOs == '-'  (eggNOG assigned zero GO terms)
    has_go -> GOs != '-'  (eggNOG assigned at least one GO term already)"""
    targets = []
    header = None
    with open(EGGNOG_ANN) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                continue
            if header is None:
                continue
            cols = line.split("\t")
            if len(cols) < len(header):
                continue
            row = dict(zip(header, cols))
            has_go = row.get("GOs", "-") != "-"
            if (mode == "has_go") == has_go:
                targets.append(row["query"])
    return targets


def best_blast_hits(transcripts):
    """transcript -> best (lowest evalue) accession, scanning the blastp file once."""
    wanted = set(transcripts)
    best = {}  # transcript -> (evalue, accession)
    with open(BLASTP) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            query = cols[0]
            if query not in wanted:
                continue
            acc = cols[1]
            try:
                evalue = float(cols[10])
            except ValueError:
                continue
            if query not in best or evalue < best[query][0]:
                best[query] = (evalue, acc)
    return {t: acc for t, (ev, acc) in best.items()}


REQUEST_FAILED = "__REQUEST_FAILED__"  # sentinel: distinct from a confirmed "no linked gene"


def elink_batch(accessions, max_retries=4):
    """POST a batch of protein accessions to elink, one id= param each so
    linksets keep per-input order (comma-joined ids are pooled). Returns
    {accession: geneid | None (no linked gene) | REQUEST_FAILED}."""
    import json
    import http.client

    params = [("dbfrom", "protein"), ("db", "gene"), ("retmode", "json"),
              ("tool", TOOL), ("email", EMAIL)]
    for acc in accessions:
        params.append(("id", acc))
    data = urllib.parse.urlencode(params).encode()

    result = None
    for attempt in range(1, max_retries + 1):
        req = urllib.request.Request(f"{EUTILS_BASE}/elink.fcgi", data=data)
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                result = json.loads(resp.read().decode())
            break
        except (urllib.error.URLError, urllib.error.HTTPError,
                http.client.HTTPException, ConnectionError, TimeoutError,
                json.JSONDecodeError, UnicodeDecodeError) as e:
            # E-utilities can return a 200 with a truncated body under load
            wait = 2 * attempt
            print(f"  WARNING: elink batch failed ({e}), attempt {attempt}/{max_retries}, "
                  f"retrying in {wait}s..." if attempt < max_retries else
                  f"  WARNING: elink batch failed ({e}), giving up after {max_retries} attempts",
                  file=sys.stderr)
            if attempt < max_retries:
                time.sleep(wait)

    if result is None:
        return {acc: REQUEST_FAILED for acc in accessions}

    linksets = result.get("linksets", [])
    resolved = {}
    # linksets are returned in submission order, one per input id
    for acc, ls in zip(accessions, linksets):
        links = None
        for ldb in ls.get("linksetdbs", []):
            if ldb.get("dbto") == "gene":
                links = ldb.get("links")
                break
        resolved[acc] = links[0] if links else None  # confirmed: no linked gene
    for acc in accessions:
        resolved.setdefault(acc, None)  # elink omits ids with zero links entirely
    return resolved


def load_acc2gene_cache():
    """accession -> GeneID (or None for a confirmed no-linked-gene) from prior runs."""
    cache = {}
    if not os.path.exists(ACC2GENE_CACHE):
        return cache
    with open(ACC2GENE_CACHE) as f:
        next(f, None)  # header
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            acc, gid = cols[0], cols[1]
            cache[acc] = gid if gid else None
    return cache


def append_to_cache(new_entries):
    """Append confirmed (acc -> geneid-or-None) results and flush immediately."""
    if not new_entries:
        return
    is_new_file = not os.path.exists(ACC2GENE_CACHE)
    with open(ACC2GENE_CACHE, "a") as f:
        if is_new_file:
            f.write("accession\tgeneid\n")
        for acc, gid in new_entries.items():
            f.write(f"{acc}\t{gid or ''}\n")
        f.flush()
        os.fsync(f.fileno())


def resolve_accessions_to_geneids(accessions, max_rounds=6):
    """Resolve accessions via elink, skipping those in ACC2GENE_CACHE and
    retrying failed requests for up to max_rounds rounds."""
    cache = load_acc2gene_cache()
    resolved = {}
    todo = []
    for acc in sorted(set(accessions)):
        if acc in cache:
            resolved[acc] = cache[acc]
        else:
            todo.append(acc)
    if len(accessions) - len(todo):
        print(f"  {len(accessions) - len(todo)}/{len(accessions)} accessions already resolved "
              f"in a prior run's cache ({ACC2GENE_CACHE}) — skipping elink for those.", file=sys.stderr)

    for round_num in range(1, max_rounds + 1):
        if not todo:
            break
        n_batches = (len(todo) + BATCH_SIZE - 1) // BATCH_SIZE
        label = "" if round_num == 1 else f" [retry round {round_num}]"
        for i in range(0, len(todo), BATCH_SIZE):
            batch = todo[i:i + BATCH_SIZE]
            batch_num = i // BATCH_SIZE + 1
            print(f"  elink batch {batch_num}/{n_batches}{label} ({len(batch)} accessions)...", file=sys.stderr)
            batch_result = elink_batch(batch)
            resolved.update(batch_result)
            newly_confirmed = {acc: gid for acc, gid in batch_result.items() if gid != REQUEST_FAILED}
            append_to_cache(newly_confirmed)
            time.sleep(SLEEP_BETWEEN_BATCHES)

        todo = [acc for acc, v in resolved.items() if v == REQUEST_FAILED]
        if todo:
            print(f"  {len(todo)} accession(s) still unresolved after round {round_num} "
                  f"(request failures, not confirmed no-gene) — retrying...", file=sys.stderr)

    still_failed = [acc for acc, v in resolved.items() if v == REQUEST_FAILED]
    if still_failed:
        print(f"  WARNING: {len(still_failed)} accession(s) never got a usable elink response "
              f"after {max_rounds} rounds — treating as unresolved (not as confirmed no-gene): "
              f"{still_failed[:10]}{'...' if len(still_failed) > 10 else ''}", file=sys.stderr)
    return resolved, still_failed


def ensure_gene2go():
    if os.path.exists(GENE2GO_LOCAL) and os.path.getsize(GENE2GO_LOCAL) > 0:
        print(f"gene2go already present at {GENE2GO_LOCAL}", file=sys.stderr)
        return
    print(f"Downloading {GENE2GO_URL} -> {GENE2GO_LOCAL} (~1.3GB, one-time)...", file=sys.stderr)
    req = urllib.request.Request(GENE2GO_URL, headers={"User-Agent": f"{TOOL} ({EMAIL})"})
    with urllib.request.urlopen(req, timeout=300) as resp, open(GENE2GO_LOCAL, "wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    print("Download complete.", file=sys.stderr)


def load_gene2go(gene_ids):
    """Stream gene2go.gz, keeping only the requested GeneIDs (any taxon)."""
    wanted = set(gene_ids)
    gene2go = defaultdict(set)
    with gzip.open(GENE2GO_LOCAL, "rt") as f:
        header = f.readline()
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            gid, go_id = cols[1], cols[2]
            if gid in wanted:
                gene2go[gid].add(go_id)
    return gene2go


def main():
    args = sys.argv[1:]
    mode = "no_go"
    if "--mode" in args:
        idx = args.index("--mode")
        mode = args[idx + 1]
        del args[idx:idx + 2]
    if mode not in DEFAULT_OUT:
        sys.exit(f"--mode must be one of {sorted(DEFAULT_OUT)}, got {mode!r}")
    limit = None
    if "--limit" in args:
        idx = args.index("--limit")
        limit = int(args[idx + 1])
        del args[idx:idx + 2]
    out_file = args[0] if args else DEFAULT_OUT[mode]

    print("Loading transcript -> gene map from GFF...", file=sys.stderr)
    tx2gene = load_tx2gene()

    print(f"Finding eggNOG '{mode}' transcripts...", file=sys.stderr)
    targets = find_target_transcripts(mode)
    print(f"Found {len(targets)} target transcripts", file=sys.stderr)
    if limit:
        targets = targets[:limit]
        print(f"--limit set: using first {len(targets)}", file=sys.stderr)

    print("Finding best BLAST hit per target transcript...", file=sys.stderr)
    best_hits = best_blast_hits(targets)
    print(f"{len(best_hits)}/{len(targets)} targets have a BLAST hit", file=sys.stderr)

    accessions = sorted(set(best_hits.values()))
    print(f"Resolving {len(accessions)} unique accessions to GeneIDs via elink...", file=sys.stderr)
    acc2gene, unresolved_accessions = resolve_accessions_to_geneids(accessions)
    n_resolved = sum(1 for v in acc2gene.values() if v and v != REQUEST_FAILED)
    n_confirmed_no_gene = sum(1 for v in acc2gene.values() if v is None)
    print(f"{n_resolved}/{len(accessions)} accessions resolved to a GeneID "
          f"({n_confirmed_no_gene} confirmed no linked gene, "
          f"{len(unresolved_accessions)} still unresolved after retries — see WARNING above)", file=sys.stderr)
    if unresolved_accessions:
        unresolved_file = os.path.join(os.path.dirname(os.path.abspath(out_file)) or ".",
                                        "harvest_blast_go_unresolved_accessions.txt")
        with open(unresolved_file, "w") as f:
            f.write("\n".join(unresolved_accessions) + "\n")
        print(f"Wrote {len(unresolved_accessions)} unresolved accessions to {unresolved_file} "
              f"for a future retry pass.", file=sys.stderr)

    ensure_gene2go()
    gene_ids_needed = {gid for gid in acc2gene.values() if gid and gid != REQUEST_FAILED}
    print(f"Looking up GO terms for {len(gene_ids_needed)} GeneIDs in gene2go...", file=sys.stderr)
    ncbi_gene2go = load_gene2go(gene_ids_needed)
    print(f"{len(ncbi_gene2go)}/{len(gene_ids_needed)} GeneIDs had at least one GO term", file=sys.stderr)

    pairs = set()
    n_unmapped_tx = 0
    for tx, acc in best_hits.items():
        our_gene = tx2gene.get(tx)
        if our_gene is None:
            our_gene = tx.split(".")[0]
            n_unmapped_tx += 1
        ncbi_gene = acc2gene.get(acc)
        if not ncbi_gene or ncbi_gene == REQUEST_FAILED:
            continue
        for go_id in ncbi_gene2go.get(ncbi_gene, ()):
            pairs.add((our_gene, go_id))

    pairs = sorted(pairs)
    with open(out_file, "w") as out:
        out.write("gene_id\tgo_term\n")
        for gene_id, go_id in pairs:
            out.write(f"{gene_id}\t{go_id}\n")

    print(f"({n_unmapped_tx} target transcripts not found in GFF, used dot-strip fallback)", file=sys.stderr)
    print(f"Wrote {len(pairs)} gene-GO pairs across {len(set(g for g, _ in pairs))} genes to {out_file}")


if __name__ == "__main__":
    main()
