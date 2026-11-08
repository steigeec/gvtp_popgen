#!/usr/bin/env python3
"""
Builds the Orthologs, Ortholog_Source, Acronym, Function_Description and
BLAST_Gene_Names columns used by compile_gene_summary.py.

For each gene, in priority order:
  1. Best BLASTP hit from an Asteraceae species (exact genus match to
     ASTERACEAE_GENERA) with evalue <= ASTERACEAE_EVALUE_MAX; the hit
     description is used directly.
  2. Else, the OrthoFinder (Results_Jul18) Arabidopsis_thaliana ortholog,
     with symbol and description from the TAIR10 peptide FASTA headers.
  3. Else, the eggNOG-mapper Description/Preferred_name (Di2_viridiplantae run).
  4. Else, Orthologs is empty.

BLAST_Gene_Names is the unrestricted BLASTP top hit, written only when
steps 1-3 find nothing.

Usage:
  python3 build_gene_orthology.py <gene_list_file> --out orthologs_<var>.txt
Output columns: Gene_ID, Orthologs, Ortholog_Source, BLAST_Gene_Names
"""
import argparse
import re
import sys
from collections import defaultdict

sys.path.insert(0, f"{BASE_DIR}/gvtp_popgen/selection/gemma")
from curate_intergenic_by_ld import load_gff_resolver

BASE_DIR = "/private/groups/meyerlab/esteiger"

BLASTP_OUT = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.blastp.out"
ORTHOFINDER_DEINANDRA_TSV = (f"{BASE_DIR}/genome/annotation/orthofinder_results/"
                             "Results_Jul18/Orthologues/Deinandra.tsv")
TAIR_PEP_FASTA = f"{BASE_DIR}/genome/annotation/Arabidopsis_thaliana.TAIR10.pep.all.fa"
EGGNOG_ANNOTATIONS = f"{BASE_DIR}/genome/annotation/gvtp_eggnog_Di2_viridiplantae.emapper.annotations"

ASTERACEAE_EVALUE_MAX = 1e-50

# Asteraceae genera present in BLASTP_OUT, matched exactly (a substring match on
# "aster" would include e.g. Cotoneaster and Acanthaster). No Madiinae genera are in the database.
ASTERACEAE_GENERA = {
    "Ageratina", "Ambrosia", "Archibaccharis", "Arctium", "Artemisia",
    "Aster", "Baccharis", "Bidens", "Carthamus", "Centaurea",
    "Chrysanthemum", "Cichorium", "Cirsium", "Cosmos", "Cynara",
    "Dendrosenecio", "Echinacea", "Erigeron", "Gerbera", "Guizotia",
    "Helianthus", "Ixeris", "Lactuca", "Parthenium", "Pluchea",
    "Saussurea", "Senecio", "Sinosenecio", "Solidago", "Sonchus",
    "Sphagneticola", "Symphyotrichum", "Tagetes", "Tanacetum",
    "Taraxacum", "Xanthium", "Youngia",
}


def load_gene2tx():
    """gene_id -> [transcript_id, ...], inverse of curate_intergenic_by_ld.py's tx2gene."""
    _, tx2gene = load_gff_resolver()
    gene2tx = defaultdict(list)
    for tx, gid in tx2gene.items():
        gene2tx[gid].append(tx)
    return gene2tx


def load_blast_hits():
    """transcript_id -> list of (evalue, stitle), best-first."""
    hits = defaultdict(list)
    with open(BLASTP_OUT) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 14:
                continue
            qseqid, evalue, stitle = cols[0], float(cols[10]), cols[13]
            hits[qseqid].append((evalue, stitle))
    for tx in hits:
        hits[tx].sort(key=lambda x: x[0])
    return hits


ORGANISM_RE = re.compile(r'\[([A-Z][a-z]+)\s')


def best_asteraceae_hit(tx_hits):
    """Best (stitle, genus) from any ASTERACEAE_GENERA species clearing
    ASTERACEAE_EVALUE_MAX, or None. Matches on the exact genus (first word after '[')."""
    for evalue, stitle in tx_hits:
        if evalue > ASTERACEAE_EVALUE_MAX:
            break  # sorted by evalue
        m = ORGANISM_RE.search(stitle)
        if m and m.group(1) in ASTERACEAE_GENERA:
            return stitle, m.group(1)
    return None


ACRONYM_TOKEN_RE = re.compile(r'\b([A-Z][A-Za-z0-9]{1,9})\b')
ISOFORM_TOKEN_RE = re.compile(r'^X\d+$')


def extract_acronym(orthologs_text):
    """Extract a gene symbol/locus id from the Orthologs hit text
    (e.g. "...mitochondrial acidic protein MAM33 [Helianthus annuus]" gives "MAM33").
    Only uses the Orthologs text, since Function_Description may come from another species.
    The token must contain a digit, to exclude all-caps family names (GPI, ABC).
    Returns None if nothing qualifies."""
    if orthologs_text is None:
        return None
    t = re.sub(r'^[A-Z]{2}_\d+\.\d+\s*', '', orthologs_text)
    t = re.sub(r'\s*\[.*?\]\s*$', '', t)
    t = re.sub(r',?\s*isoform\s+X\d+\s*$', '', t, flags=re.I)
    t = re.sub(r',?\s*putative\s*$', '', t, flags=re.I)
    candidates = [
        m.group(1) for m in ACRONYM_TOKEN_RE.finditer(t)
        if not ISOFORM_TOKEN_RE.match(m.group(1))
        and not m.group(1).startswith("LOC")
        and any(c.isdigit() for c in m.group(1))
        and len(m.group(1)) >= 3
    ]
    return candidates[-1] if candidates else None


def load_orthofinder_arabidopsis():
    """Deinandra transcript_id -> [Arabidopsis gene_id, ...]"""
    tx2at = defaultdict(list)
    with open(ORTHOFINDER_DEINANDRA_TSV) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(idx.values()):
                continue
            if cols[idx["Species"]] != "Arabidopsis_thaliana":
                continue
            deinandra_txs = [t.strip() for t in cols[idx["Deinandra"]].split(",")]
            at_genes = [g.strip() for g in cols[idx["Orthologs"]].split(",")]
            for tx in deinandra_txs:
                tx2at[tx].extend(at_genes)
    return tx2at


def load_tair_descriptions():
    """Arabidopsis gene_id (AT#G#####, no transcript suffix) -> (symbol_or_None, description)."""
    desc = {}
    header_re = re.compile(
        r"^>(\S+)\s+pep\s+.*?gene:(\S+)\s+.*?(?:gene_symbol:(\S+)\s+)?description:(.+?)\s*\[Source"
    )
    with open(TAIR_PEP_FASTA) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m = header_re.search(line)
            if not m:
                continue
            gene_id, symbol, description = m.group(2), m.group(3), m.group(4)
            if gene_id in desc:
                continue
            # gene_symbol equal to the AT number means no symbol is assigned
            real_symbol = symbol if (symbol and symbol != gene_id) else None
            desc[gene_id] = (real_symbol, description)
    return desc


def load_eggnog():
    """Deinandra transcript_id -> (symbol_or_None, description)."""
    desc = {}
    with open(EGGNOG_ANNOTATIONS) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            query, description, preferred_name = cols[0], cols[7], cols[8]
            if description == "-" and preferred_name == "-":
                continue
            real_symbol = preferred_name if preferred_name != "-" else None
            real_description = description if description != "-" else None
            desc[query] = (real_symbol, real_description)
    return desc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gene_list_file")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.gene_list_file) as f:
        gene_ids = sorted({line.strip().rstrip(".") for line in f if line.strip()})

    print(f"build_gene_orthology.py: loading reference data...", file=sys.stderr)
    gene2tx = load_gene2tx()
    blast_hits = load_blast_hits()
    of_arabidopsis = load_orthofinder_arabidopsis()
    tair_desc = load_tair_descriptions()
    eggnog_desc = load_eggnog()

    n_asteraceae = n_arabidopsis = n_eggnog = n_none = n_mechanical_acronym = 0
    with open(args.out, "w") as out:
        # Acronym/Function_Description are filled for the OrthoFinder and eggNOG tiers;
        # compile_gene_summary.py queries UniProt for the rest.
        out.write("Gene_ID\tOrthologs\tOrtholog_Source\tAcronym\tFunction_Description\tBLAST_Gene_Names\n")
        for gid in gene_ids:
            txs = gene2tx.get(gid, [gid])

            ortholog_text = None
            acronym = None
            function_description = None
            source = "none"

            # 1. Asteraceae BLAST hit (strict e-value); source names the matched genus
            for tx in txs:
                hit = best_asteraceae_hit(blast_hits.get(tx, []))
                if hit:
                    ortholog_text, genus = hit
                    source = f"BLAST:{genus}"
                    break

            # 2. OrthoFinder -> Arabidopsis, described via local TAIR10 FASTA
            if ortholog_text is None:
                at_genes = []
                for tx in txs:
                    at_genes.extend(of_arabidopsis.get(tx, []))
                for at_gene in at_genes:
                    if at_gene in tair_desc:
                        symbol, description = tair_desc[at_gene]
                        ortholog_text = f"{(symbol + ': ') if symbol else ''}{description} [Arabidopsis thaliana]"
                        acronym = symbol
                        function_description = f"{description} [Arabidopsis thaliana]"
                        source = "OrthoFinder:Arabidopsis_thaliana"
                        break

            # 3. eggNOG description for the gene itself
            if ortholog_text is None:
                for tx in txs:
                    if tx in eggnog_desc:
                        symbol, description = eggnog_desc[tx]
                        parts = [p for p in (symbol, description) if p]
                        ortholog_text = ": ".join(parts) if parts else None
                        acronym = symbol
                        function_description = description
                        source = "eggNOG"
                        break

            if source.startswith("BLAST:"):
                n_asteraceae += 1
            elif source == "OrthoFinder:Arabidopsis_thaliana":
                n_arabidopsis += 1
            elif source == "eggNOG":
                n_eggnog += 1
            else:
                n_none += 1

            # No curated symbol (mainly the BLAST tier): extract one from the Orthologs text
            if acronym is None and ortholog_text is not None:
                acronym = extract_acronym(ortholog_text)
                if acronym:
                    n_mechanical_acronym += 1

            # BLAST_Gene_Names: general top hit, unrestricted, only when nothing above found
            blast_gene_name = "NA"
            if ortholog_text is None:
                for tx in txs:
                    hits = blast_hits.get(tx, [])
                    if hits:
                        blast_gene_name = hits[0][1]
                        break

            out.write(f"{gid}\t{ortholog_text or 'NA'}\t{source}\t"
                      f"{acronym or 'NA'}\t{function_description or 'NA'}\t{blast_gene_name}\n")

    print(f"build_gene_orthology.py: {len(gene_ids)} genes -- "
          f"{n_asteraceae} via Asteraceae BLAST, {n_arabidopsis} via OrthoFinder/Arabidopsis, "
          f"{n_eggnog} via eggNOG, {n_none} with no ortholog call; "
          f"{n_mechanical_acronym} otherwise-symbol-less genes got one mechanically "
          f"extracted from their Orthologs text", file=sys.stderr)


if __name__ == "__main__":
    main()
