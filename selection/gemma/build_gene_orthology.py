#!/usr/bin/env python3
"""
Builds the Orthologs / BLAST_Gene_Names columns for compile_gene_summary.py,
replacing the old approach of using raw BLASTP top-hit as "Orthologs"
(that's just sequence similarity, not orthology -- see the cylicin-2-like /
CLPTM1 cross-kingdom false hits this produced during the gene_summary_all.tsv
literature pass).

For each gene, in priority order:
  1. Best BLASTP hit from any Asteraceae species (ASTERACEAE_GENERA
     allowlist, matched on exact genus -- not substring, since substring
     matching on "aster" alone pulls in false positives like Cotoneaster
     [Rosaceae], Verticillium dahliae [a fungus], and Acanthaster [a
     starfish]) with evalue <= ASTERACEAE_EVALUE_MAX. Restricting to the
     actual plant family our species belongs to, rather than either an
     unrestricted whole-database search (which produced spurious cross-
     kingdom hits like cylicin-2-like [mammalian sperm-specific] and a
     literal Plasmodium falciparum locus tag) or an arbitrarily narrow
     single-species/OrthoFinder-species-only search (which misses strong,
     genuine hits like Cynara cardunculus/artichoke -- second only to
     sunflower in hit volume in this database, but absent from the
     OrthoFinder comparison, which included non-Asteraceae references like
     Arabidopsis and tomato for annotation-curation-quality reasons, not
     taxonomic proximity). A very strict e-value plus a real-family
     restriction gives orthology-level confidence despite BLAST itself
     doing no orthology inference -- the hit description text is used
     directly (it's already human-readable). No species preference beyond
     the e-value ranking itself: whichever qualifying Asteraceae hit scores
     best wins, which in practice is usually (not always) sunflower.
  2. Else, the gene's OrthoFinder (Results_Jul18) ortholog call to
     Arabidopsis_thaliana, resolved via the local TAIR10 peptide FASTA
     header (which embeds a curated gene_symbol + description -- no web
     lookup needed). OrthoFinder gives real phylogenetically-inferred
     orthology, unlike raw BLAST.
  3. Else, the gene's own eggNOG-mapper Description/Preferred_name
     (Di2_viridiplantae run -- taxonomy-restricted to plants, avoiding the
     cross-kingdom false-hit problem raw whole-database BLAST has).
  4. Else, Orthologs is empty for this gene.

BLAST_Gene_Names is the general (unrestricted-species, no e-value filter)
BLASTP top hit -- i.e. what "Orthologs" used to mean -- and is only written
when steps 1-3 all come up empty, per the "only fall back to BLAST if
Orthofinder/eggNOG found nothing" rule.

Usage:
  python3 build_gene_orthology.py <gene_list_file> --out orthologs_<var>.txt
Output columns: Gene_ID, Orthologs, Ortholog_Source, BLAST_Gene_Names
"""
import argparse
import re
import sys
from collections import defaultdict

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
from curate_intergenic_by_ld import load_gff_resolver

BLASTP_OUT = "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.blastp.out"
ORTHOFINDER_DEINANDRA_TSV = ("/private/groups/meyerlab/esteiger/genome/annotation/orthofinder_results/"
                             "Results_Jul18/Orthologues/Deinandra.tsv")
TAIR_PEP_FASTA = "/private/groups/meyerlab/esteiger/genome/annotation/Arabidopsis_thaliana.TAIR10.pep.all.fa"
EGGNOG_ANNOTATIONS = "/private/groups/meyerlab/esteiger/genome/annotation/gvtp_eggnog_Di2_viridiplantae.emapper.annotations"

ASTERACEAE_EVALUE_MAX = 1e-50

# Confirmed by grepping every distinct organism bracket in BLASTP_OUT and
# checking each against known Asteraceae taxonomy -- NOT a substring match
# (that pulls in Cotoneaster [Rosaceae], Verticillium dahliae [a fungus],
# Circaeaster/Guyanagaster/Aztecaster/Acanthaster [unrelated genera, a
# starfish in the last case] purely because they contain "aster"/"dahlia").
# Deinandra's own tarweed-tribe (Madiinae) relatives (Madia, Layia,
# Hemizonia, Holocarpha, Wyethia, Arnica, Encelia) were checked for and are
# not present in this database.
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
    ASTERACEAE_EVALUE_MAX, or None. Matches on the exact genus (first word
    after '['), not a substring of the organism name -- see
    ASTERACEAE_GENERA's docstring for why."""
    for evalue, stitle in tx_hits:
        if evalue > ASTERACEAE_EVALUE_MAX:
            break  # sorted by evalue -- nothing later can qualify either
        m = ORGANISM_RE.search(stitle)
        if m and m.group(1) in ASTERACEAE_GENERA:
            return stitle, m.group(1)
    return None


ACRONYM_TOKEN_RE = re.compile(r'\b([A-Z][A-Za-z0-9]{1,9})\b')
ISOFORM_TOKEN_RE = re.compile(r'^X\d+$')


def extract_acronym(orthologs_text):
    """Best-effort mechanical extraction of an embedded gene symbol/locus id
    from the Orthologs hit text itself (e.g. "...mitochondrial acidic
    protein MAM33 [Helianthus annuus]" -> "MAM33"). Deliberately Orthologs-
    only, never Functional_Description -- that field can come from a
    UniProt match to an unrelated organism (confirmed: one row's
    Functional_Description named a maize-specific symbol, "ZmCASPL5A1",
    while its actual verified Orthologs hit was to Cynara cardunculus/
    artichoke -- extracting from it would misattribute a different
    species' gene symbol to this one). Requires the candidate token to
    contain a digit, which reliably separates real symbols (MAM33, MTERF4,
    NST1, GLO1, ALG2, At1g10780) from generic domain/family/superfamily
    names that also happen to be all-caps (GPI, ABC, MYM, SSR, TRAFAC,
    CASP) -- checked empirically against this dataset before adopting the
    rule. Returns None rather than guess when nothing qualifies."""
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
            # TAIR/Ensembl repeats the AT-number as gene_symbol when no real
            # symbol is assigned -- that's not an acronym, treat as absent.
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
        # Acronym/Function_Description are populated directly here (no live
        # lookup needed) for the OrthoFinder/eggNOG tiers, since those
        # sources already carry a curated symbol + description -- compile_
        # gene_summary.py only needs to fall back to a UniProt query for
        # genes where these are still NA (BLAST:<genus> tier, or no
        # ortholog call at all).
        out.write("Gene_ID\tOrthologs\tOrtholog_Source\tAcronym\tFunction_Description\tBLAST_Gene_Names\n")
        for gid in gene_ids:
            txs = gene2tx.get(gid, [gid])

            ortholog_text = None
            acronym = None
            function_description = None
            source = "none"

            # 1. Asteraceae-specific, strict-evalue BLAST hit (no curated
            # symbol/description available -- just the raw hit description).
            # Source label names the actual matched genus for transparency.
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

            # Mechanical fallback: if no tier supplied a curated symbol
            # (mainly the BLAST tier, which never carries one), try to pull
            # one directly out of the Orthologs hit text itself -- see
            # extract_acronym's docstring for why it's Orthologs-only.
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
