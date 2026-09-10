#!/usr/bin/env python3
"""
Rebuilds <var_name>/snp_to_gene_<var_name>.tsv: for every outlier SNP in
outliers_<var_name>.ann.INFO, the set of genes SnpEff actually assigned it
to -- consumed by annotated_manhattan.R / select_peak_labels.R to label a
Manhattan peak with the gene(s) tied to its own retained SNP(s), instead of
a coordinate-range overlap (which misses genes reached only through
SnpEff's upstream/downstream/intergenic flanking window -- the majority
case here, since most outlier SNPs are intergenic).

Recreates a script by this name referenced in those two R scripts' headers
that no longer exists on disk (the on-disk snp_to_gene_<var>.tsv files were
last written 2026-07-14 with no surviving generator). Reuses the same
GENE_TOKEN_RE / load_gff_resolver / resolve_to_gene / load_linked_positions
helpers curate_intergenic_by_ld.py and compute_gene_snp_stats.py already
use, so a SNP's gene set here matches what those scripts would credit it
with: genic and bounded upstream/downstream_gene_variant annotations are
kept unconditionally; a SNP whose ONLY annotation is intergenic_region
(unbounded nearest-flanking-gene distance) only credits that gene if the
SNP's own position is in the gene's LD-linked position set (r2 >=
--r2-threshold), per the precomputed ld_curation/ scan.

Gene IDs are read from ONLY the Gene_ID subfield (pipe-field 5, i.e.
parts[4]) of each ANN entry -- never the whole entry text. AUGUSTUS
transcript IDs (Feature_ID, pipe-field 7) commonly retain a raw
"g<number>.t<N>" stem even for genes since renamed to "IDmodified-gene-N",
and scanning the whole entry lets that stem get matched as a spurious
extra gene call that can collide with a real, unrelated gene elsewhere in
the genome that kept that same "g<number>" id (confirmed: this is exactly
how SCPL13 and UGT73E1 -- on different scaffolds -- previously ended up
labeled on the same BIO15 peak).

Usage:
  python3 build_snp_gene_map.py <var_name> \\
      --ld-curation-dir /path/to/ld_curation [--r2-threshold 0.5]
Writes <var_name>/snp_to_gene_<var_name>.tsv (Scaffold, Position, Gene_IDs
-- semicolon-joined, sorted).
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
from curate_intergenic_by_ld import GENE_TOKEN_RE, load_gff_resolver, resolve_to_gene, load_linked_positions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("var_name")
    ap.add_argument("--ld-curation-dir", required=True)
    ap.add_argument("--r2-threshold", type=float, default=0.5)
    args = ap.parse_args()
    var = args.var_name

    ann_file = os.path.join(var, f"outliers_{var}.ann.INFO")
    rows = []
    with open(ann_file) as f:
        header = f.readline()
        if not header.startswith("CHROM"):
            f.seek(0)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            rows.append((cols[0], int(cols[1]), cols[4]))

    scaffolds = sorted({r[0] for r in rows})
    linked = load_linked_positions(args.ld_curation_dir, scaffolds, args.r2_threshold)
    gene_ids_gff, tx2gene = load_gff_resolver()

    snp_genes = defaultdict(set)  # (scaffold, position) -> {gene_id, ...}

    for chrom, pos, ann in rows:
        subs = [s.split("|") for s in ann.split(",")]
        subs = [s for s in subs if len(s) >= 5]
        has_non_intergenic = any(s[1] != "intergenic_region" for s in subs)

        for s in subs:
            if s[1] != "intergenic_region":
                for tok in GENE_TOKEN_RE.findall(s[4]):
                    gid = (resolve_to_gene(tok, gene_ids_gff, tx2gene) or tok).rstrip(".")
                    snp_genes[(chrom, pos)].add(gid)
                continue

            if has_non_intergenic:
                # redundant with a bounded annotation already credited above
                continue

            for tok in GENE_TOKEN_RE.findall(s[4]):
                gid = resolve_to_gene(tok, gene_ids_gff, tx2gene)
                if gid is None:
                    continue
                positions = linked.get((chrom, gid))
                if positions and pos in positions:
                    snp_genes[(chrom, pos)].add(gid)

    out_file = os.path.join(var, f"snp_to_gene_{var}.tsv")
    with open(out_file, "w") as out:
        out.write("Scaffold\tPosition\tGene_IDs\n")
        for (scaffold, pos), gids in sorted(snp_genes.items()):
            if not gids:
                continue
            out.write(f"{scaffold}\t{pos}\t{';'.join(sorted(gids))}\n")

    n_rows = sum(1 for gids in snp_genes.values() if gids)
    print(f"build_snp_gene_map.py: wrote {n_rows} SNP->gene rows to {out_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
