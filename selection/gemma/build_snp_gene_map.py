#!/usr/bin/env python3
"""
Builds <var_name>/snp_to_gene_<var_name>.tsv: for each outlier SNP in
outliers_<var_name>.ann.INFO, the genes SnpEff assigned it to. Used by
annotated_manhattan.R and select_peak_labels.R to label peaks.

Uses the same gene-crediting rules as curate_intergenic_by_ld.py: genic and
upstream/downstream annotations are always kept; a SNP annotated only as
intergenic_region credits the flanking gene only if the SNP is in LD with it
(r2 >= --r2-threshold, from the ld_curation/ scan).

Gene IDs are read only from the ANN Gene_ID subfield (field 5); Feature_ID
can carry old AUGUSTUS "g<number>.t<N>" stems that would match unrelated genes.

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

sys.path.insert(0, f"{BASE_DIR}/gvtp_popgen/selection/gemma")
from curate_intergenic_by_ld import GENE_TOKEN_RE, load_gff_resolver, resolve_to_gene, load_linked_positions

BASE_DIR = "/private/groups/meyerlab/esteiger"


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
