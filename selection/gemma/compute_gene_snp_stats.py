#!/usr/bin/env python3
"""
Per gene in gene_list_<var>_unfiltered.txt: number of linked outlier SNPs and
their minimum p-value (lead SNP), for the SNP_Count/Min_P_Value columns of
gene_summary_all.tsv.

SNP-gene linkage follows curate_intergenic_by_ld.py: genic and
upstream/downstream annotations always link; a SNP annotated only as
intergenic_region links only if it passed the LD scan at --r2-threshold.
A SNP linked to two genes counts toward both.

p-value: column 11 (p_final) of unfiltered_outliers_for_annot_<var>.txt,
the value find_outliers.R compares against the Bonferroni threshold.

Usage:
  python3 compute_gene_snp_stats.py <var_name> --ld-curation-dir DIR [--r2-threshold 0.5]
Writes <var_name>/gene_snp_stats_<var_name>.tsv: gene_id, snp_count, min_p
"""
import argparse
import sys
from collections import defaultdict

sys.path.insert(0, f"{BASE_DIR}/gvtp_popgen/selection/gemma")
from curate_intergenic_by_ld import GENE_TOKEN_RE, load_gff_resolver, resolve_to_gene, load_linked_positions

BASE_DIR = "/private/groups/meyerlab/esteiger"

GEMMA_DIR = f"{BASE_DIR}/selection/gemma"


def load_scaffold_snp_map(path):
    """(gemma_chr, snp_name) -> (scaffold, vcf_pos)"""
    m = {}
    with open(path) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 4:
                continue
            m[(cols[0], cols[1])] = (cols[2], int(cols[3]))
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("var_name")
    ap.add_argument("--ld-curation-dir", required=True)
    ap.add_argument("--r2-threshold", type=float, default=0.5)
    args = ap.parse_args()
    var = args.var_name

    scaf_map = load_scaffold_snp_map(f"{GEMMA_DIR}/scaffold_snp_map.txt")

    outlier_file = f"{GEMMA_DIR}/{var}/unfiltered_outliers_for_annot_{var}.txt"
    snp_pfinal = {}
    with open(outlier_file) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 12:
                continue
            key = (cols[0], cols[1])
            if key not in scaf_map:
                continue
            scaffold, vcf_pos = scaf_map[key]
            snp_pfinal[(scaffold, vcf_pos)] = float(cols[10])

    ann_file = f"{GEMMA_DIR}/{var}/outliers_{var}.ann.INFO"
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

    gene_snps = defaultdict(set)

    for chrom, pos, ann in rows:
        subs = [s.split("|") for s in ann.split(",")]
        subs = [s for s in subs if len(s) >= 5]
        has_non_intergenic = any(s[1] != "intergenic_region" for s in subs)

        for s in subs:
            if s[1] != "intergenic_region":
                # Gene_ID subfield only; Feature_ID stems (e.g. "g1478.t2") match unrelated genes
                for tok in GENE_TOKEN_RE.findall(s[4]):
                    gid = (resolve_to_gene(tok, gene_ids_gff, tx2gene) or tok).rstrip(".")
                    gene_snps[gid].add((chrom, pos))
                continue

            if has_non_intergenic:
                continue

            for tok in GENE_TOKEN_RE.findall(s[4]):
                gid = resolve_to_gene(tok, gene_ids_gff, tx2gene)
                if gid is None:
                    continue
                gid = gid.rstrip(".")
                positions = linked.get((chrom, gid))
                if positions and pos in positions:
                    gene_snps[gid].add((chrom, pos))

    out_file = f"{GEMMA_DIR}/{var}/gene_snp_stats_{var}.tsv"
    n_written = 0
    with open(out_file, "w") as out:
        out.write("gene_id\tsnp_count\tmin_p\n")
        for gid, snps in sorted(gene_snps.items()):
            pvals = [snp_pfinal[s] for s in snps if s in snp_pfinal]
            if not pvals:
                continue
            out.write(f"{gid}\t{len(pvals)}\t{min(pvals):.6e}\n")
            n_written += 1

    print(f"compute_gene_snp_stats.py: wrote stats for {n_written} genes to {out_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
