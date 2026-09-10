#!/usr/bin/env python3
"""
Computes, per gene in gene_list_<var>_unfiltered.txt, how many outlier SNPs
actually linked to it and their minimum Bonferroni-tested p-value -- for
gene_summary_all.tsv's SNP_Count/Min_P_Value columns.

Minimum (not mean) because linked SNPs are often not independent -- the
LD-tagged intergenic ones in particular, by construction, come from the
same LD block as one real signal -- so averaging would dilute a gene with
one strong hit plus several correlated marginal ones. Min matches the
standard GWAS convention of characterizing a locus by its lead SNP.

"Linked" uses the exact same rule that built the gene list in the first
place (imported from curate_intergenic_by_ld.py, not re-derived): genic and
bounded upstream/downstream_gene_variant annotations link unconditionally;
a SNP whose only annotation is intergenic_region links only if it passed
the LD scan at --r2-threshold. A SNP touching two genes at once (e.g.
upstream of one, downstream of another) counts toward both.

p_final is GEMMA's own genomic-control-corrected p-value (raw p_lrt if no
correction was needed), column 11 (0-indexed 10) of
unfiltered_outliers_for_annot_<var>.txt -- the same value find_outliers.R
compares against the Bonferroni threshold.

Usage:
  python3 compute_gene_snp_stats.py <var_name> --ld-curation-dir DIR [--r2-threshold 0.5]
Writes <var_name>/gene_snp_stats_<var_name>.tsv: gene_id, snp_count, min_p
"""
import argparse
import sys
from collections import defaultdict

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
from curate_intergenic_by_ld import GENE_TOKEN_RE, load_gff_resolver, resolve_to_gene, load_linked_positions

GEMMA_DIR = "/private/groups/meyerlab/esteiger/selection/gemma"


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
                # s[4] (Gene_ID) only -- see curate_intergenic_by_ld.py's
                # matching fix for why scanning the whole entry (s joined
                # back with "|") picks up Feature_ID transcript stems
                # (e.g. "g1478.t2") as spurious extra gene calls.
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
