#!/usr/bin/env python3
"""
Combine per-SNP GEMMA p-values into one ACAT p-value per gene, as the ranking
input for GSEA_GO.R.

Steps: map GEMMA chr ids to scaffolds (scaffold_snp_map.txt); apply genomic
control as in find_outliers.R (correct if lambda > 1.05); join SNPs to genes
via snp_gene_map.tsv (build_snp_gene_map.sh); drop SNPs in the putative
inversion (DEINVI_00006:43-50 Mb); combine each gene's SNPs by equal-weight
ACAT (Liu & Xie 2020); residualize -log10(ACAT p) on log10(n_snps) by LOWESS
and write it as snp_count_adjusted_score, the column GSEA_GO.R ranks on.

Output: gene_id, acat_p, n_snps, snp_count_adjusted_score
Usage: python3 build_gene_gwas_scores.py <pheno_number> <var_name> [output.tsv]
  e.g. python3 build_gene_gwas_scores.py 6 fog
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.nonparametric.smoothers_lowess import lowess

BASE_DIR = "/private/groups/meyerlab/esteiger"

GEMMA_DIR = f"{BASE_DIR}/selection/gemma"
SCAFFOLD_MAP = f"{GEMMA_DIR}/scaffold_snp_map.txt"
SNP_GENE_MAP = f"{BASE_DIR}/genome/annotation/snp_gene_map.tsv"

LAMBDA_THRESHOLD = 1.05
INV_SCAFFOLD = "DEINVI_00006"
INV_START = 43_000_000
INV_END = 50_000_000

MIN_P = 1e-300     # clip to avoid literal p=0 (division by zero in ACAT)
SMALL_P = 1e-15    # below this, use the large-value asymptotic ACAT term


def acat_combine(pvals):
    """Equal-weight Cauchy combination test. pvals: 1D numpy array, one gene's SNPs."""
    p = np.clip(pvals, MIN_P, 1 - MIN_P)
    n = len(p)
    small = p < SMALL_P
    terms = np.empty(n)
    terms[~small] = np.tan((0.5 - p[~small]) * np.pi)
    # tan((0.5-p)*pi) ~ 1/(pi*p) as p -> 0; avoids overflow
    terms[small] = 1.0 / (np.pi * p[small])
    cct_stat = terms.mean()
    if cct_stat > 1e15:
        pval = 1.0 / (cct_stat * np.pi)
    else:
        pval = 0.5 - np.arctan(cct_stat) / np.pi
    return pval


def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: build_gene_gwas_scores.py <pheno_number> <var_name> [output.tsv]")
    pheno_num = sys.argv[1]
    var_name = sys.argv[2]
    out_file = sys.argv[3] if len(sys.argv) > 3 else \
        f"{BASE_DIR}/selection/gemma/{var_name}/{var_name}_gene_gwas_scores.tsv"

    # --- chr -> scaffold crosswalk ---
    dic = pd.read_csv(SCAFFOLD_MAP, sep=r"\s+", header=None, names=["chr", "snp", "scaffold", "pos"])
    chr2scaf = dic.drop_duplicates(subset=["chr"]).set_index("chr")["scaffold"].to_dict()

    # --- load GEMMA results, apply genomic control (same logic as find_outliers.R) ---
    gemma = pd.read_csv(f"{GEMMA_DIR}/output/pheno{pheno_num}_lmm_results.assoc.txt", sep="\t")
    print(f"Loaded {len(gemma)} SNPs from pheno{pheno_num}_lmm_results.assoc.txt ({var_name})")

    p_lrt = gemma["p_lrt"].clip(lower=MIN_P, upper=1 - MIN_P).to_numpy()
    chisq_obs = chi2.isf(p_lrt, df=1)
    lam = np.median(chisq_obs) / chi2.ppf(0.5, df=1)
    print(f"Lambda: {lam:.3f}")
    if lam > LAMBDA_THRESHOLD:
        p_final = chi2.sf(chisq_obs / lam, df=1)
        print("Applied genomic control")
    else:
        p_final = p_lrt
        print("No genomic control needed")
    gemma["p_final"] = p_final

    gemma["scaffold"] = gemma["chr"].map(chr2scaf)
    n_unmapped = gemma["scaffold"].isna().sum()
    if n_unmapped:
        print(f"WARNING: {n_unmapped} SNPs had no scaffold_snp_map.txt entry for their chr — dropped")
    gemma = gemma.dropna(subset=["scaffold"])

    # --- join to SNP->gene map ---
    snp_gene = pd.read_csv(SNP_GENE_MAP, sep="\t", header=None, names=["scaffold", "ps", "gene_id"])
    merged = gemma[["scaffold", "ps", "p_final"]].merge(snp_gene, on=["scaffold", "ps"], how="inner")
    print(f"{len(merged)} (SNP, gene) rows after joining to the SNP->gene map "
          f"({merged[['scaffold', 'ps']].drop_duplicates().shape[0]} distinct SNPs, "
          f"{merged['gene_id'].nunique()} distinct genes)")

    # --- drop SNPs inside the putative inversion ---
    in_inversion = (merged["scaffold"] == INV_SCAFFOLD) & \
                   (merged["ps"] >= INV_START) & (merged["ps"] <= INV_END)
    n_inv = in_inversion.sum()
    merged = merged[~in_inversion]
    print(f"Dropped {n_inv} (SNP, gene) rows inside the putative inversion "
          f"({INV_SCAFFOLD}:{INV_START}-{INV_END})")

    # --- per-gene ACAT combination ---
    records = []
    for gene_id, sub in merged.groupby("gene_id"):
        pvals = sub["p_final"].to_numpy()
        records.append((gene_id, acat_combine(pvals), len(pvals)))
    result = pd.DataFrame(records, columns=["gene_id", "acat_p", "n_snps"]).sort_values("acat_p")

    # --- SNP-count bias correction ---
    neglog_p = -np.log10(np.clip(result["acat_p"].to_numpy(), 1e-300, 1))
    log_n_snps = np.log10(result["n_snps"].to_numpy())
    before_rho = pd.Series(neglog_p).corr(pd.Series(log_n_snps), method="spearman")
    smoothed = lowess(neglog_p, log_n_snps, frac=0.3, return_sorted=False)
    result["snp_count_adjusted_score"] = neglog_p - smoothed
    after_rho = pd.Series(result["snp_count_adjusted_score"]).corr(pd.Series(log_n_snps), method="spearman")
    print(f"SNP-count bias check: Spearman rho(log10 n_snps, -log10 p) "
          f"{before_rho:.3f} -> {after_rho:.3f} after LOWESS residualization")

    result.to_csv(out_file, sep="\t", index=False)
    print(f"Wrote {len(result)} gene-level ACAT p-values to {out_file}")
    print(f"Min ACAT p-value: {result['acat_p'].min():.3g} "
          f"| median n_snps/gene: {result['n_snps'].median():.0f}")


if __name__ == "__main__":
    main()
