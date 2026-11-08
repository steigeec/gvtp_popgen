#!/usr/bin/env python3
"""
Genotype-based check of the heterozygote excess at the chromosome 6 centromere
in Tranquillon Mtn (Strauss), where the "Het" PCA group of breakpoint_samples.tsv
was assigned from k-means clusters.

1. Diagnostic SNPs in REGION: |ALT frequency (DS / 2) in PCA group A - group B|
   >= MIN_DIFF, oriented to count the A-group allele.
2. Per plant (untrustworthy-coordinate plants excluded): hybrid index = mean
   A-allele dosage / 2, and heterozygosity = fraction of diagnostic SNPs with
   rounded dosage 1.
3. Genotype classes (AA > 1 - CUT, BB < CUT, else AB), exact HWE test
   (Wigginton et al. 2005), and a cross-tab against the PCA groups.

Usage: centromere_hwe_check.py
Output: OUT_TSV (per plant), OUT_FIG
"""

import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, f"{BASE_DIR}/gvtp_popgen/selection/gemma")
import temporal_null as tn

BASE_DIR = "/private/groups/meyerlab/esteiger"

VCF = f"{BASE_DIR}/variant_calling/angsd/postParalog_MAF_auto_GLimp_clean.vcf.gz"
GROUPS = f"{BASE_DIR}/inversions/breakpoints/breakpoint_samples.tsv"
FIELD = f"{BASE_DIR}/field_data.txt"
OUT_TSV = f"{BASE_DIR}/inversions/centromere/centromere_hybrid_index_strauss.tsv"
OUT_FIG = f"{BASE_DIR}/inversions/centromere/centromere_hybrid_index_strauss.jpg"
CHROM, REGION = "DEINVI_00006", (45_500_000, 48_500_000)
MIN_DIFF, CUT = 0.8, 0.25


def hwe_exact(n_aa, n_ab, n_bb):
    """Wigginton et al. (2005) exact HWE test; returns (p_two_sided, p_het_excess)"""
    n = n_aa + n_ab + n_bb; n_a = 2 * n_aa + n_ab; n_b = 2 * n - n_a
    rare = min(n_a, n_b)
    probs = np.zeros(rare + 1)
    mid = int(rare * (2 * n - rare) / (2 * n))
    if (rare - mid) % 2:
        mid += 1
    probs[mid] = 1.0
    h = mid
    while h > 1:
        hom_r = (rare - h) // 2; hom_c = n - h - hom_r
        probs[h - 2] = probs[h] * h * (h - 1) / (4 * (hom_r + 1) * (hom_c + 1)); h -= 2
    h = mid
    while h <= rare - 2:
        hom_r = (rare - h) // 2; hom_c = n - h - hom_r
        probs[h + 2] = probs[h] * 4 * hom_r * hom_c / ((h + 2) * (h + 1)); h += 2
    probs /= probs.sum()
    obs = probs[n_ab]
    return probs[probs <= obs].sum(), probs[n_ab:].sum()


def main():
    f = pd.read_csv(FIELD, sep="\t", dtype={"SampleID": str}).set_index("SampleID")
    vcf_samples = subprocess.run(["bcftools", "query", "-l", VCF], capture_output=True, text=True).stdout.split()
    st = [s for s in vcf_samples if s in f.index and f.loc[s, "Location"] == "Strauss" and s not in tn.BAD_COORD_IDS]
    g = pd.read_csv(GROUPS, sep="\t").query("block == 'INV1_45.5-49.1'")
    g["SampleID"] = g["sample"].str.replace("_rmdup.bam", "", regex=False)
    pca_group = g.set_index("SampleID").group.reindex(st)
    out = subprocess.run(["bcftools", "query", "-r", f"{CHROM}:{REGION[0]}-{REGION[1]}", "-s", ",".join(st),
                          "-f", "%POS[\t%DS]\n", VCF], capture_output=True, text=True).stdout.splitlines()
    a = np.array([l.split("\t") for l in out], dtype=float)
    pos, D = a[:, 0].astype(int), a[:, 1:]                      # SNPs x plants
    A, B = (pca_group == "A").values, (pca_group == "B").values
    fa, fb = D[:, A].mean(1) / 2, D[:, B].mean(1) / 2
    diag = np.abs(fa - fb) >= MIN_DIFF
    Dd = D[diag]; Dd = np.where((fa < fb)[diag][:, None], 2 - Dd, Dd)   # count the A-group allele
    print(f"{len(st)} Tranquillon Mtn plants; {len(pos):,} SNPs in {REGION[0]/1e6:.1f}-{REGION[1]/1e6:.1f} Mb; "
          f"{diag.sum()} diagnostic (|freq A - freq B| >= {MIN_DIFF}; PCA groups A n={A.sum()}, B n={B.sum()})")
    idx = Dd.mean(0) / 2
    het = (np.rint(Dd) == 1).mean(0)
    cls = np.where(idx > 1 - CUT, "AA", np.where(idx < CUT, "BB", "AB"))
    res = pd.DataFrame(dict(SampleID=st, pca_group=pca_group.values, hybrid_index=idx, diag_het_frac=het, geno_class=cls,
                            Year=f.loc[st, "Year"].values, cov=pd.to_numeric(f.loc[st, "Cov"], errors="coerce").values))
    res.to_csv(OUT_TSV, sep="\t", index=False)

    n = res.geno_class.value_counts().reindex(["AA", "AB", "BB"]).fillna(0).astype(int)
    p_a = (2 * n.AA + n.AB) / (2 * n.sum())
    exp = n.sum() * np.array([p_a ** 2, 2 * p_a * (1 - p_a), (1 - p_a) ** 2])
    p2, pex = hwe_exact(n.AA, n.AB, n.BB)
    print(f"\ngenotype classes (hybrid index): AA {n.AA}, AB {n.AB}, BB {n.BB}; freq(A) = {p_a:.3f}")
    print(f"expected under HWE: AA {exp[0]:.1f}, AB {exp[1]:.1f}, BB {exp[2]:.1f}; observed/expected het = {n.AB / exp[1]:.2f}; "
          f"F = {1 - n.AB / exp[1]:+.3f}; exact HWE p = {p2:.3g} (het-excess one-sided p = {pex:.3g})")
    print("\nPCA group vs genotype class (plants with a PCA group):")
    print(pd.crosstab(res.pca_group.fillna("not grouped"), res.geno_class).to_string())
    print("\nmedian diagnostic heterozygosity by genotype class:", res.groupby("geno_class").diag_het_frac.median().round(2).to_dict())
    for y, x in res.groupby("Year"):
        m = x.geno_class.value_counts().reindex(["AA", "AB", "BB"]).fillna(0).astype(int)
        q = (2 * m.AA + m.AB) / (2 * m.sum()); e = 2 * q * (1 - q) * m.sum()
        print(f"  {y}: AA {m.AA}, AB {m.AB}, BB {m.BB}; freq(A) {q:.2f}; obs/exp het {m.AB / e:.2f}; HWE p {hwe_exact(m.AA, m.AB, m.BB)[0]:.3g}")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    a1.hist(idx, bins=50, color="#7f7f7f")
    for c in (CUT, 1 - CUT):
        a1.axvline(c, ls="--", color="black", lw=1)
    a1.set_xlabel("Hybrid index (A-allele fraction, diagnostic SNPs)", fontsize=13); a1.set_ylabel("Plants", fontsize=13)
    col = {"A": "#4d4d4d", "Het": "#d55e00", "B": "#56b4e9"}
    for grp, x in res.groupby(res.pca_group.fillna("not grouped")):
        a2.scatter(x.hybrid_index, x.diag_het_frac, s=12, alpha=0.7, color=col.get(grp, "#cccccc"), label=f"PCA {grp} (n={len(x)})")
    a2.set_xlabel("Hybrid index", fontsize=13); a2.set_ylabel("Heterozygosity at diagnostic SNPs", fontsize=13)
    a2.legend(frameon=False, fontsize=10)
    for ax in (a1, a2):
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT_FIG, dpi=170, bbox_inches="tight")
    print(f"Saved: {OUT_FIG}")


if __name__ == "__main__":
    main()
