#!/usr/bin/env python3
"""
Temporal test for every fog outlier peak (the chr 6 peak is shown in
DEINVI_0006_stacked.py panel J).

Peaks: retained (40000_3) fog outlier SNPs, grouped into peaks of SNPs on the
same scaffold <= 40 kb apart (as in fog_deinvi6_figure.R).
Per peak and Strauss individual: favorable-allele dosage averaged over the
peak's SNPs. Favorable = allele whose coverage-adjusted dosage increases with
long-term fog (same polarization as fog_deinvi6_figure.R).
Statistic (same as panel J): favorable-allele frequency in FOGGY_YEAR minus
the pooled other years.
  - Random-SNP null: random SNPs (one per random 2 kb region off chr 6, same
    regions as DEINVI_0006_stacked.py) polarized like the outliers
    (coverage-adjusted regression on long-term fog across all samples, as in
    germ_year_polygenic.R), matched on Strauss favorable-allele frequency.
  - Centromere-matched null: as the random-SNP null, but the random SNPs are
    drawn at a similar distance from their own scaffold's centromere as the
    peak is from its centromere (temporal_null.py: d +- max(1 Mb, 25% of d),
    scaffolds with a confident centromere call, chr 6 and +-500 kb around
    fog peaks excluded). NaN for peaks on scaffolds with no call.
  - Within-stand permutation: years shuffled within OHANA_order stands.
  - Benjamini-Hochberg FDR across peaks for both p-values.

Usage: fog_peaks_temporal.py
Output: OUT (one row per peak)
"""

import numpy as np
import pandas as pd

import temporal_null as tn

BASE_DIR = "/private/groups/meyerlab/esteiger"

GEMMA = f"{BASE_DIR}/selection/gemma"
OUTLIERS = f"{GEMMA}/fog/outliers_for_annot_40000_3_fog.txt"
DOSAGE = f"{GEMMA}/fog/dosage_fog.txt"
SCAF_MAP = f"{GEMMA}/scaffold_snp_map.txt"
ENV = f"{BASE_DIR}/env_var.txt"
FIELD = f"{BASE_DIR}/field_data.txt"
VCF_DS = f"{BASE_DIR}/variant_calling/angsd/postParalog_MAF_auto_GLimp_clean.vcf.gz"
FAI = f"{BASE_DIR}/inversions/filtered_scaffolds.fai"
NULL_CACHE = f"{GEMMA}/fog/temporal_null_random_snps_allsamples.tsv.gz"
N_NULL_REGIONS = 8000
OUT = f"{GEMMA}/fog/fog_peaks_temporal.tsv"
FOGGY_YEAR, PEAK_GAP, FREQ_TOL, N_PERM = 2021, 40_000, 0.04, 5000
rng = np.random.default_rng(3)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def shift(d, years):
    f = years == FOGGY_YEAR
    return d[f].mean() / 2 - d[~f].mean() / 2


def random_snp_dosages(samples):
    """Imputed DS dosages (all samples) for SNPs in random 2 kb regions off chr 6; cached."""
    import os, subprocess, tempfile
    if os.path.exists(NULL_CACHE):
        return pd.read_csv(NULL_CACHE, sep="\t", index_col=[0, 1])
    fai = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"])
    fai = fai[fai.chrom != "DEINVI_00006"]
    r = np.random.default_rng(7)
    chroms = r.choice(fai.chrom.values, N_NULL_REGIONS, p=fai.len.values / fai.len.sum())
    lens = fai.set_index("chrom")["len"]
    starts = [int(r.integers(1, lens[c] - 2000)) for c in chroms]
    reg = pd.DataFrame({"chrom": chroms, "start": starts, "end": [x + 2000 for x in starts]}).sort_values(["chrom", "start"])
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        reg.to_csv(fh, sep="\t", header=False, index=False); regfile = fh.name
    out = subprocess.run(["bcftools", "query", "-R", regfile, "-s", ",".join(samples),
                          "-f", "%CHROM\t%POS[\t%DS]\n", VCF_DS], capture_output=True, text=True).stdout.splitlines()
    os.remove(regfile)
    d = pd.DataFrame([l.split("\t") for l in out]).drop_duplicates([0, 1]).set_index([0, 1]).astype(float)
    d.index.names = ["chrom", "pos"]; d.columns = samples
    d.to_csv(NULL_CACHE, sep="\t", compression="gzip")
    return d


def main():
    o = pd.read_csv(OUTLIERS, sep=" ", header=None).iloc[:, :3]
    o.columns = ["chr", "snp", "pos"]
    o = o.drop_duplicates("snp").sort_values(["chr", "pos"])
    o["peak"] = ((o.chr != o.chr.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    smap = pd.read_csv(SCAF_MAP, sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    num2scaf = dict(zip(smap[0], smap[2]))

    ds = pd.read_csv(DOSAGE, sep="\t", index_col=0); ds.columns = ds.columns.astype(str)
    env = pd.read_csv(ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    common = [c for c in ds.columns if c in env.index and c not in tn.BAD_COORD_IDS]   # drop untrustworthy-coordinate IDs
    fog = env.loc[common, "fog"].astype(float).values
    cov = pd.to_numeric(env.loc[common, "Cov"], errors="coerce").values
    ok = ~np.isnan(cov)
    X = np.column_stack([np.ones(ok.sum()), cov[ok]]); H = X @ np.linalg.solve(X.T @ X, X.T)
    fr = fog[ok] - H @ fog[ok]

    strauss = [s for s in common if env.loc[s, "Location"] == "Strauss"]
    years = env.loc[strauss, "Year"].astype(int).values
    field = pd.read_csv(FIELD, sep=None, engine="python", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    stand = field["OHANA_order"].reindex(strauss).values
    groups = [np.where(stand == st)[0] for st in pd.unique(stand[~pd.isna(stand)])]

    # random SNPs, thinned to one per 2 kb region, polarized like the outliers
    null = random_snp_dosages(common)
    null.columns = null.columns.astype(str); null = null[common]
    region = pd.Series([f"{c}:{int(p) // 2000}" for c, p in null.index])
    N = null.values.T[ok]                                  # samples x SNPs
    beta = (N - H @ N).T @ fr
    P = np.where(beta[None, :] > 0, null.values.T, 2 - null.values.T)   # all samples x SNPs, favorable dose
    si = [common.index(s) for s in strauss]
    PS = P[si]
    nf = pd.Series(PS.mean(0) / 2)                         # Strauss favorable-allele frequency
    f21 = years == FOGGY_YEAR
    null_shift = PS[f21].mean(0) / 2 - PS[~f21].mean(0) / 2
    poly = (PS.mean(0) > 0) & (PS.mean(0) < 2)
    nf, null_shift, region = nf[poly].reset_index(drop=True), null_shift[poly], region[poly].reset_index(drop=True)
    print(f"random-SNP null: {len(null_shift)} SNPs polarized by long-term fog; "
          f"mean 2021 shift {null_shift.mean():+.4f} (sd {null_shift.std():.4f})")

    # centromere-matched null: one distance band per peak with a centromere call
    cen = tn.centromeres()
    pk_tab = o.groupby("peak").agg(chr=("chr", "first"), start=("pos", "min"), end=("pos", "max"))
    pk_tab["chrom"] = pk_tab.chr.map(num2scaf)
    pk_tab["dist"] = [tn.cen_dist(c, (a + b) / 2, cen) for c, a, b in zip(pk_tab.chrom, pk_tab.start, pk_tab.end)]
    pk_tab["band_id"] = pk_tab.chrom + ":" + pk_tab.start.astype(str)
    bands = {r.band_id: tn.cen_band(r.dist) for r in pk_tab.itertuples() if r.dist == r.dist}
    cnull = tn.cen_null_dosages(common, bands, pk_tab[["chrom", "start", "end"]])
    cP = tn.polarize(cnull)[strauss]
    c_nf = cP.mean(axis=1).values / 2
    c_shift = cP.loc[:, years == FOGGY_YEAR].mean(axis=1).values / 2 - cP.loc[:, years != FOGGY_YEAR].mean(axis=1).values / 2
    c_band = cnull.index.get_level_values("band").values
    c_region = np.array([f"{c}:{int(p) // 2000}" for c, p, _ in cnull.index])
    c_poly = (c_nf > 0) & (c_nf < 1)
    print(f"centromere-matched null: {len(bands)} peak bands, {cnull.shape[0]} SNPs")

    rows = []
    for pk, g in o.groupby("peak"):
        snps = [s for s in g.snp if s in ds.index]
        if not snps:
            continue
        sign = []
        for snp in snps:
            d = ds.loc[snp, common].values.astype(float)[ok]
            sign.append(1 if ((d - H @ d) @ fr) > 0 else -1)
        D = ds.loc[snps, strauss].values.astype(float)
        D = np.where(np.array(sign)[:, None] > 0, D, 2 - D)
        dose = np.nanmean(D, axis=0)
        f_fav = np.nanmean(dose) / 2
        obs = shift(dose, years)
        # frequency-match, then thin to one SNP per random 2 kb region (as in DEINVI_0006_stacked.py)
        m = np.where(np.abs(nf.values - f_fav) <= FREQ_TOL)[0]
        m = pd.Series(m).groupby(region.values[m]).apply(lambda x: x.sample(1, random_state=1)).values
        nm = null_shift[m]
        p_null = (np.sum(nm >= obs) + 1) / (len(nm) + 1)
        info = pk_tab.loc[pk]
        if info.band_id in bands:
            cm = np.where((c_band == info.band_id) & c_poly & (np.abs(c_nf - f_fav) <= FREQ_TOL))[0]
            cm = pd.Series(cm).groupby(c_region[cm]).apply(lambda x: x.sample(1, random_state=1)).values
            p_cen = (np.sum(c_shift[cm] >= obs) + 1) / (len(cm) + 1); n_cen = len(cm)
        else:
            p_cen, n_cen = np.nan, 0
        perm = np.empty(N_PERM)
        for k in range(N_PERM):
            y = years.copy()
            for gg in groups:
                y[gg] = rng.permutation(y[gg])
            perm[k] = shift(dose, y)
        p_stand = (np.sum(perm >= obs) + 1) / (N_PERM + 1)
        rows.append(dict(scaffold=num2scaf.get(g.chr.iloc[0], g.chr.iloc[0]), start=g.pos.min(), end=g.pos.max(),
                         n_snps=len(snps), fav_freq=round(f_fav, 3), shift_2021=obs, n_null=len(nm),
                         p_random_snp=p_null, p_within_stand=p_stand,
                         cen_dist_Mb=round(info.dist / 1e6, 2) if info.dist == info.dist else np.nan,
                         n_null_cen=n_cen, p_cen_matched=p_cen))
    r = pd.DataFrame(rows)
    r["q_random_snp"] = bh(r.p_random_snp); r["q_within_stand"] = bh(r.p_within_stand)
    has = r.p_cen_matched.notna()
    r.loc[has, "q_cen_matched"] = bh(r.loc[has, "p_cen_matched"])
    r = r.sort_values("p_random_snp")
    r.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(f"{len(r)} fog peaks; Strauss individuals {len(strauss)}; random SNPs in null pool {len(null_shift)}")
    print(r.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"\nshift > 0 (favorable allele more common in {FOGGY_YEAR}) in {np.sum(r.shift_2021 > 0)} of {len(r)} peaks; "
          f"mean shift {r.shift_2021.mean():+.4f}")


if __name__ == "__main__":
    main()
