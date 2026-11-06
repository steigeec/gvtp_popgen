"""
Shared pieces of the temporal (germination-year) tests' random-SNP nulls,
used by fog_peaks_temporal.py and DEINVI_0006_stacked.py (panel J).

- polarize(): favorable-allele dosage, exactly as for the outlier SNPs
  (allele whose coverage-adjusted dosage increases with long-term fog
  across all samples; germ_year_polygenic.R / fog_deinvi6_figure.R).
- cen_dist() / cen_band(): distance of a position from its scaffold's
  centromere (0 inside the call, else bp to the nearer edge; NaN if the
  scaffold has no confident call -- inversions/centromere/centromeres.tsv,
  from call_centromeres.py), and the band of distances a peak's null SNPs
  are drawn from: d +- max(CEN_BAND_MIN, CEN_BAND_FRAC * d), floored at 0.
- cen_null_dosages(): imputed DS dosages for SNPs in random 2 kb regions
  drawn uniformly within each requested distance band, on scaffolds with a
  confident centromere call, excluding EXCLUDE_CHROM and anything within
  EXCLUDE_FLANK of a fog outlier peak; cached in CEN_NULL_CACHE together
  with the band each region was drawn for.
"""

import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

VCF_DS = "/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog_MAF_auto_GLimp_clean.vcf.gz"
ENV = "/private/groups/meyerlab/esteiger/env_var.txt"
CEN = "/private/groups/meyerlab/esteiger/inversions/centromere/centromeres.tsv"
CEN_NULL_CACHE = "/private/groups/meyerlab/esteiger/selection/gemma/fog/temporal_null_cendist_snps.tsv.gz"
CEN_BAND_MIN, CEN_BAND_FRAC = 1_000_000, 0.25
CEN_REGIONS_PER_BAND, REGION_BP = 1500, 2000
EXCLUDE_CHROM, EXCLUDE_FLANK = "DEINVI_00006", 500_000
# 13 individuals with untrustworthy coordinates (Location "Strauss" but genetic
# cluster Point_Conception / VAFB_Lion_Head; see gvtp_popgen/wind/*) -- excluded
# from polarization and from every test
BAD_COORD_IDS = {"4548", "4549", "4550", "4551", "4552", "4553", "4554",
                 "4555", "4556", "4557", "4559", "4561", "4562"}


def centromeres():
    return pd.read_csv(CEN, sep="\t").query("confident").set_index("chr")


def cen_dist(chrom, pos, cen=None):
    cen = centromeres() if cen is None else cen
    if chrom not in cen.index:
        return np.nan
    c = cen.loc[chrom]
    return float(max(c.cen_start - pos, pos - c.cen_end, 0))


def cen_band(d):
    w = max(CEN_BAND_MIN, CEN_BAND_FRAC * d)
    return max(0.0, d - w), d + w


def polarize(D):
    """D: SNP x sample DataFrame of ALT dosages (all samples). Returns favorable dosages
    for the samples with coverage and fog values (same orientation rule as the outliers)."""
    env = pd.read_csv(ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    env["Cov"] = pd.to_numeric(env["Cov"], errors="coerce")
    env = env.dropna(subset=["Cov", "fog"])
    env = env[~env.index.isin(BAD_COORD_IDS)]
    D = D.copy(); D.columns = D.columns.astype(str)
    S = [s for s in D.columns if s in env.index]
    X = np.column_stack([np.ones(len(S)), env.loc[S, "Cov"].values])
    H = X @ np.linalg.solve(X.T @ X, X.T)
    fog = env.loc[S, "fog"].values.astype(float)
    M = D[S].values.T
    beta = (M - H @ M).T @ (fog - H @ fog)
    return D[S].where(pd.Series(beta > 0, index=D.index), 2 - D[S], axis=0)


def _allowed_intervals(c, L, lo, hi):
    """Intervals on one scaffold whose distance from the centromere call is in [lo, hi]."""
    iv = []
    if lo == 0:
        iv.append((c.cen_start, c.cen_end))
    iv.append((max(0, c.cen_start - hi), max(0, c.cen_start - lo)))
    iv.append((min(L, c.cen_end + lo), min(L, c.cen_end + hi)))
    return [(a, b) for a, b in iv if b - a > REGION_BP]


def cen_null_dosages(samples, bands, peaks):
    """bands: {band_id: (lo, hi)}; peaks: DataFrame chrom/start/end of fog peaks (excluded +-flank).
    Returns SNP x sample dosages indexed by (chrom, pos, band)."""
    if os.path.exists(CEN_NULL_CACHE):
        d = pd.read_csv(CEN_NULL_CACHE, sep="\t", index_col=[0, 1, 2])
        if set(d.index.get_level_values("band")) >= set(bands):
            return d
    cen = centromeres()
    cen = cen[cen.index != EXCLUDE_CHROM]
    rng = np.random.default_rng(13)
    regs = []
    for bid, (lo, hi) in bands.items():
        iv = [(ch, a, b) for ch, c in cen.iterrows() for a, b in _allowed_intervals(c, c.len, lo, hi)]
        w = np.array([b - a for _, a, b in iv], float)
        n = 0
        while n < CEN_REGIONS_PER_BAND:
            ch, a, b = iv[rng.choice(len(iv), p=w / w.sum())]
            s = int(rng.integers(a, b - REGION_BP))
            near = peaks[(peaks.chrom == ch) & (peaks.end + EXCLUDE_FLANK > s) & (peaks.start - EXCLUDE_FLANK < s + REGION_BP)]
            if len(near):
                continue
            regs.append((ch, s, s + REGION_BP, bid)); n += 1
    reg = pd.DataFrame(regs, columns=["chrom", "start", "end", "band"]).sort_values(["chrom", "start"])
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        reg[["chrom", "start", "end"]].to_csv(fh, sep="\t", header=False, index=False); regfile = fh.name
    out = subprocess.run(["bcftools", "query", "-R", regfile, "-s", ",".join(samples),
                          "-f", "%CHROM\t%POS[\t%DS]\n", VCF_DS], capture_output=True, text=True).stdout.splitlines()
    os.remove(regfile)
    d = pd.DataFrame([l.split("\t") for l in out]).drop_duplicates([0, 1])
    d[1] = d[1].astype(int)
    rows = []
    for (ch, s, e, bid) in reg.itertuples(index=False):
        sel = d[(d[0] == ch) & (d[1] >= s) & (d[1] <= e)]
        if len(sel):
            rows.append(sel.assign(band=bid))
    d = pd.concat(rows)
    d = d.set_index([0, 1, "band"]).astype(float)
    d.index.names = ["chrom", "pos", "band"]; d.columns = samples
    d.to_csv(CEN_NULL_CACHE, sep="\t", compression="gzip")
    return d
