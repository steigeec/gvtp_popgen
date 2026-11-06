#!/usr/bin/env python3
"""
Tests of two hypotheses for fog-favorable alleles within Tranquillon Mtn
(Location "Strauss"), at two levels:
  - chr 6: the chr 6 fog-peak haplotype (mean favorable dosage over its SNPs)
  - polygenic: the other fog peaks (each peak = mean favorable dosage over its
    SNPs; score = mean over peaks, so every peak counts once; chr 6 left out
    so the two levels are independent)
Favorable = allele whose coverage-adjusted dosage increases with long-term
fog across all samples (temporal_null.polarize; same rule as the outliers).
Untrustworthy-coordinate individuals excluded (temporal_null.BAD_COORD_IDS).

Spatial units: stands (field_data.txt OHANA_order), with SPLIT_STAND split
by long-term fog grid cell (it spans several cells); units < MIN_UNIT plants
dropped.

H1 local adaptation (space, long-term fog = 2000-2022 Jun-Sep mean):
  statistic = n-weighted least-squares slope of unit favorable-allele
  frequency on unit long-term fog (frequency per fog day/month).
H2 differential germination / early mortality (time, germination-year fog
  = June fog of the sampling year; plants sampled May-June):
  statistic = slope of plant favorable-allele dose / 2 on June fog centred
  within its unit (unit fixed effects), so only changes over time within a
  place count.
Nulls (one-sided, H1/H2 predict positive slopes):
  - permutation: H1 long-term fog shuffled among units; H2 June fog shuffled
    among plants within units (~ shuffling years within units) -- N_PERM.
  - random SNPs: fog-polarized random SNPs (temporal null cache, one per 2 kb
    region), frequency-matched (+-FREQ_TOL within the site). chr 6: each SNP
    is one null draw. Polygenic: each draw averages one matched random SNP
    per peak (N_SET draws). The statistics are linear in dosage, so set
    statistics are averages of per-SNP statistics.
  - chr 6 only, centromere-matched random SNPs: as above, but drawn at a
    similar distance from their own scaffold's centromere as the chr 6 peak
    is from the chr 6 centromere (temporal_null.cen_band; same cache and
    band as fog_peaks_temporal.py / DEINVI_0006_stacked.py), since the peak
    sits ~0.6 Mb from the centromere, where low recombination can make
    allele frequencies more variable.

Usage: fog_hypothesis_tests.py
Output: OUT_TSV (tests), OUT_FIG (plot), OUT_BY_YEAR (H2 per-cohort values), printed summary
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import temporal_null as tn

G = "/private/groups/meyerlab/esteiger/selection/gemma"
OUTLIERS = f"{G}/fog/outliers_for_annot_40000_3_fog.txt"
DOSAGE = f"{G}/fog/dosage_fog.txt"
FIELD = "/private/groups/meyerlab/esteiger/field_data.txt"
NULL_CACHE = f"{G}/fog/temporal_null_random_snps_allsamples.tsv.gz"
OUT_TSV = f"{G}/fog/fog_hypothesis_tests.tsv"
OUT_FIG = f"{G}/fog/fog_hypothesis_tests.jpg"
OUT_BY_YEAR = f"{G}/fog/fog_hypothesis_tests_by_year.tsv"   # H2 per-cohort values (main figure panel H)
CHR6_SPAN = (6, 49_031_257, 49_146_391)
SPLIT_STAND, MIN_UNIT, PEAK_GAP = 88, 5, 40_000
FREQ_TOL, N_PERM, N_SET = 0.04, 10000, 10000
LABEL_SIZE, TICK_SIZE = 18, 14
COL = {"chr 6 haplotype": "#d55e00", "polygenic (other peaks)": "#0072b2"}
rng = np.random.default_rng(31)


def wslope(x, y, w):
    """n-weighted LS slope of each row of y (features x units) on x"""
    w = w / w.sum(); xc = x - np.sum(w * x)
    return ((y - (y * w).sum(-1, keepdims=True)) * w * xc).sum(-1) / np.sum(w * xc ** 2)


def main():
    env = pd.read_csv(tn.ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    field = pd.read_csv(FIELD, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    ds = pd.read_csv(DOSAGE, sep="\t", index_col=0); ds.index = ds.index.astype(str); ds.columns = ds.columns.astype(str)
    o = pd.read_csv(OUTLIERS, sep=" ", header=None).iloc[:, :3]; o.columns = ["chr", "snp", "pos"]
    o = o.drop_duplicates("snp").sort_values(["chr", "pos"]); o["snp"] = o.snp.astype(str)
    o["peak"] = ((o.chr != o.chr.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    c6 = o[(o.chr == CHR6_SPAN[0]) & o.pos.between(*CHR6_SPAN[1:])].peak.unique()
    assert len(c6) == 1

    fav = tn.polarize(ds.loc[o.snp])                              # drops BAD_COORD_IDS
    st = [s for s in fav.columns if env.loc[s, "Location"] == "Strauss" and s in field.index]
    unit = field.loc[st, "OHANA_order"].astype(int).astype(str)
    sp = unit == str(SPLIT_STAND)
    unit[sp] = [f"{SPLIT_STAND}_fog{v:.2f}" for v in env.loc[unit.index[sp], "fog"]]
    big = unit.groupby(unit).transform("size") >= MIN_UNIT
    st = list(unit.index[big]); unit = unit[st]
    units = np.array(sorted(unit.unique())); U = np.array([np.where(units == u)[0][0] for u in unit])
    n_u = np.bincount(U).astype(float)
    lt = env.loc[st, "fog"].astype(float).values
    u_lt = np.bincount(U, weights=lt) / n_u                       # unit long-term fog
    june = env.loc[st, "fog_germ_yr"].astype(float).values
    year = env.loc[st, "Year"].astype(int).values
    jc = june - (np.bincount(U, weights=june) / n_u)[U]           # June fog centred within unit
    print(f"{len(st)} plants, {len(units)} units; unit long-term fog {u_lt.min():.2f}-{u_lt.max():.2f}; "
          f"June fog by year: " + ", ".join(f"{y} {june[year == y].mean():.2f}" for y in sorted(set(year))))

    def unit_freq(F):                                             # F: features x plants (allele dose/2)
        return np.stack([np.bincount(U, weights=f) / n_u for f in F])

    def h2_slope(F):
        Fc = F - (unit_freq(F))[:, U]
        return Fc @ jc / (jc @ jc)

    # per-peak favorable allele frequency (dose / 2), plants in st
    peak_F = {pk: fav.loc[g.snp, st].values.mean(0) / 2 for pk, g in o.groupby("peak")}
    levels = {"chr 6 haplotype": [c6[0]], "polygenic (other peaks)": [p for p in peak_F if p != c6[0]]}

    # random-SNP pool: per SNP site frequency, unit frequencies, H1 and H2 slopes
    nul = tn.polarize(pd.read_csv(NULL_CACHE, sep="\t", index_col=[0, 1]))[st]
    region = np.array([f"{c}:{int(p) // 2000}" for c, p in nul.index])
    NF = nul.values / 2
    nf = NF.mean(1); poly = (nf > 0) & (nf < 1)
    NF, nf, region = NF[poly], nf[poly], region[poly]
    NU = unit_freq(NF); n_h1 = wslope(u_lt, NU, n_u); n_h2 = h2_slope(NF)
    print(f"random-SNP pool: {len(nf)} SNPs")

    def matched(f, nf=nf, region=region):
        m = np.where(np.abs(nf - f) <= FREQ_TOL)[0]
        return pd.Series(m).groupby(region[m]).apply(lambda x: x.sample(1, random_state=1)).values

    # centromere-matched random SNPs for the chr 6 peak
    smap = pd.read_csv(f"{G}/scaffold_snp_map.txt", sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    chrom6 = dict(zip(smap[0], smap[2]))[CHR6_SPAN[0]]
    dist = tn.cen_dist(chrom6, (CHR6_SPAN[1] + CHR6_SPAN[2]) / 2)
    band_id = f"{chrom6}:{CHR6_SPAN[1]}"
    fog_peaks = o.assign(chrom=o.chr.map(dict(zip(smap[0], smap[2])))).groupby("peak").agg(
        chrom=("chrom", "first"), start=("pos", "min"), end=("pos", "max"))
    C = tn.cen_null_dosages(list(env.index), {band_id: tn.cen_band(dist)}, fog_peaks)
    C = tn.polarize(C[C.index.get_level_values("band") == band_id])[st]
    CF = C.values / 2; cf = CF.mean(1); cpoly = (cf > 0) & (cf < 1)
    c_region = np.array([f"{c}:{int(p) // 2000}" for c, p, _ in C.index])[cpoly]
    CF, cf = CF[cpoly], cf[cpoly]
    c_h1, c_h2 = wslope(u_lt, unit_freq(CF), n_u), h2_slope(CF)
    lo, hi = tn.cen_band(dist)
    print(f"centromere-matched pool for chr 6 (peak {dist / 1e6:.2f} Mb from its centromere; band "
          f"{lo / 1e6:.2f}-{hi / 1e6:.2f} Mb): {len(cf)} SNPs")

    rows, plot = [], {}
    for lvl, pks in levels.items():
        F = np.mean([peak_F[p] for p in pks], axis=0)[None, :]    # plant-level score (freq scale)
        uf = unit_freq(F)[0]
        h1, h2 = wslope(u_lt, uf[None], n_u)[0], h2_slope(F)[0]
        # permutations
        p1 = np.array([wslope(rng.permutation(u_lt), uf[None], n_u)[0] for _ in range(N_PERM)])
        Fc = F[0] - uf[U]
        p2 = np.empty(N_PERM)
        groups = [np.where(U == u)[0] for u in range(len(units))]
        for k in range(N_PERM):
            jp = jc.copy()
            for g in groups: jp[g] = rng.permutation(jc[g])
            p2[k] = Fc @ jp / (jc @ jc)
        # random-SNP nulls
        c1 = c2 = None
        if len(pks) == 1:
            m = matched(F.mean()); r1, r2 = n_h1[m], n_h2[m]
            mc = matched(F.mean(), cf, c_region); c1, c2 = c_h1[mc], c_h2[mc]
        else:
            pools = [matched(peak_F[p].mean()) for p in pks]
            draw = np.stack([rng.choice(pl, N_SET) for pl in pools])     # peaks x draws
            r1, r2 = n_h1[draw].mean(0), n_h2[draw].mean(0)
        for hyp, obs, perm, rnd, crn in [("H1 long-term fog (space)", h1, p1, r1, c1),
                                         ("H2 June fog (time, within units)", h2, p2, r2, c2)]:
            rows.append(dict(level=lvl, hypothesis=hyp, n_peaks=len(pks), slope=obs,
                             p_permutation=(np.sum(perm >= obs) + 1) / (len(perm) + 1),
                             p_random_snps=(np.sum(rnd >= obs) + 1) / (len(rnd) + 1), n_null=len(rnd),
                             null_95pct=np.percentile(rnd, 95),
                             p_cen_matched=(np.sum(crn >= obs) + 1) / (len(crn) + 1) if crn is not None else np.nan,
                             n_null_cen=len(crn) if crn is not None else 0))
        # plot data: H1 unit frequencies; H2 unit-adjusted year frequencies with SE
        adj = F[0] - uf[U] + np.average(uf, weights=n_u)
        yrs = sorted(set(year))
        plot[lvl] = dict(uf=uf, h1=h1, mean=np.average(uf, weights=n_u),
                         yr=[(y, june[year == y].mean(), adj[year == y].mean(),
                              1.96 * adj[year == y].std(ddof=1) / np.sqrt((year == y).sum())) for y in yrs], h2=h2)
    res = pd.DataFrame(rows)
    res.to_csv(OUT_TSV, sep="\t", index=False)
    pd.DataFrame([dict(level=lvl, year=int(y), june_fog=jf, unit_adjusted_freq=f, ci95_half=ci)
                  for lvl, d in plot.items() for y, jf, f, ci in d["yr"]]).to_csv(OUT_BY_YEAR, sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # ── figure ──
    fig, (a, b) = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw=dict(wspace=0.3))
    xs = np.linspace(u_lt.min() - 0.05, u_lt.max() + 0.05, 2)
    for lvl, d in plot.items():
        c = COL[lvl]; mk = "D" if lvl.startswith("chr") else "o"
        a.scatter(u_lt, d["uf"], s=n_u * 1.2 + 15, marker=mk, color=c, edgecolor="black", alpha=0.85, label=lvl, zorder=3)
        xm = np.average(u_lt, weights=n_u)
        a.plot(xs, d["mean"] + d["h1"] * (xs - xm), color=c, lw=2)
        yy = np.array(d["yr"])
        b.errorbar(yy[:, 1], yy[:, 2], yerr=yy[:, 3], fmt=mk, ms=11, color=c, mec="black", capsize=6, lw=2, label=lvl)
        jm = np.mean(june); jx = np.linspace(min(yy[:, 1]) - 0.3, max(yy[:, 1]) + 0.3, 2)
        b.plot(jx, d["mean"] + d["h2"] * (jx - jm), color=c, lw=1.5, ls="--")
        if lvl.startswith("chr"):
            for y, x, f, _ in d["yr"]:
                b.annotate(str(int(y)), (x, f), xytext=(8, 8), textcoords="offset points", fontsize=TICK_SIZE - 2)
    for ax, (hyp, xl, t) in zip((a, b), [("H1 long-term fog (space)", "Long-term fog of unit (days/month)",
                                         "Local adaptation: units within Tranquillon Mtn"),
                                        ("H2 June fog (time, within units)", "June fog of germination year (days/month)",
                                         "Cohort response: adjusted for unit")]):
        txt = []
        for lvl in COL:
            r = res[(res.level == lvl) & (res.hypothesis == hyp)].iloc[0]
            cen_txt = f"; centromere-matched {r.p_cen_matched:.3f}" if r.p_cen_matched == r.p_cen_matched else ""
            txt.append(f"{lvl.split(' (')[0]}: slope {r.slope:+.4f}\n   p: permutation {r.p_permutation:.3f}; "
                       f"random SNPs {r.p_random_snps:.3f}{cen_txt}")
        ax.text(0.02, 0.98, "\n".join(txt), transform=ax.transAxes, va="top", fontsize=TICK_SIZE - 4)
        ax.set_xlabel(xl, fontsize=LABEL_SIZE); ax.set_title(t, fontsize=LABEL_SIZE - 2)
        ax.set_ylabel("Favorable-allele frequency", fontsize=LABEL_SIZE)
        ax.tick_params(labelsize=TICK_SIZE); ax.spines[["top", "right"]].set_visible(False)
    for ax in (a, b): ax.set_ylim(0, ax.get_ylim()[1] * 1.3)   # headroom for the annotation text
    a.legend(frameon=False, fontsize=TICK_SIZE - 2, loc="upper left", bbox_to_anchor=(0.0, 0.82))
    for ax, l in [(a, "A"), (b, "B")]:
        ax.text(-0.14, 1.04, l, transform=ax.transAxes, fontsize=24, va="bottom")
    fig.savefig(OUT_FIG, dpi=170, bbox_inches="tight")
    print(f"Saved: {OUT_FIG}\n-> {OUT_TSV}")


if __name__ == "__main__":
    main()
