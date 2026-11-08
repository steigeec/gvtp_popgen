#!/usr/bin/env python3
"""
Per-peak tables and an exploratory figure for the fog peaks within
Tranquillon Mtn (Strauss), using the same data, polarization and spatial
units as fog_hypothesis_tests.py (736 plants after excluding
untrustworthy-coordinate plants; stands from field_data.txt with stand 88
split by long-term fog cell; units < 5 plants dropped; favorable allele = the
one whose coverage-adjusted dosage rises with long-term fog across all samples).

Tables, one row set per peak (frequency = mean favorable dosage / 2 over the
peak's SNPs), plotted in DEINVI_0006_stacked.py:
  OUT_UNITS: unit frequency and unit long-term fog, with the n-weighted
     least-squares line (slope_longterm, intercept_longterm).
  OUT_YEARS: unit-adjusted frequency per germination cohort (+-95% CI, 1.96 SE
     across plants) and that cohort's mean June fog, with the within-unit
     June-fog line (freq = mean_freq + slope_june * (fog - june_mean)).
OUT_PLANTS ("chr 6, plant by plant"): every plant's chr 6 favorable-allele
  frequency (dosage / 2; jittered vertically for display) against the June fog
  of its germination year at its location, with a separate least-squares line
  per germination year (coloured by year).

Usage: fog_peak_figures.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import temporal_null as tn

BASE_DIR = "/private/groups/meyerlab/esteiger"

G = f"{BASE_DIR}/selection/gemma"
OUTLIERS = f"{G}/fog/outliers_for_annot_40000_3_fog.txt"
DOSAGE = f"{G}/fog/dosage_fog.txt"
FIELD = f"{BASE_DIR}/field_data.txt"
OUT_UNITS = f"{G}/fog/fog_peaks_by_unit.tsv"     # read by DEINVI_0006_stacked.py (panel H)
OUT_YEARS = f"{G}/fog/fog_peaks_by_year.tsv"     # read by DEINVI_0006_stacked.py (panel J)
OUT_PLANTS = f"{G}/fog/fog_chr6_plants_june_fog.jpg"
CHR6_SPAN = (6, 49_031_257, 49_146_391)
SPLIT_STAND, MIN_UNIT, PEAK_GAP = 88, 5, 40_000
LABEL_SIZE, TICK_SIZE = 18, 14
YEAR_COL = {2020: "#999999", 2021: "#000000", 2022: "#d55e00"}
rng = np.random.default_rng(3)


def wfit(x, y, w):
    """n-weighted least-squares slope and intercept"""
    w = w / w.sum(); mx, my = np.sum(w * x), np.sum(w * y)
    b = np.sum(w * (x - mx) * (y - my)) / np.sum(w * (x - mx) ** 2)
    return b, my - b * mx


def main():
    env = pd.read_csv(tn.ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    field = pd.read_csv(FIELD, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    ds = pd.read_csv(DOSAGE, sep="\t", index_col=0); ds.index = ds.index.astype(str); ds.columns = ds.columns.astype(str)
    o = pd.read_csv(OUTLIERS, sep=" ", header=None).iloc[:, :3]; o.columns = ["chr", "snp", "pos"]
    o = o.drop_duplicates("snp").sort_values(["chr", "pos"]); o["snp"] = o.snp.astype(str)
    o["peak"] = ((o.chr != o.chr.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    c6 = o[(o.chr == CHR6_SPAN[0]) & o.pos.between(*CHR6_SPAN[1:])].peak.unique()[0]

    fav = tn.polarize(ds.loc[o.snp])
    st = [s for s in fav.columns if env.loc[s, "Location"] == "Strauss" and s in field.index]
    unit = field.loc[st, "OHANA_order"].astype(int).astype(str)
    sp = unit == str(SPLIT_STAND)
    unit[sp] = [f"{SPLIT_STAND}_fog{v:.2f}" for v in env.loc[unit.index[sp], "fog"]]
    st = list(unit.index[unit.groupby(unit).transform("size") >= MIN_UNIT]); unit = unit[st]
    units = np.array(sorted(unit.unique())); U = np.array([np.where(units == u)[0][0] for u in unit])
    n_u = np.bincount(U).astype(float)
    u_lt = np.bincount(U, weights=env.loc[st, "fog"].astype(float).values) / n_u
    june = env.loc[st, "fog_germ_yr"].astype(float).values
    year = env.loc[st, "Year"].astype(int).values
    jc = june - (np.bincount(U, weights=june) / n_u)[U]
    years = sorted(set(year))
    print(f"{len(st)} plants, {len(units)} units, {o.peak.nunique()} peaks")

    # ---------------- tables: every peak separately ----------------
    smap = pd.read_csv(f"{G}/scaffold_snp_map.txt", sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    scaf = dict(zip(smap[0], smap[2]))
    u_rows, y_rows, rows = [], [], []
    for pk, g in o.groupby("peak"):
        F = fav.loc[g.snp, st].values.mean(0) / 2                    # plant favorable-allele frequency
        uf = np.bincount(U, weights=F) / n_u
        b1, a1 = wfit(u_lt, uf, n_u)
        Fc = F - uf[U]; b2 = Fc @ jc / (jc @ jc)
        mean = np.average(uf, weights=n_u)
        adj = Fc + mean
        info = dict(peak=pk, scaffold=scaf[g.chr.iloc[0]], start=g.pos.min(), end=g.pos.max(), chr6=pk == c6)
        u_rows += [dict(info, unit=u, n=int(n), longterm_fog=x, freq=f, slope_longterm=b1, intercept_longterm=a1)
                   for u, n, x, f in zip(units, n_u, u_lt, uf)]
        # per year: June fog, unit-adjusted mean frequency, 95% CI half-width (1.96 SE across plants)
        y_rows += [dict(info, year=y, n=int((year == y).sum()), june_fog=june[year == y].mean(),
                        unit_adjusted_freq=adj[year == y].mean(),
                        ci95_half=1.96 * adj[year == y].std(ddof=1) / np.sqrt((year == y).sum()),
                        slope_june=b2, mean_freq=mean, june_mean=june.mean()) for y in years]
        rows.append(dict(peak=pk, start=g.pos.min(), slope_longterm=b1, slope_june=b2))
    r = pd.DataFrame(rows)
    for col, lbl in [("slope_longterm", "long-term fog"), ("slope_june", "June fog")]:
        rk = int((r[col] > r.loc[r.peak == c6, col].iloc[0]).sum()) + 1
        print(f"{lbl}: {(r[col] > 0).sum()} of {len(r)} peaks have a positive slope; chr 6 ranks {rk} of {len(r)}")
    pd.DataFrame(u_rows).to_csv(OUT_UNITS, sep="\t", index=False)
    pd.DataFrame(y_rows).to_csv(OUT_YEARS, sep="\t", index=False)
    print(f"Saved: {OUT_UNITS}\nSaved: {OUT_YEARS}")

    # ---------------- figure 2: chr 6, every plant ----------------
    F6 = fav.loc[o[o.peak == c6].snp, st].values.mean(0) / 2
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for y in years:
        m = year == y
        ax.scatter(june[m] + rng.uniform(-0.12, 0.12, m.sum()), F6[m] + rng.uniform(-0.03, 0.03, m.sum()), s=14,
                   color=YEAR_COL[y], alpha=0.45, edgecolor="none")
        b_, a_ = np.polyfit(june[m], F6[m], 1)
        xx = np.linspace(june[m].min(), june[m].max(), 2)
        ax.plot(xx, a_ + b_ * xx, color=YEAR_COL[y], lw=3, label=f"{y} (n = {m.sum()}; slope {b_:+.3f})")
        print(f"{y}: June fog {june[m].min():.0f}-{june[m].max():.0f} days/month; mean freq {F6[m].mean():.3f}; slope {b_:+.4f} per day/month")
    ax.set_xlabel("June fog of germination year (days/month)", fontsize=LABEL_SIZE)
    ax.set_ylabel("chr 6 favorable-allele frequency\n(per plant)", fontsize=LABEL_SIZE)
    ax.tick_params(labelsize=TICK_SIZE); ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=TICK_SIZE - 3, title="Germination year", title_fontsize=TICK_SIZE - 3)
    fig.savefig(OUT_PLANTS, dpi=170, bbox_inches="tight")
    print(f"Saved: {OUT_PLANTS}")


if __name__ == "__main__":
    main()
