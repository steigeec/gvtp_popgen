#!/usr/bin/env python3
"""
Does the chr 6 fog-peak haplotype (and each other fog peak) track long-term
fog WITHIN Tranquillon Mtn (Location "Strauss"), i.e. without the
between-locality contrast (e.g. Pt. Arguello vs the rest)?

Long-term fog varies little within the site (6.9-8.4 days/month, 10 distinct
values) and mostly between stands (OHANA_order), so the unit is the stand.
Statistic: n-weighted Pearson correlation, across stands, between the
stand's favorable-allele frequency (mean favorable dosage / 2 over the
peak's SNPs) and the stand's mean fog. Stands with < MIN_STAND plants dropped.
  - Stand permutation: stand fog values shuffled among stands (N_PERM).
  - Random-SNP null: the same correlation for fog-polarized random SNPs
    (temporal null cache; one per 2 kb region), frequency-matched to the
    peak within the site (+-FREQ_TOL). Fog is spatially structured, so
    random SNPs polarized by fog also correlate with it somewhat; this null
    calibrates for that background.
Caveat: these plants were part of the GWAS that found the peaks, so this is
not an independent replication; it asks whether the signal exists within
the site at all, or comes only from contrasts between localities.
Untrustworthy-coordinate individuals excluded (temporal_null.BAD_COORD_IDS).

Usage: within_strauss_fog.py
Output: OUT (one row per peak)
"""

import numpy as np
import pandas as pd

import temporal_null as tn

G = "/private/groups/meyerlab/esteiger/selection/gemma"
OUTLIERS = f"{G}/fog/outliers_for_annot_40000_3_fog.txt"
DOSAGE = f"{G}/fog/dosage_fog.txt"
FIELD = "/private/groups/meyerlab/esteiger/field_data.txt"
NULL_CACHE = f"{G}/fog/temporal_null_random_snps_allsamples.tsv.gz"
OUT = f"{G}/fog/within_strauss_fog.tsv"
MIN_STAND, FREQ_TOL, N_PERM, PEAK_GAP = 5, 0.04, 10000, 40_000
SPLIT_STAND = "88"
rng = np.random.default_rng(21)


def wcorr(x, y, w):
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    return np.sum(w * (x - mx) * (y - my)) / np.sqrt(np.sum(w * (x - mx) ** 2) * np.sum(w * (y - my) ** 2))


def wcorr_rows(Y, x, w):
    """weighted correlation of each row of Y (features x stands) with x"""
    w = w / w.sum(); xc = x - np.sum(w * x)
    Yc = Y - (Y * w).sum(1, keepdims=True)
    return (Yc * w * xc).sum(1) / np.sqrt((Yc ** 2 * w).sum(1) * np.sum(w * xc ** 2))


def main():
    ds = pd.read_csv(DOSAGE, sep="\t", index_col=0); ds.index = ds.index.astype(str); ds.columns = ds.columns.astype(str)
    env = pd.read_csv(tn.ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    field = pd.read_csv(FIELD, sep=None, engine="python", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    o = pd.read_csv(OUTLIERS, sep=" ", header=None).iloc[:, :3]; o.columns = ["chr", "snp", "pos"]
    o = o.drop_duplicates("snp").sort_values(["chr", "pos"])
    o["peak"] = ((o.chr != o.chr.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    smap = pd.read_csv(f"{G}/scaffold_snp_map.txt", sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    n2s = dict(zip(smap[0], smap[2]))

    fav = tn.polarize(ds.loc[o.snp.astype(str)])                  # drops BAD_COORD_IDS
    st = [s for s in fav.columns if env.loc[s, "Location"] == "Strauss"]
    stand = field["OHANA_order"].reindex(st)
    stand = stand.map(lambda v: f"{int(v)}" if pd.notna(v) else v)
    # Stand SPLIT_STAND is a ~1.7 km strip spanning several fog grid cells
    # (6.89-7.86 days/month); for these fog tests only, it is split into one
    # stand per fog cell (field_data.txt itself is not changed)
    split = stand == SPLIT_STAND
    stand[split] = [f"{SPLIT_STAND}_fog{v:.2f}" for v in env.loc[stand.index[split], "fog"].astype(float)]
    print("stand", SPLIT_STAND, "split by fog cell:", stand[split].value_counts().to_dict())
    keep = stand.groupby(stand).transform("size") >= MIN_STAND
    st = [s for s, k in zip(st, keep) if k and pd.notna(stand[s])]
    stand = stand[st]
    stands = np.array(sorted(stand.unique()))
    sfog = env.loc[st, "fog"].astype(float).groupby(stand.values).mean().reindex(stands).values
    n = stand.value_counts().reindex(stands).values.astype(float)
    print(f"{len(st)} plants in {len(stands)} stands (>= {MIN_STAND} plants); stand fog "
          f"{sfog.min():.2f}-{sfog.max():.2f} days/month ({len(np.unique(sfog))} distinct)")
    perm_fog = np.array([rng.permutation(sfog) for _ in range(N_PERM)])

    def stand_freq(D):                                            # D: features x plants (favorable dose)
        return np.column_stack([D[:, (stand.values == s)].mean(1) / 2 for s in stands])

    nul = tn.polarize(pd.read_csv(NULL_CACHE, sep="\t", index_col=[0, 1]))[st]
    region = np.array([f"{c}:{int(p) // 2000}" for c, p in nul.index])
    NF = stand_freq(nul.values)
    nf = np.average(NF, axis=1, weights=n)
    poly = (nf > 0) & (nf < 1) & (NF.std(1) > 0)
    NF, nf, region = NF[poly], nf[poly], region[poly]
    nr = wcorr_rows(NF, sfog, n)

    rows = []
    for pk, g in o.groupby("peak"):
        F = stand_freq(fav.loc[g.snp.astype(str), st].values.mean(0, keepdims=True))[0]
        f = np.average(F, weights=n); r = wcorr(sfog, F, n)
        rp = np.array([wcorr(pf, F, n) for pf in perm_fog])
        p_perm = (np.sum(rp >= r) + 1) / (N_PERM + 1)
        m = np.where(np.abs(nf - f) <= FREQ_TOL)[0]
        m = pd.Series(m).groupby(region[m]).apply(lambda x: x.sample(1, random_state=1)).values
        p_null = (np.sum(nr[m] >= r) + 1) / (len(m) + 1)
        rows.append(dict(scaffold=n2s[g.chr.iloc[0]], start=g.pos.min(), end=g.pos.max(), n_snps=len(g),
                         fav_freq=round(f, 3), r_stand_fog=r, p_stand_perm=p_perm, n_null=len(m),
                         null_median_r=np.median(nr[m]), p_random_snp=p_null))
    res = pd.DataFrame(rows).sort_values("r_stand_fog", ascending=False)
    res.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(f"random fog-polarized SNPs: median stand-level r with fog {np.median(nr):+.3f}")
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    c6 = res[res.scaffold == "DEINVI_00006"].iloc[0]
    print(f"\nchr 6 haplotype: r = {c6.r_stand_fog:+.3f}, stand permutation p = {c6.p_stand_perm:.3f}, "
          f"vs random SNPs p = {c6.p_random_snp:.3f}; rank {int((res.r_stand_fog > c6.r_stand_fog).sum()) + 1} of {len(res)}")
    print("per-stand chr 6 haplotype frequency (by stand fog):")
    F6 = stand_freq(fav.loc[o[(o.chr == 6) & o.pos.between(49_031_257, 49_146_391)].snp.astype(str), st].values.mean(0, keepdims=True))[0]
    print(pd.DataFrame({"stand": stands, "n": n.astype(int), "fog": sfog.round(2), "hap_freq": F6.round(3)})
          .sort_values("fog").to_string(index=False))


if __name__ == "__main__":
    main()
