#!/usr/bin/env python3
"""
Test whether GWAS outlier loci lie closer to centromeres than expected.

Centromeres: confident calls from call_centromeres.py. Outlier SNPs are grouped
into peaks (<= PEAK_GAP apart) and peaks into loci (<= LOCUS_GAP apart); a locus
is represented by the mid-point of its span. Relative distance = bp to the
nearer call edge / arm length (0 = centromere, 1 = chromosome end).
Null: each locus replaced by a random tested SNP on the same scaffold, N_PERM times.

Inputs: outliers_for_annot_40000_3_<var>.txt ("retained"),
        unfiltered_outliers_for_annot_<var>.txt ("unfiltered"), scaffold_snp_map.txt
Usage: outlier_centromere_enrichment.py [var ...]   (default: fog)
Output: OUT_DIR/outlier_centromere_<var>.tsv (per locus), printed summary
"""

import os
import sys

import numpy as np
import pandas as pd

BASE_DIR = "/private/groups/meyerlab/esteiger"

GEMMA = f"{BASE_DIR}/selection/gemma"
CEN = f"{BASE_DIR}/inversions/centromere/centromeres.tsv"
OUT_DIR = f"{BASE_DIR}/inversions/centromere"
PEAK_GAP, LOCUS_GAP, NEAR_MB, N_PERM = 40_000, 2_000_000, 5, 10000
VARS = sys.argv[1:] or ["fog"]
rng = np.random.default_rng(5)


def peaks(path, smap):
    o = pd.read_csv(path, sep=" ", header=None).iloc[:, :2]
    o.columns = ["gemma_chr", "snp"]
    o = o.drop_duplicates("snp").merge(smap, on=["gemma_chr", "snp"])
    o = o.sort_values(["chr", "pos"])
    o["peak"] = ((o.chr != o.chr.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    return o.groupby("peak").agg(chr=("chr", "first"), start=("pos", "min"), end=("pos", "max"),
                                 n_snps=("pos", "size")).reset_index(drop=True)


def to_loci(pk):
    pk = pk.sort_values(["chr", "start"]).copy()
    pk["locus"] = ((pk.chr != pk.chr.shift()) | (pk.start - pk.end.shift() > LOCUS_GAP)).cumsum()
    return pk.groupby("locus").agg(chr=("chr", "first"), start=("start", "min"), end=("end", "max"),
                                   n_peaks=("start", "size"), n_snps=("n_snps", "sum")).reset_index(drop=True)


def rel_dist(chrom, pos, cen):
    c = cen.loc[chrom]
    pos = np.asarray(pos, float)
    left = pos < c.cen_start; right = pos > c.cen_end
    d = np.where(left, c.cen_start - pos, np.where(right, pos - c.cen_end, 0.0))
    arm = np.where(left, c.left_arm, np.where(right, c.right_arm, 1.0))
    return d, d / arm


def main():
    cen = pd.read_csv(CEN, sep="\t").query("confident").set_index("chr")
    smap = pd.read_csv(f"{GEMMA}/scaffold_snp_map.txt", sep=r"\s+", header=None,
                       names=["gemma_chr", "snp", "chr", "pos"])
    bg = {c: g.pos.values for c, g in smap[smap.chr.isin(cen.index)].groupby("chr")}
    allbg = np.concatenate([rel_dist(c, p, cen)[1] for c, p in bg.items()])
    alld = np.concatenate([rel_dist(c, p, cen)[0] for c, p in bg.items()])
    print(f"centromeres: {len(cen)} scaffolds; background SNPs {len(allbg):,}: mean relative distance "
          f"{allbg.mean():.3f}, within {NEAR_MB} Mb of a centromere {np.mean(alld <= NEAR_MB * 1e6):.3f}")
    summary = []
    for var in VARS:
        for kind, fname in [("retained", f"outliers_for_annot_40000_3_{var}.txt"),
                            ("unfiltered", f"unfiltered_outliers_for_annot_{var}.txt")]:
            path = f"{GEMMA}/{var}/{fname}"
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                continue
            pk = peaks(path, smap)
            if kind == "retained":
                pk = pk[pk.n_snps >= 3]
            n_all = len(pk)
            loci = to_loci(pk[pk.chr.isin(cen.index)])
            if len(loci) < 3:
                continue
            loci["mid"] = (loci.start + loci.end) / 2
            loci["dist_bp"], loci["rel_dist"] = 0.0, 0.0
            for c, g in loci.groupby("chr"):
                d, r = rel_dist(c, g.mid.values, cen)
                loci.loc[g.index, "dist_bp"], loci.loc[g.index, "rel_dist"] = d, r
            obs_rel, obs_near = loci.rel_dist.mean(), np.mean(loci.dist_bp <= NEAR_MB * 1e6)
            null_rel, null_near = np.empty(N_PERM), np.empty(N_PERM)
            chroms = loci.chr.values
            draws = {c: rng.choice(bg[c], size=(N_PERM, (chroms == c).sum())) for c in np.unique(chroms)}
            dr = {c: rel_dist(c, draws[c].ravel(), cen) for c in draws}
            D = np.hstack([dr[c][0].reshape(N_PERM, -1) for c in draws])
            R = np.hstack([dr[c][1].reshape(N_PERM, -1) for c in draws])
            null_rel, null_near = R.mean(1), (D <= NEAR_MB * 1e6).mean(1)
            p_rel = (np.sum(null_rel <= obs_rel) + 1) / (N_PERM + 1)
            p_near = (np.sum(null_near >= obs_near) + 1) / (N_PERM + 1)
            loci.to_csv(f"{OUT_DIR}/outlier_centromere_{var}_{kind}.tsv", sep="\t", index=False)
            summary.append(dict(var=var, outliers=kind, peaks_all=n_all, loci_used=len(loci),
                                obs_mean_rel_dist=obs_rel, null_mean_rel_dist=null_rel.mean(), p_closer=p_rel,
                                obs_frac_within_5Mb=obs_near, null_frac_within_5Mb=null_near.mean(),
                                p_within_5Mb=p_near))
            if var == "fog":
                pd.set_option("display.width", 200)
                print(f"\n{var} {kind}: {n_all} peaks, {len(loci)} loci on scaffolds with a centromere call")
                print(loci.assign(mid_Mb=loci.mid / 1e6, dist_Mb=loci.dist_bp / 1e6)[
                    ["chr", "mid_Mb", "n_peaks", "n_snps", "dist_Mb", "rel_dist"]].to_string(
                    index=False, float_format=lambda x: f"{x:.2f}"))
    s = pd.DataFrame(summary)
    s.to_csv(f"{OUT_DIR}/outlier_centromere_summary.tsv", sep="\t", index=False)
    print("\n" + s.to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
