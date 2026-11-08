#!/usr/bin/env python3
"""
Sliding-window mapping-quality check along DEINVI_00006, testing whether the
fog GWAS peak (GNS1/SNF1, ~49.03-49.15 Mb) shows signs of ambiguous mapping.

From the merged BAMs of the 339 all-locality samples, summed across samples on
a 10 kb grid:
  total_depth     samtools depth -Q 0
  highqual_depth  samtools depth -Q 20
  frac_lowmapq    1 - highqual_depth / total_depth
Test: one-sided Mann-Whitney U (peak > background) of frac_lowmapq, with peak
and centromere windows excluded from the background. Runs in N_CHUNKS
parallel chunks per depth pass.

Usage: compute_deinvi6_mapq.py
Output: OUT (pos, total_depth, highqual_depth, frac_lowmapq; tsv.gz),
        TEST_OUT (test result); both also printed to stdout.
"""
from concurrent.futures import ThreadPoolExecutor
import subprocess

import pandas as pd
from scipy.stats import mannwhitneyu

BASE_DIR = "/private/groups/meyerlab/esteiger"

BAM_DIR = f"{BASE_DIR}/mapping/merged"
PC1 = f"{BASE_DIR}/inversions/all_localities/DEINVI_00006_GL_maf5_flipped.pc_1.tsv.gz"
FAI = f"{BASE_DIR}/inversions/filtered_scaffolds.fai"
OUT = f"{BASE_DIR}/inversions/all_localities/DEINVI_00006_mapq_windows.tsv.gz"
TEST_OUT = f"{BASE_DIR}/inversions/all_localities/DEINVI_00006_mapq_peak_test.txt"
CHROM = "DEINVI_00006"
WIN = 10_000
N_CHUNKS = 16
HIGHQUAL_MAPQ = 20
CENTROMERE = (45_000_000, 48_500_000)
PEAK = (49_031_000, 49_147_000)


def run_chunk(args):
    i, start, end, samples, min_mapq, label = args
    awk = (
        "awk -v w=%d '{"
        "s=0; for(i=3;i<=NF;i++) s+=$i; "
        "win=int($2/w)*w; sum[win]+=s"
        "} END {for (k in sum) print k, sum[k]}'" % WIN
    )
    bams = " ".join(f"{BAM_DIR}/{s}" for s in samples)
    cmd = f"samtools depth -a -Q {min_mapq} -r {CHROM}:{start}-{end} {bams} | {awk}"
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    print(f"[{label}] chunk {i} ({start}-{end}) done, {len(out.splitlines())} windows", flush=True)
    return out


def depth_track(samples, length, min_mapq, label):
    bounds = [1 + i * length // N_CHUNKS for i in range(N_CHUNKS)] + [length]
    jobs = [(i, bounds[i], bounds[i + 1], samples, min_mapq, label) for i in range(N_CHUNKS)]
    with ThreadPoolExecutor(max_workers=N_CHUNKS) as ex:
        chunk_outputs = list(ex.map(run_chunk, jobs))
    totals = {}
    for out in chunk_outputs:
        for line in out.splitlines():
            win, s = line.split()
            totals[int(win)] = totals.get(int(win), 0) + int(s)
    return pd.Series(totals, name=label).sort_index()


def summarize(df, lo, hi, label):
    d = df[(df["pos"] >= lo) & (df["pos"] <= hi)]
    return (f"{label}: n_windows={len(d)}, median frac_lowmapq={d.frac_lowmapq.median():.4f}, "
            f"max frac_lowmapq={d.frac_lowmapq.max():.4f}")


def test_peak_vs_background(df):
    """One-sided Mann-Whitney U of peak vs. background frac_lowmapq (background excludes peak and centromere)."""
    in_peak = (df["pos"] >= PEAK[0]) & (df["pos"] <= PEAK[1])
    in_cen = (df["pos"] >= CENTROMERE[0]) & (df["pos"] <= CENTROMERE[1])
    peak_vals = df.loc[in_peak, "frac_lowmapq"].dropna()
    background_vals = df.loc[~in_peak & ~in_cen, "frac_lowmapq"].dropna()
    u, p = mannwhitneyu(peak_vals, background_vals, alternative="greater")
    lines = [
        summarize(df, *PEAK, "Peak"),
        summarize(df, *CENTROMERE, "Centromere"),
        summarize(df, df["pos"].min(), df["pos"].max(), "Genome-wide"),
        "",
        "One-sided Mann-Whitney U test: peak frac_lowmapq > background "
        "(background = genome-wide, excluding peak and centromere windows)",
        f"  n_peak={len(peak_vals)}, n_background={len(background_vals)}, U={u:.1f}, p={p:.4g}",
    ]
    return "\n".join(lines)


def main():
    import os
    if os.path.exists(OUT):
        # delete OUT to force recomputing depth from the BAMs
        print(f"{OUT} already exists -- loading it rather than recomputing depth", flush=True)
        df = pd.read_csv(OUT, sep="\t")
    else:
        samples = pd.read_csv(PC1, sep="\t", nrows=0).columns[1:].tolist()
        print(f"{len(samples)} samples", flush=True)
        length = int(pd.read_csv(FAI, sep="\t", header=None, index_col=0).loc[CHROM, 1])
        print(f"{CHROM} length: {length}", flush=True)

        total = depth_track(samples, length, 0, "total")
        highqual = depth_track(samples, length, HIGHQUAL_MAPQ, "highqual")

        df = pd.concat([total, highqual], axis=1).fillna(0).reset_index().rename(columns={"index": "pos"})
        df["total_depth"] = df["total"]
        df["highqual_depth"] = df["highqual"]
        df["frac_lowmapq"] = 1 - df["highqual_depth"] / df["total_depth"].replace(0, pd.NA)
        df = df[["pos", "total_depth", "highqual_depth", "frac_lowmapq"]].sort_values("pos")
        df.to_csv(OUT, sep="\t", index=False, compression="gzip")
        print(f"Saved: {OUT} ({len(df)} windows)")

    print()
    result = test_peak_vs_background(df)
    print(result)
    with open(TEST_OUT, "w") as fh:
        fh.write(result + "\n")
    print(f"\nSaved: {TEST_OUT}")


if __name__ == "__main__":
    main()
