#!/usr/bin/env python3
"""
Sliding-window "ambiguous mapping" check along DEINVI_00006: an independent,
filter-free test of whether the fog GWAS peak (GNS1/SNF1, ~49.03-49.15 Mb)
could be a collapsed-paralog/mismapping artifact, since it's centromere-
adjacent (centromere: 45.0-48.5 Mb) and extra scrutiny is warranted there.

Why not just reuse the paralog-filtered VCF's depth (as an earlier version
of this analysis did, in DEINVI_0006_coverage_windows.tsv.gz / panel D of
chr6_synteny_panels.py -- removed 2026-10-23, per user decision, since it
wasn't actually useful for this question): that VCF's sites already passed
a paralog/mismapping filter BEFORE this depth could be computed, so a truly
mismapped region would just lose sites there, not show elevated depth among
the sites that survive -- the check was circular for exactly the question
being asked here.

This instead computes two filter-INDEPENDENT depth tracks directly from the
merged BAMs for all 339 all-locality samples (the same sample set used in
chr6_synteny_panels.py's PC1/heterozygosity panels), summed across samples,
on the same 10kb grid:
  total_depth    -- samtools depth -Q 0  (every read, any mapping quality)
  highqual_depth -- samtools depth -Q 20 (only MAPQ >= 20 reads)
  frac_lowmapq   = 1 - highqual_depth / total_depth
A region where reads pile up ambiguously (e.g. a collapsed duplication) will
show elevated frac_lowmapq, since multi-mapping reads get MAPQ ~0 regardless
of local depth -- this is the standard signature checked for this purpose,
and the earlier (removed) coverage-only check couldn't have caught it even
in principle.

Parallelized across N_CHUNKS scaffold chunks per depth pass (matches the
approach benchmarked for the earlier BAM-based coverage computation, ~45 min
single-threaded for one pass across all 339 BAMs x 70.5 Mb); two passes here
(-Q 0 and -Q 20) roughly double that.

Statistical test (2026-10-23, per user decision -- this isn't going in a
figure, so the comparison needs to be a real test, not an eyeballed panel):
a one-sided Mann-Whitney U test of frac_lowmapq in the peak's windows against
all other genome-wide windows (centromere windows also excluded from the
background, so a globally elevated centromere can't itself drive the
result). One-sided, alternative="greater", since the only direction that
would actually support a mismapping concern is the peak showing HIGHER
ambiguous-mapping than background -- testing both directions here would just
be diluting a pre-specified, directional question.

Usage: compute_deinvi6_mapq.py
Output: OUT (pos, total_depth, highqual_depth, frac_lowmapq), tab-separated,
        gzip-compressed; TEST_OUT (the Mann-Whitney result, plain text, for
        the archive); both also printed to stdout.
"""
from concurrent.futures import ThreadPoolExecutor
import subprocess

import pandas as pd
from scipy.stats import mannwhitneyu

BAM_DIR = "/private/groups/meyerlab/esteiger/mapping/merged"
PC1 = "/private/groups/meyerlab/esteiger/inversions/all_localities/DEINVI_00006_GL_maf5_flipped.pc_1.tsv.gz"
FAI = "/private/groups/meyerlab/esteiger/inversions/filtered_scaffolds.fai"
OUT = "/private/groups/meyerlab/esteiger/inversions/all_localities/DEINVI_00006_mapq_windows.tsv.gz"
TEST_OUT = "/private/groups/meyerlab/esteiger/inversions/all_localities/DEINVI_00006_mapq_peak_test.txt"
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
    """One-sided Mann-Whitney U: is the peak's frac_lowmapq GREATER than background?
    Background excludes both the peak and the centromere, so an elevated
    centromere baseline can't inflate the comparison either way."""
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
        # Re-running just to redo the statistical test (e.g. after editing it)
        # shouldn't repeat the ~15-20 min BAM computation -- load the window
        # data straight from disk if it's already there. Delete OUT first to
        # force a genuine recompute (e.g. if the BAMs or sample list changed).
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
