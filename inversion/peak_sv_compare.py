#!/usr/bin/env python3
"""
Compare short-read structural-variant signals across the chr 6 fog outlier
peak between favorable-haplotype homozygotes and high-coverage non-carriers.

Region: REGION on DEINVI_00006 (fog peak 49.031-49.146 Mb plus flanks).
Signals per sample (MAPQ >= 20, primary alignments):
  1. depth in BIN bp bins, normalized to that sample's NORM_REGIONS depth;
     carriers vs the median of non-carriers (deletion ~0, duplication > 1.7)
  2. clip clusters: >= MIN_CLIP_READS reads soft-clipped (>= 20 bp) at the
     same position (+-CLIP_TOL bp) -- with where their SA parts align
  3. abnormal pairs: mate on another scaffold, same-strand mates, or mates
     > ABN_INSERT bp apart, binned per BIN
  4. reads whose mate is unmapped, per BIN
Carrier-specific features: present in both carriers, absent (<= 1 read) in
>= 80% of non-carriers.

Usage: peak_sv_compare.py
"""

import collections
import re
import subprocess

import numpy as np
import pandas as pd

BAM_DIR = "/private/groups/meyerlab/esteiger/mapping/merged"
CHROM = "DEINVI_00006"
REGION = (48_950_000, 49_250_000)
PEAK = (49_031_000, 49_147_000)
NORM_REGIONS = [(48_700_000, 48_800_000), (49_400_000, 49_500_000)]
BIN, CLIP_TOL, MIN_CLIP_READS, ABN_INSERT = 1000, 5, 3, 2000
CARRIERS = ["5059_rmdup.bam", "7991_rmdup.bam"]
N_NON = 10
FAV = "/private/groups/meyerlab/esteiger/selection/gemma/fog/fog_deinvi6_by_year.tsv"
ENV = "/private/groups/meyerlab/esteiger/env_var.txt"
OUT = "/private/groups/meyerlab/esteiger/inversions/breakpoints/peak_sv_compare"


def alen(c):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", c))


def reads(bam, a, b):
    return subprocess.run(["samtools", "view", "-F", "0x904", "-q", "20", f"{BAM_DIR}/{bam}", f"{CHROM}:{a}-{b}"],
                          capture_output=True, text=True).stdout.splitlines()


def depth_bins(bam):
    bed = f"{OUT}_{bam}.bed"
    with open(bed, "w") as f:
        for a in range(REGION[0], REGION[1], BIN):
            f.write(f"{CHROM}\t{a}\t{a + BIN}\n")
        for a, b in NORM_REGIONS:
            f.write(f"{CHROM}\t{a}\t{b}\n")
    rows = [l.split("\t") for l in subprocess.run(["samtools", "bedcov", "-Q", "20", bed, f"{BAM_DIR}/{bam}"],
                                                   capture_output=True, text=True).stdout.splitlines()]
    v = np.array([int(r[3]) / (int(r[2]) - int(r[1])) for r in rows])
    n = len(NORM_REGIONS)
    return v[:-n] / v[-n:].mean()


def signals(bam):
    clips = collections.Counter(); sa_to = collections.defaultdict(collections.Counter)
    abn = collections.Counter(); unm = collections.Counter()
    for l in reads(bam, *REGION):
        f = l.split("\t"); flag = int(f[1]); s = int(f[3]); e = s + alen(f[5])
        a = re.match(r"^(\d+)S", f[5]); z = re.search(r"(\d+)S$", f[5])
        for pos, c in ((s, int(a.group(1)) if a else 0), (e, int(z.group(1)) if z else 0)):
            if c >= 20 and REGION[0] <= pos <= REGION[1]:
                key = pos // CLIP_TOL * CLIP_TOL
                clips[key] += 1
                for t in f[11:]:
                    if t.startswith("SA:Z:"):
                        x = t[5:].split(";")[0].split(",")
                        sa_to[key][f"{x[0].replace('DEINVI_', '')}:{int(x[1]) / 1e6:.2f}"] += 1
        b = s // BIN * BIN
        if flag & 1:
            if flag & 8:
                unm[b] += 1
            elif f[6] != "=" or bool(flag & 16) == bool(flag & 32) or abs(int(f[7]) - s) > ABN_INSERT:
                abn[b] += 1
    return clips, sa_to, abn, unm


def main():
    g = pd.read_csv(FAV, sep="\t", dtype={"SampleID": str})
    g["copies"] = g.pk_score.round().clip(0, 2)
    env = pd.read_csv(ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    g["cov"] = pd.to_numeric(env["Cov"], errors="coerce").reindex(g.SampleID).values
    non = (g[g.copies == 0].sort_values("cov", ascending=False).head(N_NON).SampleID + "_rmdup.bam").tolist()
    print(f"carriers: {CARRIERS}\nnon-carriers ({len(non)}): {non}\n")

    # 1. depth
    dep = {s: depth_bins(s) for s in CARRIERS + non}
    non_med = np.median(np.array([dep[s] for s in non]), axis=0)
    starts = np.arange(REGION[0], REGION[1], BIN)
    print("1. DEPTH: 1 kb bins where BOTH carriers deviate from the non-carrier median "
          "(ratio < 0.4 or > 1.7, non-carrier depth >= 0.3):")
    hits = 0
    for i, a in enumerate(starts):
        r = [dep[s][i] / non_med[i] if non_med[i] > 0 else np.nan for s in CARRIERS]
        if non_med[i] >= 0.3 and (all(x < 0.4 for x in r) or all(x > 1.7 for x in r)):
            hits += 1
            print(f"   {a/1e6:.3f}-{(a+BIN)/1e6:.3f} Mb: non-carriers {non_med[i]:.2f}, carriers "
                  + ", ".join(f"{dep[s][i]:.2f}" for s in CARRIERS) + (" (in peak)" if PEAK[0] <= a <= PEAK[1] else ""))
    if not hits:
        print("   none")

    # 2-4. read signals
    sig = {s: signals(s) for s in CARRIERS + non}
    print("\n2. CLIP CLUSTERS present in both carriers (>= %d reads each), absent (<=1) in >= 80%% of non-carriers:" % MIN_CLIP_READS)
    keys = set.intersection(*[{k for k, v in sig[s][0].items() if v >= MIN_CLIP_READS} for s in CARRIERS])
    found = 0
    for k in sorted(keys):
        absent = sum(sig[s][0].get(k, 0) <= 1 for s in non)
        if absent >= 0.8 * len(non):
            found += 1
            sa = collections.Counter()
            for s in CARRIERS:
                sa.update(sig[s][1][k])
            print(f"   {k/1e6:.6f} Mb: carriers " + ", ".join(str(sig[s][0][k]) for s in CARRIERS)
                  + f" reads; non-carriers with >1 read: {len(non) - absent}; SA parts -> {sa.most_common(3) or 'none (clip unaligned)'}"
                  + (" (in peak)" if PEAK[0] <= k <= PEAK[1] else ""))
    if not found:
        print("   none")

    for idx, name in [(2, "3. ABNORMAL PAIRS"), (3, "4. MATE-UNMAPPED READS")]:
        print(f"\n{name}: 1 kb bins where both carriers exceed 3x the non-carrier median and have >= 4 reads:")
        found = 0
        for a in starts:
            nm = np.median([sig[s][idx].get(a, 0) for s in non])
            c = [sig[s][idx].get(a, 0) for s in CARRIERS]
            if all(x >= 4 and x > 3 * max(nm, 1) for x in c):
                found += 1
                print(f"   {a/1e6:.3f} Mb: carriers {c}, non-carrier median {nm:.0f}" + (" (in peak)" if PEAK[0] <= a <= PEAK[1] else ""))
        if not found:
            print("   none")


if __name__ == "__main__":
    main()
