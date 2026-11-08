#!/usr/bin/env python3
"""
Long-read check of the INV1 region (DEINVI_00006 ~42-51 Mb) in the DEINVI
reference individual, using the chr-6 BAMs written by longread_chr6.sh.

1. Junctions: opposite-strand split reads (both parts >= MIN_PART bp, MAPQ
   >= MIN_MAPQ, > MIN_DIST apart, >= 1 end in REGION), clustered at CLUSTER_BP.
2. Breakpoint scan: every SCAN_STEP bp across REGION, reads that break there
   vs. reads that span it, flagged against the 99th percentile of BASELINE_N
   random positions outside REGION.

Usage: longread_inv1.py <bam> [<bam> ...]
"""

import re
import subprocess
import sys
from collections import Counter

import numpy as np

CHROM = "DEINVI_00006"
REGION = (42_000_000, 51_000_000)
MIN_PART, MIN_MAPQ, MIN_DIST = 5000, 50, 20000
CLUSTER_BP, MIN_READS = 2000, 3
SCAN_STEP, BASELINE_N = 5000, 100
CANDIDATES = {45.30: "45.3 (left edge cand.)", 46.12: "46.12 (short-read A)", 47.80: "47.80 (short-read A)",
              48.87: "48.85-48.90 (short-read B)", 49.30: "49.3 (right edge cand.)"}


def alen(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cigar))


def clips(cigar):
    a = re.match(r"^(\d+)[SH]", cigar)
    b = re.search(r"(\d+)[SH]$", cigar)
    return (int(a.group(1)) if a else 0, int(b.group(1)) if b else 0)


def junction_coord(pos, cigar):
    """Alignment end next to the larger clip (where the read leaves this locus)."""
    left, right = clips(cigar)
    return pos + alen(cigar) if right >= left else pos


def view(bam, region=None, q=20):
    cmd = ["samtools", "view", "-F", "0x904", "-q", str(q), bam] + ([region] if region else [CHROM])
    return subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()


def junctions(bam):
    jc = []
    for line in view(bam, q=MIN_MAPQ):
        f = line.split("\t")
        sa = [t for t in f[11:] if t.startswith("SA:Z:")]
        if not sa or alen(f[5]) < MIN_PART:
            continue
        strand = "-" if int(f[1]) & 16 else "+"
        j1 = junction_coord(int(f[3]), f[5])
        for part in sa[0][5:].split(";"):
            if not part:
                continue
            rn, sp, ss, sc, sq, _ = part.split(",")
            if rn != CHROM or int(sq) < MIN_MAPQ or alen(sc) < MIN_PART or ss == strand:
                continue
            j2 = junction_coord(int(sp), sc)
            if abs(j2 - j1) < MIN_DIST:
                continue
            a, b = sorted((j1, j2))
            if REGION[0] <= a <= REGION[1] or REGION[0] <= b <= REGION[1]:
                jc.append((a, b))
    cl = Counter((round(a / CLUSTER_BP), round(b / CLUSTER_BP)) for a, b in jc)
    return len(jc), sorted(((a * CLUSTER_BP, b * CLUSTER_BP), n) for (a, b), n in cl.items() if n >= MIN_READS)


def break_span(bam, p):
    brk = span = 0
    for line in view(bam, f"{CHROM}:{p - 40000}-{p + 40000}"):
        f = line.split("\t")
        s = int(f[3]); e = s + alen(f[5]); cl, cr = clips(f[5])
        if s <= p - 5000 and e >= p + 5000:
            span += 1
        elif (cr >= 1000 and abs(e - p) <= 1000) or (cl >= 1000 and abs(s - p) <= 1000):
            brk += 1
    return brk, span


def main(bams):
    rng = np.random.default_rng(1)
    outside = np.concatenate([rng.integers(1_000_000, 12_000_000, BASELINE_N // 2),
                              rng.integers(53_000_000, 69_000_000, BASELINE_N // 2)])
    for bam in bams:
        print(f"\n=================== {bam}")
        n, cl = junctions(bam)
        print(f"strict opposite-strand split reads with an end in {REGION[0]/1e6:.0f}-{REGION[1]/1e6:.0f} Mb: {n}")
        print(f"junction clusters (>= {MIN_READS} reads, {CLUSTER_BP/1e3:.0f} kb resolution):")
        for (a, b), k in cl:
            print(f"   {a/1e6:8.3f} <-> {b/1e6:8.3f} Mb   {k} reads")
        base = np.array([break_span(bam, int(p)) for p in outside])
        thr = np.percentile(base[:, 0], 99)
        print(f"baseline at {BASELINE_N} random positions outside the region: breaking reads median "
              f"{np.median(base[:,0]):.0f}, 99th pct {thr:.0f}; spanning reads median {np.median(base[:,1]):.0f}")
        print(f"positions in {REGION[0]/1e6:.0f}-{REGION[1]/1e6:.0f} Mb with breaking reads > 99th pct of baseline:")
        hits = []
        for p in range(REGION[0], REGION[1] + 1, SCAN_STEP):
            b, s = break_span(bam, p)
            if b > thr:
                hits.append((p, b, s))
        for p, b, s in hits:
            print(f"   {p/1e6:8.3f} Mb: {b:3d} breaking, {s:3d} spanning ({100*b/max(b+s,1):.0f}% break)")
        print("candidate positions:")
        for c, lab in CANDIDATES.items():
            b, s = break_span(bam, int(c * 1e6))
            print(f"   {lab:28s}: {b:3d} breaking, {s:3d} spanning")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
