#!/usr/bin/env python3
"""
De novo search for centromeric satellite arrays on one scaffold (default
DEINVI_00006) from high-copy k-mers.

Per WIN-bp window: fraction of positions whose canonical K-mer occurs
>= HIGH_COPY times on the scaffold, GC content and N content. For the top
windows, the most frequent k-mer and the modal spacing between its occurrences
(tandem monomer length).

Usage: centromere_scan.py [scaffold]
Output: OUT_DIR/<scaffold>_kmer_windows.tsv, printed summary
"""
import gzip, sys
import numpy as np
import pandas as pd

BASE_DIR = "/private/groups/meyerlab/esteiger"

GENOME = f"{BASE_DIR}/genome/DEINVI_genome.fasta.gz"
OUT_DIR = f"{BASE_DIR}/inversions/centromere"
CHROM = sys.argv[1] if len(sys.argv) > 1 else "DEINVI_00006"
K, WIN, HIGH_COPY, TOP = 21, 100_000, 200, 15


def read_seq():
    seq, on = [], False
    with gzip.open(GENOME, "rt") as f:
        for line in f:
            if line[0] == ">":
                if on: break
                on = line[1:].split()[0] == CHROM
            elif on:
                seq.append(line.strip())
    return "".join(seq).upper()


def kmer_codes(s):
    lut = np.full(256, 255, np.uint8)
    for i, c in enumerate("ACGT"): lut[ord(c)] = i
    b = lut[np.frombuffer(s.encode(), np.uint8)]
    bad = b == 255; b = np.where(bad, 0, b).astype(np.uint64)
    n = len(b) - K + 1
    fwd = np.zeros(n, np.uint64); rev = np.zeros(n, np.uint64)
    for j in range(K):
        fwd = (fwd << np.uint64(2)) | b[j:j + n]
        rev = rev | ((np.uint64(3) - b[j:j + n]) << np.uint64(2 * j))
    badk = np.convolve(bad.astype(np.int8), np.ones(K, np.int8), "valid") > 0
    return np.minimum(fwd, rev), badk


def main():
    s = read_seq(); L = len(s)
    print(f"{CHROM}: {L/1e6:.2f} Mb")
    code, badk = kmer_codes(s)
    u, inv, cnt = np.unique(code[~badk], return_inverse=True, return_counts=True)
    copies = np.zeros(len(code), np.int64); copies[~badk] = cnt[inv]
    hc = copies >= HIGH_COPY
    seqa = np.frombuffer(s.encode(), np.uint8)
    rows = []
    for a in range(0, L, WIN):
        w = seqa[a:a + WIN]
        nN = np.sum(w == ord("N")); acgt = len(w) - nN
        gc = np.sum((w == ord("G")) | (w == ord("C"))) / max(acgt, 1)
        rows.append(dict(start=a, end=min(a + WIN, L), high_copy_frac=hc[a:a + WIN].mean() if a < len(hc) else np.nan,
                         gc=gc, N_frac=nN / len(w)))
    df = pd.DataFrame(rows)
    # assembly gaps
    isN = (seqa == ord("N")).astype(np.int8)
    d = np.diff(np.r_[0, isN, 0]); gs, ge = np.where(d == 1)[0], np.where(d == -1)[0]
    gaps = [(a, b - a) for a, b in zip(gs, ge) if b - a >= 10]
    print(f"assembly gaps (>=10 N): {len(gaps)}; " + ", ".join(f"{a/1e6:.2f} Mb ({n} bp)" for a, n in gaps[:40]))
    # monomer length in the top windows
    top = df.sort_values("high_copy_frac", ascending=False).head(TOP).sort_values("start")
    mono = []
    for _, r in top.iterrows():
        a, b = int(r.start), int(min(r.end, len(code)))
        c = code[a:b][~badk[a:b]]
        if not len(c): mono.append(np.nan); continue
        vals, n = np.unique(c, return_counts=True); best = vals[np.argmax(n)]
        p = np.where(code[a:b] == best)[0]; sp = np.diff(p)
        mono.append(int(np.bincount(sp).argmax()) if len(sp) else np.nan)
    top["top_kmer_spacing_bp"] = mono
    df.to_csv(f"{OUT_DIR}/{CHROM}_kmer_windows.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(f"\nscaffold-wide high-copy fraction: median {df.high_copy_frac.median():.3f}, "
          f"95th pct {df.high_copy_frac.quantile(.95):.3f}")
    print(f"\ntop {TOP} windows by high-copy fraction:"); print(top.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    df["Mb"] = (df.start // 1_000_000)
    print("\nper-Mb summary:"); print(df.groupby("Mb")[["high_copy_frac", "gc", "N_frac"]].mean().round(3).to_string())


if __name__ == "__main__":
    main()
