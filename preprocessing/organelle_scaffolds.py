#!/usr/bin/env python3
"""
Classify scaffolds as plastid or mitochondrial from NUCmer alignments to
published organelle genomes (organelle_nucmer.sh) and read depth.

A scaffold is organellar when alignments of >= MIN_ALN_LEN bp to that
organelle's genomes cover >= MIN_COV of its length and its depth is >= MIN_DEPTH
times nuclear depth. Depth = mapped reads / length (samtools idxstats); nuclear
depth = the same over all scaffolds >= 1 Mb; the ratio is computed per sample
and its median across samples used. Borderline scaffolds are assigned manually
(MANUAL_*).

Input:  BASE_DIR/genome/organelle_refs/organelle_vs_genome.coords
        BASE_DIR/genome/organelle_refs/idxstats/*.idxstats
          (samtools idxstats of each merged, deduplicated BAM)
Output: BASE_DIR/genome/chloroplast.txt, BASE_DIR/genome/mito.txt
Conda env: mummer (samtools for the idxstats files; this script needs only Python 3)
  conda create -n mummer -c conda-forge -c bioconda mummer4=4.0.0 samtools
Usage:  python3 organelle_scaffolds.py
"""
import collections
import glob
import statistics

BASE_DIR = "/private/groups/meyerlab/esteiger"
REF_DIR = f"{BASE_DIR}/genome/organelle_refs"
COORDS = f"{REF_DIR}/organelle_vs_genome.coords"
IDXSTATS = f"{REF_DIR}/idxstats/*.idxstats"
OUT = {"plastid": f"{BASE_DIR}/genome/chloroplast.txt", "mito": f"{BASE_DIR}/genome/mito.txt"}

MIN_ALN_LEN = 500
MIN_COV = 0.2
MIN_DEPTH = 6
NUCLEAR_MIN_LEN = 1_000_000   # chromosome-scale scaffolds define nuclear depth

# Borderline scaffolds assigned by inspecting alignment identity and depth
MANUAL_PLASTID = ["DEINVI_03031", "DEINVI_10394"]
MANUAL_NUCLEAR = ["DEINVI_05241", "DEINVI_05983", "DEINVI_08486",
                  "DEINVI_10033", "DEINVI_12868", "DEINVI_12981", "DEINVI_26210"]


def covered(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_s is None or s > cur_e + 1:
            if cur_s is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    return total + (cur_e - cur_s + 1 if cur_s is not None else 0)


def main():
    # show-coords -rclTH columns: S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ COVR COVQ REF QUERY
    aln = collections.defaultdict(lambda: collections.defaultdict(list))
    qlen = {}
    for line in open(COORDS):
        c = line.rstrip("\n").split("\t")
        if int(c[5]) < MIN_ALN_LEN:
            continue
        s, e = sorted((int(c[2]), int(c[3])))
        aln[c[12]][c[11].split("|")[0]].append((s, e))
        qlen[c[12]] = int(c[8])

    rel = collections.defaultdict(list)
    for f in glob.glob(IDXSTATS):
        depth, nuc_reads, nuc_len = {}, 0, 0
        for line in open(f):
            name, length, mapped, _ = line.split("\t")
            if name == "*" or int(length) == 0:
                continue
            depth[name] = int(mapped) / int(length)
            if int(length) >= NUCLEAR_MIN_LEN:
                nuc_reads += int(mapped)
                nuc_len += int(length)
        if nuc_reads == 0:
            continue
        nuclear = nuc_reads / nuc_len
        for k in qlen:
            rel[k].append(depth.get(k, 0) / nuclear)
    rel_depth = {k: statistics.median(v) for k, v in rel.items()}

    calls = {"plastid": set(), "mito": set()}
    for q in qlen:
        cov = {org: covered(aln[q][org]) / qlen[q] for org in ("plastid", "mito")}
        org = max(cov, key=cov.get)
        if cov[org] >= MIN_COV and rel_depth[q] >= MIN_DEPTH:
            calls[org].add(q)
    calls["plastid"] = (calls["plastid"] | set(MANUAL_PLASTID)) - set(MANUAL_NUCLEAR)

    for org, names in calls.items():
        with open(OUT[org], "w") as fh:
            fh.write("\n".join(sorted(names)) + "\n")
        print(f"{org}: {len(names)} scaffolds, median relative depth "
              f"{statistics.median(rel_depth[q] for q in names):.1f}x nuclear -> {OUT[org]}")


if __name__ == "__main__":
    main()
