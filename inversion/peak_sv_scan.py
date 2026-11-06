#!/usr/bin/env python3
"""
Unbiased short-read structural-variant scan across the chr 6 fog outlier peak
region, using ALL samples (no grouping used for discovery).

Step 1 (discovery): for every sample, find positions in REGION where reads
are soft-clipped (>= MIN_CLIP bp, MAPQ >= 20, primary). A candidate site is a
position (+-TOL bp) where >= MIN_IND individuals each have >= 2 clipped reads.
Step 2 (genotyping): at every candidate site, per sample, count reads clipped
there vs reads spanning it unclipped (>= 30 bp each side). An individual
carries the variant if >= 2 clipped reads and clipped fraction >= 0.15; the
clipped fraction (~0.5 heterozygote, ~1 homozygote) gives copy number.
Each site is annotated by where its clipped parts align (SA tags): elsewhere
on DEINVI_00006 within 50 kb = local rearrangement; another scaffold =
transposon-like insertion; none = unaligned clip (novel sequence/insertion).
Only AFTER discovery are carriers compared with the favorable fog haplotype
(Strauss, from fog_deinvi6_figure.R) and INV1 PC groups (breakpoint_samples.py).

Usage: peak_sv_scan.py   (runs samples in parallel with N_PROC processes)
Output: OUT_DIR/peak_sv_sites.tsv, OUT_DIR/peak_sv_genotypes.tsv
"""

import collections
import re
import subprocess
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import chisquare, spearmanr

BAM_DIR = "/private/groups/meyerlab/esteiger/mapping/merged"
OUT_DIR = "/private/groups/meyerlab/esteiger/inversions/breakpoints"
SAMPLES = f"{OUT_DIR}/breakpoint_samples.tsv"
FAV = "/private/groups/meyerlab/esteiger/selection/gemma/fog/fog_deinvi6_by_year.tsv"
META = "/private/groups/meyerlab/esteiger/inversions/metadata.tsv"
CHROM = "DEINVI_00006"
REGION = (48_950_000, 49_250_000)
PEAK = (49_031_000, 49_147_000)
MIN_CLIP, TOL, MIN_IND, N_PROC = 15, 5, 5, 16


def alen(c):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", c))


def sample_reads(s):
    out = subprocess.run(["samtools", "view", "-F", "0x904", "-q", "20", f"{BAM_DIR}/{s}",
                          f"{CHROM}:{REGION[0]}-{REGION[1]}"], capture_output=True, text=True).stdout.splitlines()
    recs = []
    for l in out:
        f = l.split("\t"); st = int(f[3]); en = st + alen(f[5])
        a = re.match(r"^(\d+)S", f[5]); z = re.search(r"(\d+)S$", f[5])
        sa = next((t[5:].split(";")[0] for t in f[11:] if t.startswith("SA:Z:")), "")
        recs.append((st, en, int(a.group(1)) if a else 0, int(z.group(1)) if z else 0, sa))
    return s, recs


def genotype(s, recs, sites):
    """0/1/2 copies (NaN if < 4 informative reads) at each site, from clipped vs spanning reads."""
    if not recs:
        return s, [np.nan] * len(sites)
    st, en, cl, cr = (np.array([r[i] for r in recs]) for i in range(4))
    clean = (cl < 10) & (cr < 10)
    out = []
    for p in sites:
        span = np.sum(clean & (st <= p - 30) & (en >= p + 35))
        clip = np.sum(((cl >= MIN_CLIP) & (st >= p - 2) & (st <= p + TOL + 2)) |
                      ((cr >= MIN_CLIP) & (en >= p - 2) & (en <= p + TOL + 2)))
        n = span + clip
        out.append(np.nan if n < 4 else (0 if clip < 2 or clip / n < 0.15 else (2 if clip / n > 0.85 else 1)))
    return s, out


def main():
    samples = sorted(pd.read_csv(SAMPLES, sep="\t")["sample"].unique())
    with Pool(N_PROC) as pool:
        data = dict(pool.map(sample_reads, samples))
    print(f"{len(samples)} samples read")

    # step 1: discovery
    ind_at = collections.defaultdict(set); sa_at = collections.defaultdict(collections.Counter)
    for s, recs in data.items():
        cnt = collections.Counter(); sas = collections.defaultdict(list)
        for st, en, cl, cr, sa in recs:
            for pos, c in ((st, cl), (en, cr)):
                if c >= MIN_CLIP and REGION[0] <= pos <= REGION[1]:
                    k = pos // TOL * TOL; cnt[k] += 1
                    if sa: sas[k].append(sa)
        for k, n in cnt.items():
            if n >= 2:
                ind_at[k].add(s)
                for sa in sas[k]:
                    x = sa.split(",")
                    if x[0] == CHROM and abs(int(x[1]) - k) <= 50_000:
                        sa_at[k]["local"] += 1
                    else:
                        sa_at[k]["other scaffold" if x[0] != CHROM else "chr6 >50kb"] += 1
    cand = sorted(k for k, v in ind_at.items() if len(v) >= MIN_IND)
    # merge neighbouring bins (same site split across a TOL boundary): keep the better-supported one
    sites = []
    for k in cand:
        if sites and k - sites[-1] <= TOL:
            if len(ind_at[k]) > len(ind_at[sites[-1]]): sites[-1] = k
        else:
            sites.append(k)
    print(f"{len(sites)} candidate sites (>= {MIN_IND} individuals with >= 2 clipped reads)")

    # step 2: genotype every site in every sample (vectorized per sample, in parallel)
    with Pool(N_PROC) as pool:
        G = dict(pool.starmap(genotype, [(s, data[s], sites) for s in samples]))
    G = pd.DataFrame(G, index=sites).T
    G.to_csv(f"{OUT_DIR}/peak_sv_genotypes.tsv", sep="\t")

    # annotate
    meta = pd.read_csv(META, sep="\t", index_col=0)
    strauss = [s for s in G.index if meta["Clusters"].get(s) == "Strauss"]
    fav = pd.read_csv(FAV, sep="\t", dtype={"SampleID": str})
    fav = fav.assign(sample=fav.SampleID + "_rmdup.bam", copies=fav.pk_score.round().clip(0, 2)).drop_duplicates("sample").set_index("sample")["copies"]
    inv1 = pd.read_csv(SAMPLES, sep="\t").query("block == 'INV1_45.5-49.1'").set_index("sample")["group"]
    rows = []
    for p in sites:
        g = G[p].dropna(); gs = G.loc[strauss, p].dropna()
        n0, n1, n2 = (gs == 0).sum(), (gs == 1).sum(), (gs == 2).sum(); N = n0 + n1 + n2
        q = (n1 + 2 * n2) / (2 * N) if N else np.nan
        hwe_p = chisquare([n0, n1, n2], [N * (1 - q) ** 2, 2 * N * q * (1 - q), N * q * q])[1] if 0 < q < 1 else np.nan
        both = gs.index.intersection(fav.index)
        rho = spearmanr(gs[both], fav[both])[0] if len(both) > 20 and gs[both].nunique() > 1 else np.nan
        kinds = sa_at[p]
        kind = ("transposon-like insertion" if kinds["other scaffold"] > kinds["local"] else
                "local rearrangement" if kinds["local"] > 0 else "unaligned clip")
        car = g[g > 0].index
        inv = {k: f"{sum(inv1.get(s) == k for s in car)}/{sum(inv1.get(s) == k for s in g.index)}" for k in ("A", "Het", "B")}
        rows.append(dict(pos=p, in_peak=PEAK[0] <= p <= PEAK[1], kind=kind, n_genotyped=len(g),
                         carriers=int((g > 0).sum()), carrier_frac=round((g > 0).mean(), 3),
                         strauss_freq=round(q, 3), strauss_hwe_p=hwe_p, rho_fav=round(rho, 2),
                         INV1_A=inv["A"], INV1_Het=inv["Het"], INV1_B=inv["B"]))
    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT_DIR}/peak_sv_sites.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
