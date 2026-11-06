#!/usr/bin/env python3
"""
Inversion vs. low-recombination test for the INV1 core on DEINVI_00006.

If the haplotype structure is caused by an inversion, recombination is
suppressed only in heterozygotes: WITHIN each homozygous group, LD should
decay with distance roughly as it does in ordinary parts of the chromosome.
If recombination is simply low there for everyone, LD stays high within
each group too.

Groups: Strauss individuals split into 3 groups by sign-aligned mean PC1 over
the INV1 core (Strauss-only winPCA run): hom1 / het / hom2. Groups are
subsampled to the same size (N_PER_GROUP); "pooled" = random individuals from
all three groups (structure present, the usual LD situation).
LD: r^2 between imputed ALT dosages (DS) from the GWAS VCF, for SNPs with
within-group MAF >= MIN_MAF, in distance bins; median r^2 per bin.
Caveat: imputation uses haplotype information and can inflate LD in strongly
structured regions; confirm a borderline result with ngsLD on genotype
likelihoods for the same individuals.

Usage: ld_within_groups.py
"""

import subprocess

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

VCF = "/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog_MAF_auto_GLimp_clean.vcf.gz"
PC1 = "/private/groups/meyerlab/esteiger/inversions/all_localities/DEINVI_00006_GL_maf5_Strauss_flipped.pc_1.tsv.gz"
OUT = "/private/groups/meyerlab/esteiger/inversions/breakpoints/ld_within_groups.tsv"
CHROM = "DEINVI_00006"
CORE = (46_300_000, 49_100_000)
REGIONS = {"INV1 core 46.3-49.1": CORE, "control 2-8 Mb": (2_000_000, 8_000_000),
           "control 55-61 Mb": (55_000_000, 61_000_000)}
BINS = [(1, 5), (5, 10), (10, 25), (25, 50), (50, 100), (100, 300)]   # kb
N_PER_GROUP, MIN_MAF, MAX_SNPS, MAX_PAIRS = 36, 0.10, 6000, 20000
rng = np.random.default_rng(1)


def groups():
    p = pd.read_csv(PC1, sep="\t", index_col=0)
    w = p[(p.index >= CORE[0]) & (p.index <= CORE[1])].dropna(how="all")
    M = w.values; ref = M[len(M) // 2]
    sgn = np.sign([np.corrcoef(ref, m)[0, 1] for m in M])
    x = (M * sgn[:, None]).mean(0)
    km = KMeans(3, n_init=20, random_state=0).fit(x[:, None])
    order = np.argsort(km.cluster_centers_.ravel())
    name = {order[0]: "hom1", order[1]: "het", order[2]: "hom2"}
    lab = pd.Series([name[k] for k in km.labels_], index=[c.replace("_rmdup.bam", "") for c in p.columns])
    return {g: list(lab.index[lab == g]) for g in ["hom1", "het", "hom2"]}


def dosages(region, samples):
    out = subprocess.run(["bcftools", "query", "-r", f"{CHROM}:{region[0]}-{region[1]}", "-s", ",".join(samples),
                          "-f", "%POS[\t%DS]\n", VCF], capture_output=True, text=True).stdout.splitlines()
    a = np.array([l.split("\t") for l in out], dtype=float)
    return a[:, 0].astype(int), a[:, 1:]


def ld_decay(pos, D):
    f = D.mean(1) / 2
    keep = (np.minimum(f, 1 - f) >= MIN_MAF) & (D.std(1) > 0)
    pos, D = pos[keep], D[keep]
    if len(pos) > MAX_SNPS:
        idx = np.sort(rng.choice(len(pos), MAX_SNPS, replace=False)); pos, D = pos[idx], D[idx]
    Z = (D - D.mean(1, keepdims=True)) / D.std(1, keepdims=True)
    res = {}
    for lo, hi in BINS:
        i, j = [], []
        for a in range(len(pos)):
            b0, b1 = np.searchsorted(pos, [pos[a] + lo * 1000, pos[a] + hi * 1000])
            if b1 > b0:
                jj = np.arange(b0, b1); i += [a] * len(jj); j += list(jj)
        i, j = np.array(i), np.array(j)
        if len(i) > MAX_PAIRS:
            s = rng.choice(len(i), MAX_PAIRS, replace=False); i, j = i[s], j[s]
        r = (Z[i] * Z[j]).mean(1) if len(i) else np.array([np.nan])
        res[f"{lo}-{hi} kb"] = np.median(r ** 2)
    return len(pos), res


def main():
    g = groups()
    print("Strauss groups in the INV1 core:", {k: len(v) for k, v in g.items()})
    sets = {k: list(rng.choice(v, N_PER_GROUP, replace=False)) for k, v in g.items()}
    allg = g["hom1"] + g["het"] + g["hom2"]
    sets["pooled (all groups)"] = list(rng.choice(allg, N_PER_GROUP, replace=False))
    rows = []
    for rname, reg in REGIONS.items():
        pos, D = dosages(reg, allg)
        idx = {s: k for k, s in enumerate(allg)}
        for sname, ss in sets.items():
            n, res = ld_decay(pos, D[:, [idx[s] for s in ss]])
            rows.append({"region": rname, "group": sname, "n_snps": n, **{k: round(v, 3) for k, v in res.items()}})
            print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("\nmedian r^2 by distance bin (n = %d individuals per group):" % N_PER_GROUP)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
