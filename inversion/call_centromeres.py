#!/usr/bin/env python3
"""
Genome-wide centromere calls from the centromere-associated LINE family.

In sunflower, centromeres carry a LINE-like repeat found on almost all
chromosomes. In Deinandra, LINE family CEN_FAMILY (RepeatModeler name) forms
one tight cluster on nearly every chromosome-scale scaffold, and on
DEINVI_00006 that cluster is the 45.5-48.0 Mb LINE-rich region.

Per scaffold >= MIN_LEN: CEN_FAMILY bp per WIN window, smoothed over SMOOTH
windows; centromere = the contiguous run of windows around the densest one
whose smoothed density is >= FRAC_OF_PEAK of the peak. A call is "confident"
if >= MIN_CONC of the scaffold's CEN_FAMILY bp lie within +-FLANK of the
call and the peak window holds >= MIN_PEAK_BP. Arm lengths are measured from
the call's edges to the scaffold ends.

Usage: call_centromeres.py
Output: OUT (one row per scaffold), OUT_WIN (window densities)
"""

import subprocess

import numpy as np
import pandas as pd

RM_OUT = "/private/groups/meyerlab/esteiger/inversions/gene_density/DEINVI_genome.fasta.LTRfix.out"
FAI = "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.fai.gz"
OUT = "/private/groups/meyerlab/esteiger/inversions/centromere/centromeres.tsv"
OUT_WIN = "/private/groups/meyerlab/esteiger/inversions/centromere/cen_family_windows.tsv"
CEN_FAMILY = "rnd-1_family-362"
MIN_LEN, WIN, SMOOTH = 10e6, 500_000, 3
FRAC_OF_PEAK, FLANK, MIN_CONC, MIN_PEAK_BP = 0.2, 2e6, 0.4, 5_000


def main():
    fai = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chr", "len"], compression="gzip")
    fai = fai[fai.len >= MIN_LEN]
    out = subprocess.run(["awk", f'$10=="{CEN_FAMILY}" {{print $5"\\t"$6"\\t"$7}}', RM_OUT],
                         capture_output=True, text=True).stdout
    rm = pd.DataFrame([l.split("\t") for l in out.splitlines()], columns=["chr", "s", "e"])
    rm[["s", "e"]] = rm[["s", "e"]].astype(int)
    rows, wins = [], []
    for chrom, L in zip(fai.chr, fai.len):
        r = rm[rm.chr == chrom]
        nwin = int(np.ceil(L / WIN))
        bp = np.bincount((r.s // WIN).values, weights=(r.e - r.s + 1).values, minlength=nwin)[:nwin]
        sm = np.convolve(bp, np.ones(SMOOTH) / SMOOTH, mode="same")
        wins += [dict(chr=chrom, start=i * WIN, cen_bp=bp[i], smoothed=sm[i]) for i in range(nwin)]
        if bp.sum() == 0:
            rows.append(dict(chr=chrom, len=L, cen_start=np.nan, cen_end=np.nan, cen_mid=np.nan,
                             concentration=np.nan, peak_bp=0, confident=False))
            continue
        pk = int(np.argmax(sm)); thr = FRAC_OF_PEAK * sm[pk]
        lo = pk
        while lo > 0 and sm[lo - 1] >= thr: lo -= 1
        hi = pk
        while hi < nwin - 1 and sm[hi + 1] >= thr: hi += 1
        cs, ce = lo * WIN, min((hi + 1) * WIN, L)
        near = r[(r.s >= cs - FLANK) & (r.s <= ce + FLANK)]
        conc = (near.e - near.s + 1).sum() / (r.e - r.s + 1).sum()
        rows.append(dict(chr=chrom, len=L, cen_start=cs, cen_end=ce, cen_mid=(cs + ce) / 2,
                         concentration=round(conc, 3), peak_bp=int(bp.max()),
                         confident=bool(conc >= MIN_CONC and bp.max() >= MIN_PEAK_BP)))
    res = pd.DataFrame(rows)
    res["left_arm"] = res.cen_start
    res["right_arm"] = res.len - res.cen_end
    res.to_csv(OUT, sep="\t", index=False)
    pd.DataFrame(wins).to_csv(OUT_WIN, sep="\t", index=False)
    pd.set_option("display.width", 200)
    show = res.assign(len_Mb=res.len / 1e6, cen_Mb=[f"{a/1e6:.1f}-{b/1e6:.1f}" if a == a else "-" for a, b in
                                                     zip(res.cen_start, res.cen_end)])
    print(show[["chr", "len_Mb", "cen_Mb", "concentration", "peak_bp", "confident"]].to_string(
        index=False, float_format=lambda x: f"{x:.1f}"))
    print(f"\n{res.confident.sum()} of {len(res)} scaffolds >= {MIN_LEN/1e6:.0f} Mb with a confident call")


if __name__ == "__main__":
    main()
