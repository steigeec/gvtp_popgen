#!/usr/bin/env python3
"""
Centromeric haplotype frequencies and FST at the chromosome 6 centromere.

1. Haplotype frequencies per population from breakpoint_samples.tsv (block
   INV1_45.5-49.1: A / Het / B), freq(A) = (2 n_A + n_Het) / 2n, with the
   COPG2-rearrangement sub-lineage of B (j1_genotypes.tsv) shown within B.
   Untrustworthy-coordinate plants excluded; CDFW and CDFW_MR pooled (Gaviota).
2. FST at the chr 6 centromere call (centromeres.tsv) for each population pair
   with ANGSD window-FST tracks (selection/window_Fst/<a>.<b>_window), compared
   with all genome-wide windows and with the other chromosomes' centromere
   calls (median FST of windows with midpoint in the call). Windows with
   < MIN_SITES sites are dropped.

Usage: centromere_haplotype_fst.py
Output: OUT_FIG, OUT_FREQ, OUT_FST
"""

import itertools
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, f"{BASE_DIR}/gvtp_popgen/selection/gemma")
import temporal_null as tn

BASE_DIR = "/private/groups/meyerlab/esteiger"

BP = f"{BASE_DIR}/inversions/breakpoints"
FIELD = f"{BASE_DIR}/field_data.txt"
FST_DIR = f"{BASE_DIR}/selection/window_Fst"
OUT_DIR = f"{BASE_DIR}/inversions/centromere"
OUT_FIG, OUT_FREQ, OUT_FST = (f"{OUT_DIR}/centromere_haplotype_fst.jpg", f"{OUT_DIR}/centromere_haplotype_freq.tsv",
                              f"{OUT_DIR}/centromere_fst.tsv")
CHROM, PEAK = "DEINVI_00006", (49_031_257, 49_146_391)
POP_FILE = {"strauss": "Tranquillon Mtn", "cdfw": "Gaviota", "conception": "Pt. Conception",
            "lionhead": "Lion's Head", "arguello": "Pt. Arguello"}
LOC = {"Strauss": "Tranquillon Mtn", "CDFW": "Gaviota", "CDFW_MR": "Gaviota", "Point_Conception": "Pt. Conception",
       "VAFB_Lion_Head": "Lion's Head", "VAFB_Arguello": "Pt. Arguello", "Dangermond": "Dangermond"}
MIN_SITES = 100
LABEL_SIZE, TICK_SIZE = 18, 14


def read_fst(a, b):
    for x, y in [(a, b), (b, a)]:
        try:
            # header names 4 columns but rows have 5 (FST unnamed), so read without it
            d = pd.read_csv(f"{FST_DIR}/{x}.{y}_window", sep="\t", header=None, skiprows=1,
                            names=["region", "chr", "midPos", "Nsites", "Fst"])
            return d[d.Nsites >= MIN_SITES]
        except FileNotFoundError:
            continue
    return None


def main():
    # ---- 1. haplotype frequencies ----
    g = pd.read_csv(f"{BP}/breakpoint_samples.tsv", sep="\t").query("block == 'INV1_45.5-49.1'")
    g["SampleID"] = g["sample"].str.replace("_rmdup.bam", "", regex=False)
    f = pd.read_csv(FIELD, sep="\t", dtype={"SampleID": str}).set_index("SampleID")
    g = g[~g.SampleID.isin(tn.BAD_COORD_IDS) & g.SampleID.isin(f.index)]
    g["pop"] = f.loc[g.SampleID, "Location"].map(LOC).values
    j = pd.read_csv(f"{BP}/j1_genotypes.tsv", sep="\t")
    j["SampleID"] = j["sample"].str.replace("_rmdup.bam", "", regex=False)
    j["c_copies"] = j.geno.map({"hom ref": 0, "het": 1, "hom event": 2})
    g = g.merge(j[["SampleID", "c_copies"]], on="SampleID", how="left")
    rows = []
    for pop, x in g.groupby("pop"):
        n = len(x); nA, nH, nB = (x.group == "A").sum(), (x.group == "Het").sum(), (x.group == "B").sum()
        c = x.c_copies.dropna()
        rows.append(dict(population=pop, n=n, n_A=nA, n_Het=nH, n_B=nB, freq_A=(2 * nA + nH) / (2 * n),
                         freq_B=(2 * nB + nH) / (2 * n), freq_COPG2_C=c.sum() / (2 * len(c)) if len(c) else np.nan,
                         n_COPG2_genotyped=len(c)))
    fr = pd.DataFrame(rows).set_index("population").reindex([p for p in list(POP_FILE.values()) + ["Dangermond"]
                                                              if p in set(g["pop"])])
    fr.to_csv(OUT_FREQ, sep="\t")
    print(fr.round(3).to_string())

    # ---- 2. FST at centromeres ----
    cen = tn.centromeres()
    frows, plot = [], {}
    for a, b in itertools.combinations(POP_FILE, 2):
        d = read_fst(a, b)
        if d is None:
            print(f"no FST file for {a}-{b}"); continue
        allw = d.Fst.values
        cm = {}
        for chrom, c in cen.iterrows():
            w = d[(d.chr == chrom) & d.midPos.between(c.cen_start, c.cen_end)]
            if len(w) >= 5:
                cm[chrom] = w.Fst.median()
        pk = d[(d.chr == CHROM) & d.midPos.between(*PEAK)].Fst.median()
        c6 = cm.pop(CHROM, np.nan); other = np.array(list(cm.values()))
        lbl = f"{POP_FILE[a]} - {POP_FILE[b]}"
        frows.append(dict(pair=lbl, genome_median=np.median(allw), chr6_centromere=c6,
                          chr6_cen_genome_pct=(allw < c6).mean() * 100,
                          other_centromeres_median=np.median(other), chr6_cen_rank_among_centromeres=f"{int((other < c6).sum()) + 1} of {len(other) + 1}",
                          fog_peak=pk, fog_peak_genome_pct=(allw < pk).mean() * 100))
        plot[lbl] = dict(all=allw, other=other, c6=c6, pk=pk)
    fst = pd.DataFrame(frows)
    fst.to_csv(OUT_FST, sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(fst.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # ---- figure ----
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(18, 6.5), gridspec_kw=dict(width_ratios=[1, 2.4], wspace=0.55))
    fr = fr[fr.n >= 5]                                    # Dangermond (1 plant) not plotted; kept in OUT_FREQ
    pops = list(fr.index); x = np.arange(len(pops))
    a1.bar(x, fr.freq_A, color="#4d4d4d", label="A")
    bc = (fr.freq_B - fr.freq_COPG2_C.fillna(0)).clip(lower=0)
    a1.bar(x, bc, bottom=fr.freq_A, color="#bdbdbd", label="B")
    a1.bar(x, fr.freq_COPG2_C.fillna(0), bottom=fr.freq_A + bc, color="white", edgecolor="black", hatch="//",
           label="B, COPG2 sub-lineage")
    for i, n in enumerate(fr.n):
        a1.text(i, 1.01, f"n={n}", ha="center", va="bottom", fontsize=TICK_SIZE - 3)
    a1.set_xticks(x, pops, rotation=40, ha="right", fontsize=TICK_SIZE - 1)
    a1.set_ylim(0, 1.1); a1.set_ylabel("Centromeric haplotype frequency", fontsize=LABEL_SIZE)
    a1.legend(frameon=False, fontsize=TICK_SIZE - 3, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    pairs = list(plot); xx = np.arange(len(pairs))
    a2.boxplot([plot[p]["all"] for p in pairs], positions=xx, widths=0.5, showfliers=False,
               boxprops=dict(color="grey"), medianprops=dict(color="grey"), whiskerprops=dict(color="grey"),
               capprops=dict(color="grey"))
    rng = np.random.default_rng(1)
    for i, p in enumerate(pairs):
        o = plot[p]["other"]
        a2.scatter(i + rng.uniform(-0.12, 0.12, len(o)), o, s=30, facecolor="white", edgecolor="black", zorder=3,
                   label="other centromeres" if i == 0 else None)
        a2.scatter(i, plot[p]["c6"], s=200, marker="D", color="black", zorder=4, label="chr 6 centromere" if i == 0 else None)
        a2.scatter(i, plot[p]["pk"], s=200, marker="^", facecolor="#d55e00", edgecolor="black", zorder=4,
                   label="chr 6 fog peak" if i == 0 else None)
    short = {"Tranquillon Mtn": "TM", "Gaviota": "Gav", "Pt. Conception": "PtC", "Lion's Head": "LH", "Pt. Arguello": "PtA"}
    a2.set_xticks(xx, ["-".join(short[q] for q in p.split(" - ")) for p in pairs], fontsize=TICK_SIZE - 1, rotation=35, ha="right")
    a2.set_ylabel("FST", fontsize=LABEL_SIZE)
    a2.legend(frameon=False, fontsize=TICK_SIZE - 3, loc="upper left")
    for ax in (a1, a2):
        ax.tick_params(axis="y", labelsize=TICK_SIZE); ax.spines[["top", "right"]].set_visible(False)
    for ax, l in [(a1, "A"), (a2, "B")]:
        ax.text(-0.12 if l == "A" else -0.06, 1.04, l, transform=ax.transAxes, fontsize=24, va="bottom")
    fig.savefig(OUT_FIG, dpi=180, bbox_inches="tight")
    print(f"Saved: {OUT_FIG}")


if __name__ == "__main__":
    main()
