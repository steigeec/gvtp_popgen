#!/usr/bin/env python3
"""
Summarize read-level inversion-breakpoint evidence on DEINVI_00006 from
breakpoint_reads.sh output, contrasting the PC1/PC2 groups chosen by
breakpoint_samples.py.

Evidence (MAPQ >= 20, ends > 20 kb apart): PAIR = same-strand read pairs,
SPLIT = split reads with parts on opposite strands. Both ends are binned (BIN
bp), and for each linked pair of bins the number of individuals per group with
>= 1 supporting read is counted.

Outputs (OUT_DIR):
  links_<block>.tsv     links supported by >= MIN_IND individuals of any group
  group_specific.tsv    links where one homozygous group carries support more
                        often than the other (one-sided Fisher's exact test,
                        Bonferroni over links x 2 directions per block)
  links_<block>.jpg     link positions (end 1 vs end 2) per group
"""

import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = "/private/groups/meyerlab/esteiger"

OUT_DIR = f"{BASE_DIR}/inversions/breakpoints"
BIN = 50_000
MIN_IND = 3                 # min. supporting individuals in some group
ALPHA = 0.01                # Bonferroni-corrected significance for group-specific links
MIN_FRAC = 0.2              # ... and >= this share of the supporting group's individuals
CANDIDATES = [15.8, 21.2, 24.1, 27.4, 29.0, 30.4, 37.5, 45.3, 49.3]   # Mb


def load_sample(sample):
    f = f"{OUT_DIR}/reads/{sample}.tsv"
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        return pd.DataFrame(columns=["type", "lo", "hi", "orient"])
    d = pd.read_csv(f, sep="\t", header=None,
                    names=["type", "pos", "strand", "opos", "ostrand", "qname"])
    pairs = d[d.type == "PAIR"]
    pairs = pairs[pairs.groupby("qname").qname.transform("size") == 2]   # both mates MAPQ >= 20
    pairs = pairs.drop_duplicates("qname")
    splits = d[d.type == "SPLIT"].drop_duplicates(["qname", "opos"])
    ev = pd.concat([pairs, splits])
    ev["lo"] = ev[["pos", "opos"]].min(axis=1) // BIN * BIN
    ev["hi"] = ev[["pos", "opos"]].max(axis=1) // BIN * BIN
    ev["orient"] = ev.strand + ev.ostrand
    return ev[["type", "lo", "hi", "orient"]]


def main():
    sel = pd.read_csv(f"{OUT_DIR}/breakpoint_samples.tsv", sep="\t")
    ev = {s: load_sample(s) for s in sel["sample"].unique()}
    missing = [s for s, e in ev.items() if not os.path.exists(f"{OUT_DIR}/reads/{s}.tsv")]
    if missing:
        print(f"WARNING: {len(missing)} samples have no output yet")

    spec_rows = []
    for block, bsel in sel.groupby("block"):
        groups = sorted(bsel.group.unique())
        per = []
        for _, r in bsel.iterrows():
            e = ev[r["sample"]]
            if len(e) == 0:
                continue
            c = e.groupby(["lo", "hi"]).agg(reads=("type", "size"),
                                            split=("type", lambda t: (t == "SPLIT").sum()))
            c["group"], c["sample"] = r.group, r["sample"]
            per.append(c.reset_index())
        per = pd.concat(per)
        n_ind = per.pivot_table(index=["lo", "hi"], columns="group", values="sample",
                                aggfunc="nunique", fill_value=0)
        reads = per.pivot_table(index=["lo", "hi"], columns="group", values="reads",
                                aggfunc="sum", fill_value=0).add_suffix("_reads")
        split = per.groupby(["lo", "hi"]).split.sum().rename("split_reads")
        links = n_ind.join(reads).join(split)
        for g in ["A", "Het", "B"]:
            if g not in links:
                links[g] = 0
        links = links[links[groups].max(axis=1) >= MIN_IND].reset_index()
        links["lo_Mb"], links["hi_Mb"] = links.lo / 1e6, links.hi / 1e6
        links["span_Mb"] = links.hi_Mb - links.lo_Mb
        links.to_csv(f"{OUT_DIR}/links_{block}.tsv", sep="\t", index=False)

        # group-specific links
        n_grp = bsel.group.value_counts()
        for g in ["A", "Het", "B"]:
            links[f"{g}_frac"] = links[g] / n_grp.get(g, np.nan)
        n_tests = 2 * len(links)
        for g, other in [("A", "B"), ("B", "A")]:
            if g not in groups or other not in groups:
                continue
            ng, no = n_grp[g], n_grp[other]
            pv = [fisher_exact([[a, ng - a], [b, no - b]], alternative="greater")[1]
                  for a, b in zip(links[g], links[other])]
            sig = links.assign(p_bonf=np.minimum(np.array(pv) * n_tests, 1))
            sig = sig[(sig.p_bonf < ALPHA) & (sig[f"{g}_frac"] >= MIN_FRAC)].copy()
            sig["block"], sig["supporting_group"] = block, g
            sig["n_A"], sig["n_Het"], sig["n_B"] = n_grp.get("A", 0), n_grp.get("Het", 0), n_grp.get("B", 0)
            spec_rows.append(sig)

        # plot: end 1 vs end 2, one panel per group
        fig, axs = plt.subplots(1, len(groups), figsize=(5.5 * len(groups), 5.2), sharex=True, sharey=True)
        for ax, g in zip(np.atleast_1d(axs), groups):
            l = links[links[g] >= MIN_IND]
            ax.scatter(l.lo_Mb, l.hi_Mb, s=l[g] * 6, c="black", alpha=0.5, edgecolor="none")
            for c in CANDIDATES:
                ax.axvline(c, color="tab:red", lw=0.5, alpha=0.5)
                ax.axhline(c, color="tab:red", lw=0.5, alpha=0.5)
            ax.set_title(f"{block}: group {g} (n={int((bsel.group == g).sum())})")
            ax.set_xlabel("End 1 (Mb)")
            ax.set_xlim(9, 55); ax.set_ylim(9, 55)
        np.atleast_1d(axs)[0].set_ylabel("End 2 (Mb)")
        plt.tight_layout()
        plt.savefig(f"{OUT_DIR}/links_{block}.jpg", dpi=130)
        plt.close()

    spec = pd.concat(spec_rows) if spec_rows else pd.DataFrame()
    spec.to_csv(f"{OUT_DIR}/group_specific.tsv", sep="\t", index=False)
    print(f"group-specific links (Bonferroni p < {ALPHA}, >= {MIN_FRAC:.0%} of supporting group):")
    if len(spec):
        for g in ["A", "Het", "B"]:
            spec[g] = spec[g].astype(int).astype(str) + "/" + spec[f"n_{g}"].astype(int).astype(str)
        cols = ["block", "supporting_group", "lo_Mb", "hi_Mb", "span_Mb", "A", "Het", "B", "split_reads", "p_bonf"]
        print(spec.sort_values(["block", "supporting_group", "lo_Mb"])[cols].to_string(index=False))
    else:
        print("  none")


if __name__ == "__main__":
    main()
