#!/usr/bin/env python3
"""
Locate the genes between the GNS1 (g4524.t4) and SNF1 (g14399.t1) orthologues
in sunflower and lettuce, and find their Deinandra counterparts (syntenic-block
GENESPACE hit, else best BLAST hit, flagged separately). Also reports, along
lettuce chr 5 (WINDOW_MB around the interval), which Deinandra chromosome each
syntenic gene maps to and the Spearman rho of positions per chromosome.

Usage: gns1_snf1_interval.py
Output: OUT_DIR/gns1_snf1_intervening_<species>.tsv, OUT_DIR/gns1_snf1_lettuce_segment.tsv,
        printed summary
"""

import pandas as pd
from scipy.stats import spearmanr

BASE_DIR = "/private/groups/meyerlab/esteiger"

GS = f"{BASE_DIR}/inversions/synteny/genespace"
OUT_DIR = f"{BASE_DIR}/inversions/synteny"
GENES = {"GNS1": "g4524.t4", "SNF1": "g14399.t1"}
SPECIES = ["Helianthus_annuus", "Lactuca_sativa"]
WINDOW_MB = 3


def partner(syn, gene):
    h = syn[(syn.id2 == gene) & syn.blkID.notna()].sort_values("bitScore", ascending=False)
    return h.iloc[0]


def main():
    dbed = pd.read_csv(f"{GS}/bed/Deinandra.bed", sep="\t", header=None, names=["chr", "s", "e", "id"])
    dpos = dbed.set_index("id")
    a, b = dpos.loc[GENES["GNS1"]], dpos.loc[GENES["SNF1"]]
    print(f"D. increscens: GNS1 {a.chr}:{a.s:,}, SNF1 {b.chr}:{b.s:,} -> {abs(b.s - a.s)/1e3:.1f} kb apart")

    for sp in SPECIES:
        syn = pd.read_csv(f"{GS}/syntenicHits/{sp}_vs_Deinandra.synHits.txt.gz", sep="\t", low_memory=False)
        p1, p2 = partner(syn, GENES["GNS1"]), partner(syn, GENES["SNF1"])
        assert str(p1.chr1) == str(p2.chr1), "GNS1 and SNF1 partners on different chromosomes"
        lo, hi = sorted([p1.start1, p2.start1])
        bed = pd.read_csv(f"{GS}/bed/{sp}.bed", sep="\t", header=None, names=["chr", "s", "e", "id"])
        bed.chr = bed.chr.astype(str)
        mid = bed[(bed.chr == str(p1.chr1)) & (bed.s > lo) & (bed.s < hi)].copy()
        print(f"\n{sp}: GNS1 -> {p1.id1} ({p1.chr1}:{p1.start1:,}); SNF1 -> {p2.id1} ({p2.chr1}:{p2.start1:,}); "
              f"{(hi - lo)/1e6:.2f} Mb apart, {len(mid)} genes between")

        # Deinandra location of each intervening gene: syntenic hit first, else best BLAST hit
        s = syn[syn.blkID.notna() & syn.id1.isin(mid.id)].sort_values("bitScore", ascending=False).drop_duplicates("id1")
        blast = pd.read_csv(f"{GS}/syntenicHits/{sp}_vs_Deinandra.allBlast.txt.gz", sep="\t", low_memory=False)
        bb = blast[blast.id1.isin(mid.id)].sort_values("bitScore", ascending=False).drop_duplicates("id1")
        s["source"] = "syntenic"
        bb = bb[~bb.id1.isin(s.id1)].assign(source="best BLAST hit only")
        loc = pd.concat([s, bb])[["id1", "chr2", "start2", "id2", "bitScore", "source"]]
        mid = mid.merge(loc, left_on="id", right_on="id1", how="left").drop(columns="id1")
        mid["source"] = mid.source.fillna("no Deinandra hit")
        mid.to_csv(f"{OUT_DIR}/gns1_snf1_intervening_{sp}.tsv", sep="\t", index=False)
        print(mid.groupby(["source", "chr2"], dropna=False).size().rename("genes").to_string())

    # gene-by-gene order along lettuce chr 5 around the interval
    syn = pd.read_csv(f"{GS}/syntenicHits/Lactuca_sativa_vs_Deinandra.synHits.txt.gz", sep="\t", low_memory=False)
    p1, p2 = partner(syn, GENES["GNS1"]), partner(syn, GENES["SNF1"])
    lo, hi = min(p1.start1, p2.start1) - WINDOW_MB * 1e6, max(p1.start1, p2.start1) + WINDOW_MB * 1e6
    seg = syn[(syn.chr1 == p1.chr1) & syn.start1.between(lo, hi) & syn.blkID.notna()].sort_values("start1")
    seg = seg[["start1", "id1", "chr2", "start2", "isAnchor"]]
    seg.to_csv(f"{OUT_DIR}/gns1_snf1_lettuce_segment.tsv", sep="\t", index=False)
    print(f"\nLettuce {p1.chr1} {lo/1e6:.1f}-{hi/1e6:.1f} Mb: {len(seg)} syntenic hits, Deinandra chromosomes:")
    for chrom, g in seg.groupby("chr2"):
        rho = spearmanr(g.start1, g.start2).statistic if len(g) > 2 else float("nan")
        print(f"  {chrom}: {len(g)} hits, {g.start2.min()/1e6:.2f}-{g.start2.max()/1e6:.2f} Mb, "
              f"order vs lettuce Spearman rho = {rho:+.2f}")
    inside = seg[seg.start1.between(min(p1.start1, p2.start1), max(p1.start1, p2.start1))]
    print("  between the GNS1 and SNF1 orthologues: " +
          ", ".join(f"{k} {v}" for k, v in inside.chr2.value_counts().items()))


if __name__ == "__main__":
    main()
