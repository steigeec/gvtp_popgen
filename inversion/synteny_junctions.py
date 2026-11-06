#!/usr/bin/env python3
"""
Is each chromosome-scale scaffold's centromere a junction between two
ancestral chromosome segments (as on DEINVI_00006)?

Reference: lettuce (Lactuca sativa), which lies outside the genome
triplication shared by sunflower and Deinandra, so each Deinandra region has
one lettuce counterpart; sunflower (Helianthus annuus) is reported alongside
as support. Syntenic gene pairs = GENESPACE syntenic hits inside a syntenic
block (blkID assigned).

For each scaffold with a confident centromere call (call_centromeres.py):
starting from the syntenic gene nearest each edge of the call and moving
outwards, follow the run of genes whose partner lies on the SAME partner
chromosome, tolerating up to MAX_STRAY consecutive genes with another partner
(single misplaced genes, small transpositions). The centromere is classed as
  junction      the partner chromosome of the two flanking runs differs, and
                each flanking run has >= MIN_RUN genes
  continuous    the same partner chromosome flanks both sides, each run
                >= MIN_RUN genes
  possible      either case, but a flanking run has < MIN_RUN genes (e.g. a
                single gene of the left-arm partner just past the centromere,
                then a switch), so the junction may lie just off the call
Also reported: syntenic genes on each arm, the length (genes, Mb) of each
flanking run, and the gap without syntenic genes spanning the centromere.

Usage: synteny_junctions.py
Output: OUT (one row per scaffold x species), printed table
"""

import numpy as np
import pandas as pd

GS = "/private/groups/meyerlab/esteiger/inversions/synteny/genespace/syntenicHits"
CEN = "/private/groups/meyerlab/esteiger/inversions/centromere/centromeres.tsv"
OUT = "/private/groups/meyerlab/esteiger/inversions/synteny/synteny_junctions.tsv"
SPECIES = {"Lactuca_sativa": "lettuce", "Helianthus_annuus": "sunflower"}
LETTUCE = {f"CM0225{18 + i}.2": str(1 + i) for i in range(9)}   # CM022518.2 = lettuce chr 1, ...
MAX_STRAY, MIN_RUN = 2, 5


def hits(sp, cen):
    d = pd.read_csv(f"{GS}/{sp}_vs_Deinandra.synHits.txt.gz", sep="\t", low_memory=False,
                    usecols=["chr1", "chr2", "start2", "end2", "blkID"])
    d = d[d.blkID.notna() & d.chr2.isin(cen.index)].copy()
    d["x"] = (d.start2 + d.end2) / 2
    c = d.chr1.astype(str)
    d["partner"] = c.map(LETTUCE) if sp == "Lactuca_sativa" else c.str.extract(r"(\d+)\D*$")[0]
    return d


def flank_run(g):
    """g: hits ordered outwards from the centromere -> (partner, n genes, Mb spanned, start, end)"""
    p0, n, stray, first, last = g.partner.iloc[0], 0, 0, g.x.iloc[0], g.x.iloc[0]
    for p, x in zip(g.partner, g.x):
        if p == p0:
            n += 1; last = x; stray = 0
        else:
            stray += 1
            if stray > MAX_STRAY:
                break
    return p0, n, abs(last - first) / 1e6, first


def main():
    cen = pd.read_csv(CEN, sep="\t").query("confident").set_index("chr")
    rows = []
    for sp, name in SPECIES.items():
        d = hits(sp, cen)
        for chrom, c in cen.iterrows():
            g = d[d.chr2 == chrom]
            L = g[g.x < c.cen_start].sort_values("x", ascending=False)
            R = g[g.x > c.cen_end].sort_values("x")
            if not len(L) or not len(R):
                rows.append(dict(scaffold=chrom, species=name, cls="one arm without synteny",
                                 left_arm_genes=len(L), right_arm_genes=len(R)))
                continue
            pl, nl, ml, el = flank_run(L)
            pr, nr, mr, er = flank_run(R)
            cls = "possible" if min(nl, nr) < MIN_RUN else ("continuous" if pl == pr else "junction")
            rows.append(dict(scaffold=chrom, species=name, cen_Mb=f"{c.cen_start/1e6:.1f}-{c.cen_end/1e6:.1f}",
                             left_arm_genes=len(L), right_arm_genes=len(R),
                             left_partner=pl, left_run_genes=nl, left_run_Mb=round(ml, 1),
                             right_partner=pr, right_run_genes=nr, right_run_Mb=round(mr, 1),
                             gap_Mb=round((er - el) / 1e6, 1), cls=cls))
    r = pd.DataFrame(rows)
    r.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(r.to_string(index=False))
    w = r.pivot(index="scaffold", columns="species", values="cls")
    print("\nCentromere class by reference species:\n" + w.to_string())
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
