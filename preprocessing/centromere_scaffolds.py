#!/usr/bin/env python3
"""
Gene density and % LINE content along the 12 largest scaffolds (D. increscens
2n = 24, so n = 12 chromosomes), with a grey box marking each scaffold's
centromere call, to show where the centromeres sit genome-wide.

Tracks use the same inputs and conventions as the DEINVI_00006 figure
(gvtp_popgen/inversion/scaffold_windows.py): WIN_SIZE windows; genes /
repeats counted in a window only if they lie entirely inside it;
RepeatMasker hits flagged as overlapping (trailing "*") excluded.
Dashed grey line = genome-wide value (all scaffolds).

Grey box = centromere call from gvtp_popgen/inversion/call_centromeres.py
(inversions/centromere/centromeres.tsv; "confident" calls only): the densest
cluster of the centromere-associated LINE family (rnd-1_family-362), which
forms one tight cluster on nearly every chromosome-scale scaffold. As a
check, the scaffold's peak window of overall % LINE (smoothed over
SMOOTH windows) is reported for each scaffold.

Layout: N_COLS columns of scaffolds, in descending size order down each
column in turn; each column's x-axis spans its longest scaffold. Each scaffold's stack of tracks is labelled by a vertical grey strip
to its left carrying the scaffold name.

With --manhattan, each scaffold also gets the fog GEMMA Manhattan on top,
built exactly as annotated_manhattan.R exports it for the stacked figures:
p_lrt with genomic control (lambda applied if > 1.05), Bonferroni line at
0.05 / n SNPs, and peak-filtered points (non-outliers plus only the
Bonferroni SNPs retained in 40000_3 peaks). For drawing speed only, a random
LOW_P_KEEP of SNPs with -log10 p < LOW_P are plotted (they form the solid
band at the bottom either way). The --manhattan version also adds a
lettuce synteny track at the bottom of each stack: one point per GENESPACE
syntenic gene pair inside a syntenic block (Lactuca_sativa_vs_Deinandra
synHits, blkID assigned), at the gene's position on the scaffold, with y =
its partner's position along the lettuce chromosomes syntenic with that
scaffold (only chromosomes with >= MIN_SYN_GENES syntenic genes and >=
MIN_SYN_FRAC of the scaffold's syntenic genes; laid end to end in numerical
order, alternate chromosomes shaded; collinear runs appear as diagonals, their
slope giving orientation). Track height is proportional to the lettuce
sequence shown (SYN_UNITS_PER_GB, at least SYN_MIN_UNITS), and the figure
height follows from the tallest column (PX_PER_UNIT). (lettuce lies outside the genome
triplication shared by sunflower and Deinandra, so it stands in for the
ancestral chromosomes; see inversion/synteny_junctions.py).

Usage: centromere_scaffolds.py [--manhattan]
Output: OUT (or OUT_MANHATTAN)
"""

import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FAI = "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked.fai"
GFF = "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
RM_OUT = "/private/groups/meyerlab/esteiger/inversions/gene_density/DEINVI_genome.fasta.LTRfix.out"
CENTROMERES = "/private/groups/meyerlab/esteiger/inversions/centromere/centromeres.tsv"
OUT = "/private/groups/meyerlab/esteiger/inversions/centromere/centromere_scaffolds.jpg"
OUT_MANHATTAN = "/private/groups/meyerlab/esteiger/inversions/centromere/centromere_scaffolds_manhattan.jpg"
GEMMA_DIR = "/private/groups/meyerlab/esteiger/selection/gemma"
ASSOC = f"{GEMMA_DIR}/output/pheno6_lmm_results.assoc.txt"            # fog = pheno6
RETAINED = f"{GEMMA_DIR}/fog/outliers_for_annot_40000_3_fog.txt"
SNP_MAP = f"{GEMMA_DIR}/scaffold_snp_map.txt"
LOW_P, LOW_P_KEEP = 1.0, 0.1
N_CHROM, N_COLS, WIN_SIZE, SMOOTH = 12, 3, 500_000, 3
SPACER = 0.25        # gap between scaffold stacks (relative to track heights)
HEADROOM = 1.02      # y-axis headroom above the highest value
DPI = 180
WIDTH_PX_MANHATTAN = 2458           # exact output width of the --manhattan figure; height follows from content
PX_PER_UNIT = 142                   # pixels per height-ratio unit (--manhattan figure)
TRACK_UNITS = [0.75, 0.6, 0.5]      # GWAS, gene density, LINE track heights (--manhattan figure)
SYN_UNITS_PER_GB, SYN_MIN_UNITS = 2.0, 1.0      # lettuce synteny track height
MIN_SYN_GENES, MIN_SYN_FRAC = 10, 0.02          # lettuce chromosome shown on a scaffold's synteny track if both met
MARGINS = dict(left=0.08, right=0.995, top=0.995, bottom=0.022)   # figure fractions, --manhattan figure
LABEL_SIZE, TICK_SIZE, NAME_SIZE = 12, 9, 12
SMALL_LABEL = 7.5      # y-axis titles of the narrow --manhattan tracks
SHADE_ALPHA = 0.25
LETTUCE_BED = "/private/groups/meyerlab/esteiger/inversions/synteny/genespace/bed/Lactuca_sativa.bed"
SYNHITS = "/private/groups/meyerlab/esteiger/inversions/synteny/genespace/syntenicHits/Lactuca_sativa_vs_Deinandra.synHits.txt.gz"
LETTUCE = {f"CM0225{18 + i}.2": 1 + i for i in range(9)}   # CM022518.2 = lettuce chr 1, ...
STRIP_COLOR, STRIP_W, STRIP_GAP = "#d9d9d9", 0.016, 0.004   # scaffold label strip (figure fractions)


def lettuce_synteny(scaffolds):
    """per scaffold: dict(x = gene positions (Mb), y = partner positions along the scaffold's
    selected lettuce chromosomes laid end to end (Mb), chroms, offsets, lengths (Mb), total (Mb))"""
    bed = pd.read_csv(LETTUCE_BED, sep="\t", header=None, names=["chr", "s", "e", "id"])
    clen = (bed[bed.chr.isin(LETTUCE)].groupby("chr").e.max() / 1e6).rename(index=LETTUCE)
    d = pd.read_csv(SYNHITS, sep="\t", low_memory=False, usecols=["chr1", "start1", "end1", "chr2", "start2", "end2", "blkID"])
    d = d[d.blkID.notna() & d.chr2.isin(scaffolds) & d.chr1.astype(str).isin(LETTUCE)].copy()
    d["lchr"] = d.chr1.astype(str).map(LETTUCE)
    out = {}
    for s, g in d.groupby("chr2"):
        vc = g.lchr.value_counts()
        keep = sorted(vc[(vc >= MIN_SYN_GENES) & (vc / vc.sum() >= MIN_SYN_FRAC)].index)
        lens = clen[keep].values
        offs = np.concatenate([[0], np.cumsum(lens)[:-1]])
        g = g[g.lchr.isin(keep)]
        y = g.lchr.map(dict(zip(keep, offs))).values + (g.start1.values + g.end1.values) / 2e6
        out[s] = dict(x=(g.start2.values + g.end2.values) / 2e6, y=y, chroms=keep, offs=offs, lens=lens,
                      total=lens.sum(), dropped={int(k): int(v) for k, v in vc.items() if k not in keep})
    return out


def syn_units(info):
    return max(SYN_MIN_UNITS, SYN_UNITS_PER_GB * info["total"] / 1000)


def windows(df, L):
    """per WIN_SIZE window: count and summed length of features lying entirely inside it"""
    edges = np.arange(0, L + WIN_SIZE, WIN_SIZE); edges[-1] = L
    idx = np.searchsorted(edges, df.start.values, side="right") - 1
    inside = df.end.values <= edges[np.minimum(idx + 1, len(edges) - 1)]
    idx = idx[inside]
    n = np.bincount(idx, minlength=len(edges) - 1)[:len(edges) - 1]
    bp = np.bincount(idx, weights=(df.end.values - df.start.values + 1)[inside], minlength=len(edges) - 1)[:len(edges) - 1]
    return edges, n, bp


def fog_manhattan(scaffolds):
    """per scaffold: positions and -log10 p (peak-filtered, as in annotated_manhattan.R), and the Bonferroni line"""
    from scipy.stats import chi2
    g = pd.read_csv(ASSOC, sep="\t", usecols=["chr", "rs", "ps", "p_lrt"])
    x = chi2.isf(g.p_lrt.values, 1)
    lam = np.nanmedian(x) / chi2.ppf(0.5, 1)
    g["p"] = chi2.sf(x / lam, 1) if lam > 1.05 else g.p_lrt
    thresh = 0.05 / len(g)
    retained = set(pd.read_csv(RETAINED, sep=" ", header=None)[1].astype(str))
    keep = (g.p >= thresh) | g.rs.astype(str).isin(retained)
    smap = pd.read_csv(SNP_MAP, sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    num2scaf = dict(zip(smap[0], smap[2]))
    g["scaf"] = g.chr.map(num2scaf)
    g = g[keep & g.scaf.isin(scaffolds)]
    g["logp"] = -np.log10(g.p)
    rng = np.random.default_rng(1)
    g = g[(g.logp >= LOW_P) | (rng.random(len(g)) < LOW_P_KEEP)]
    print(f"lambda {lam:.3f}; Bonferroni -log10 p {-np.log10(thresh):.2f}; {keep.sum() - len(g)} low-p SNPs not drawn")
    return {s: (d.ps.values / 1e6, d.logp.values) for s, d in g.groupby("scaf")}, -np.log10(thresh)


def main(with_manhattan=False):
    fai = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"])
    scafs = fai.sort_values("len", ascending=False).head(N_CHROM).reset_index(drop=True)
    genome_bp = fai.len.sum()

    genes = []
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) >= 9 and p[2] == "gene":
                genes.append((p[0], int(p[3]), int(p[4])))
    genes = pd.DataFrame(genes, columns=["chrom", "start", "end"])

    rec = []
    with open(RM_OUT) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip()
            if i < 3 or not line or line.endswith("*"):
                continue
            p = line.split()
            if len(p) >= 11 and p[10].split("/")[0] == "LINE":
                rec.append((p[4], int(p[5]), int(p[6])))
    line_rm = pd.DataFrame(rec, columns=["chrom", "start", "end"])
    ref_genes = len(genes) / (genome_bp / 1e6)
    ref_line = (line_rm.end - line_rm.start + 1).sum() / genome_bp * 100

    cen = pd.read_csv(CENTROMERES, sep="\t").set_index("chr")
    if with_manhattan:
        man, bonf = fog_manhattan(set(scafs.chrom))
        man_top = max(v[1].max() for v in man.values()) * HEADROOM
        syn = lettuce_synteny(set(scafs.chrom))
        for s_, v in syn.items():
            print(f"{s_}: lettuce chromosomes shown {v['chroms']}; below cutoff (genes) {v['dropped']}")

    n_rows = int(np.ceil(N_CHROM / N_COLS))
    cols = [scafs.iloc[i * n_rows:(i + 1) * n_rows] for i in range(N_COLS)]
    tracks = (["man"] if with_manhattan else []) + ["genes", "line"] + (["syn"] if with_manhattan else [])
    # per-scaffold height ratios (the synteny track varies with the lettuce sequence shown);
    # shorter columns get a blank filler row at the bottom so all columns share one scale
    def scaf_ratios(chrom):
        return (TRACK_UNITS + [syn_units(syn[chrom])] if with_manhattan else [1, 0.8]) + [SPACER]
    col_units = [sum(sum(scaf_ratios(c)) for c in col.chrom) for col in cols]
    if with_manhattan:      # fixed width, height from content: exact canvas, explicit margins, no tight bbox
        h_px = round(max(col_units) * PX_PER_UNIT / (MARGINS["top"] - MARGINS["bottom"]))
        fig = plt.figure(figsize=(WIDTH_PX_MANHATTAN / DPI, h_px / DPI))
        outer = fig.add_gridspec(1, N_COLS, wspace=0.42, **MARGINS)
    else:
        fig = plt.figure(figsize=(7.5 * N_COLS, 2.4 * n_rows))
        outer = fig.add_gridspec(1, N_COLS, wspace=0.28)
    stacks = []      # (column, first axis, last axis, scaffold) for the label strips
    report = []
    for ci, col in enumerate(cols):
        xmax = col.len.max() / 1e6
        rr = [scaf_ratios(c) for c in col.chrom]
        filler = max(col_units) - col_units[ci]
        hr = [r for x in rr for r in x] + ([filler] if filler > 1e-9 else [])
        inner = outer[ci].subgridspec(len(hr), 1, height_ratios=hr, hspace=0.1)
        starts = np.cumsum([0] + [len(x) for x in rr])
        for ri, (chrom, L) in enumerate(zip(col.chrom, col.len)):
            axs = {t: fig.add_subplot(inner[starts[ri] + j]) for j, t in enumerate(tracks)}
            ax_g, ax_l = axs["genes"], axs["line"]
            if with_manhattan:
                ax_m = axs["man"]
                if chrom in man:
                    ax_m.scatter(*man[chrom], s=4, color="black", alpha=0.6, edgecolor="none", rasterized=True)
                ax_m.axhline(bonf, color="grey", ls="--", lw=1)
                ax_m.set_ylim(0, man_top)
                ax_m.set_ylabel("−log p", fontsize=SMALL_LABEL)
            if with_manhattan:
                ax_s, sy = axs["syn"], syn[chrom]
                ax_s.scatter(sy["x"], sy["y"], s=3, color="black", alpha=0.7, edgecolor="none", rasterized=True)
                for i, (o, n) in enumerate(zip(sy["offs"], sy["lens"])):
                    if i % 2:
                        ax_s.axhspan(o, o + n, color="#f2f2f2", lw=0, zorder=0)
                ax_s.set_ylim(0, sy["total"])
                ax_s.set_ylabel("Lettuce chr", fontsize=SMALL_LABEL + 1)
            e_g, n_g, _ = windows(genes[genes.chrom == chrom], L)
            e_l, _, bp_l = windows(line_rm[line_rm.chrom == chrom], L)
            dens = n_g / (np.diff(e_g) / 1e6)
            pct = bp_l / np.diff(e_l) * 100
            ax_g.stairs(dens, e_g / 1e6, fill=True, color="black")
            ax_l.stairs(pct, e_l / 1e6, fill=True, color="black")
            ax_g.axhline(ref_genes, color="grey", ls="--", lw=1, zorder=3)
            ax_l.axhline(ref_line, color="grey", ls="--", lw=1, zorder=3)
            box = None
            if chrom in cen.index and bool(cen.loc[chrom, "confident"]):
                box = (cen.loc[chrom, "cen_start"] / 1e6, cen.loc[chrom, "cen_end"] / 1e6)
                for ax in axs.values():
                    ax.axvspan(*box, color="#808080", alpha=SHADE_ALPHA, lw=0, zorder=0)
            sm = np.convolve(pct, np.ones(SMOOTH) / SMOOTH, mode="same")
            pk = (e_l[np.argmax(sm)] + WIN_SIZE / 2) / 1e6
            report.append(dict(scaffold=chrom, length_Mb=round(L / 1e6, 1),
                               centromere_Mb=f"{box[0]:.1f}-{box[1]:.1f}" if box else "none",
                               LINE_peak_Mb=round(pk, 2), peak_pct=round(sm.max(), 1),
                               peak_in_or_within_2Mb=bool(box and box[0] - 2 <= pk <= box[1] + 2)))
            for ax in axs.values():
                ax.set_xlim(0, xmax)
                ax.spines[["top", "right"]].set_visible(False)
                ax.tick_params(labelsize=TICK_SIZE)
                ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(2))
            if with_manhattan:
                sy = syn[chrom]
                axs["syn"].set_yticks(sy["offs"] + sy["lens"] / 2, [str(c) for c in sy["chroms"]])
                axs["syn"].tick_params(axis="y", length=0)
            ax_g.set_ylim(0, max(dens.max(), ref_genes) * HEADROOM)
            ax_l.set_ylim(0, max(pct.max(), ref_line) * HEADROOM)
            for t in tracks[:-1]:
                axs[t].tick_params(labelbottom=False)
            stacks.append((ci, axs[tracks[0]], axs[tracks[-1]], chrom))
            # short one-line titles: the --manhattan tracks are narrow
            ax_g.set_ylabel("Genes" if with_manhattan else "Genes\nper Mb",
                            fontsize=SMALL_LABEL if with_manhattan else LABEL_SIZE - 2)
            ax_l.set_ylabel("LINE" if with_manhattan else "% LINE",
                            fontsize=SMALL_LABEL if with_manhattan else LABEL_SIZE - 2)
            if ri == n_rows - 1:
                axs[tracks[-1]].set_xlabel("Position (Mb)", fontsize=LABEL_SIZE)
        fig.align_ylabels([a for a in fig.axes if a.get_subplotspec().get_gridspec() is inner])

    # vertical scaffold-name strips, just left of each column's y-axis labels
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    for ci in range(N_COLS):
        col_axes = [a for c, top, bot, _ in stacks if c == ci for a in (top, bot)]
        left = min(inv.transform(a.get_tightbbox(fig.canvas.get_renderer()))[0, 0] for a in col_axes)
        for c, top, bot, chrom in stacks:
            if c != ci:
                continue
            y0, y1 = bot.get_position().y0, top.get_position().y1
            sx = left - STRIP_GAP - STRIP_W
            fig.patches.append(matplotlib.patches.Rectangle((sx, y0), STRIP_W, y1 - y0, transform=fig.transFigure,
                                                            facecolor=STRIP_COLOR, edgecolor="none", zorder=0))
            fig.text(sx + STRIP_W / 2, (y0 + y1) / 2, chrom, rotation=90, ha="center", va="center",
                     fontsize=NAME_SIZE, fontweight="bold")
    out = OUT_MANHATTAN if with_manhattan else OUT
    fig.savefig(out, dpi=DPI, bbox_inches=None if with_manhattan else "tight")
    pd.set_option("display.width", 200)
    print(pd.DataFrame(report).to_string(index=False))
    print(f"genome-wide: {ref_genes:.1f} genes/Mb, {ref_line:.2f}% LINE\nSaved: {out}")


if __name__ == "__main__":
    main(with_manhattan="--manhattan" in sys.argv[1:])
