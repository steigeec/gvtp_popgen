#!/usr/bin/env python3
"""
Overlap of the per-variable candidate gene sets in gene_summary_all.tsv
(BIO15, BIO18, fog):
  1. pairwise hypergeometric tests (upper tail, P(X >= observed overlap))
     against the genome-wide gene universe (GFF "gene" features);
  2. fraction of shared genes inside vs outside the DEINVI_00006:44-49 Mb
     putative inversion.

Usage: python3 gene_overlap_across_variables.py
Writes:
  - gene_overlap_summary.tsv (set sizes, pairwise/triple overlaps, hypergeometric p-values)
  - gene_overlap_shared_genes.tsv (every gene shared by >=2 variables, which variables, inversion flag)
  - gene_overlap_venn.png (3-set proportional-ish Venn diagram)
"""
import csv
import re
from collections import defaultdict

try:
    from scipy.stats import hypergeom
except ImportError:
    hypergeom = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

BASE_DIR = "/private/groups/meyerlab/esteiger"

GEMMA_DIR = f"{BASE_DIR}/selection/gemma"
GENE_SUMMARY = f"{GEMMA_DIR}/gene_summary_all.tsv"
GFF = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
VARS = ["BIO15", "BIO18", "fog"]
INVERSION_SCAFFOLD = "DEINVI_00006"
INVERSION_START, INVERSION_END = 44_000_000, 49_000_000

COORD_RE = re.compile(r'^([^:]+):(\d+)-(\d+)$')


def in_inversion(coord):
    m = COORD_RE.match(coord)
    if not m:
        return False
    scaf, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    return scaf == INVERSION_SCAFFOLD and not (end < INVERSION_START or start > INVERSION_END)


def count_genome_genes(gff_path):
    n = 0
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2 and fields[2] == "gene":
                n += 1
    return n


def load_gene_sets(path):
    sets_ = defaultdict(set)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sets_[row["Variable"]].add(row["Gene_Coordinates"])
    return sets_


def main():
    gene_universe_n = count_genome_genes(GFF)
    print(f"Genome-wide gene universe (GFF 'gene' features): {gene_universe_n}")

    gene_sets = load_gene_sets(GENE_SUMMARY)
    for v in VARS:
        print(f"{v}: {len(gene_sets[v])} candidate genes")

    # --- pairwise overlap + hypergeometric enrichment ---
    pair_rows = []
    for i in range(len(VARS)):
        for j in range(i + 1, len(VARS)):
            a, b = VARS[i], VARS[j]
            overlap = gene_sets[a] & gene_sets[b]
            k, n_, x = len(gene_sets[a]), len(gene_sets[b]), len(overlap)
            expected = k * n_ / gene_universe_n
            pval = hypergeom.sf(x - 1, gene_universe_n, k, n_) if hypergeom else float("nan")
            pair_rows.append((a, b, k, n_, x, expected, pval))
            print(f"{a} vs {b}: observed overlap={x}, expected by chance={expected:.2f}, "
                  f"hypergeometric P(X>=observed)={pval:.3e}")

    triple = gene_sets["BIO15"] & gene_sets["BIO18"] & gene_sets["fog"]
    any_pair_union = ((gene_sets["BIO15"] & gene_sets["BIO18"]) |
                       (gene_sets["BIO15"] & gene_sets["fog"]) |
                       (gene_sets["BIO18"] & gene_sets["fog"]))
    print(f"\nShared by all 3 variables: {len(triple)}")
    print(f"Shared by >=2 variables (union of pairwise intersections): {len(any_pair_union)}")

    inv_triple = sum(1 for g in triple if in_inversion(g))
    inv_shared = sum(1 for g in any_pair_union if in_inversion(g))
    print(f"  ...of which inside {INVERSION_SCAFFOLD}:{INVERSION_START}-{INVERSION_END} inversion: "
          f"{inv_triple}/{len(triple)} (triple), {inv_shared}/{len(any_pair_union)} (>=2 variables)")

    # --- write summary table ---
    with open(f"{GEMMA_DIR}/gene_overlap_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Variable_A", "Variable_B", "N_A", "N_B", "Observed_Overlap",
                    "Expected_By_Chance", "Hypergeometric_P"])
        for row in pair_rows:
            w.writerow(row)
        w.writerow([])
        w.writerow(["Shared_by_all_3_variables", len(triple)])
        w.writerow(["...inside_inversion_region", inv_triple])
        w.writerow(["Shared_by_>=2_variables", len(any_pair_union)])
        w.writerow(["...inside_inversion_region", inv_shared])
        w.writerow(["Gene_universe_size (GFF gene features)", gene_universe_n])
    print(f"\nWrote {GEMMA_DIR}/gene_overlap_summary.tsv")

    # --- write shared-gene detail table ---
    all_shared = sorted(any_pair_union)
    with open(f"{GEMMA_DIR}/gene_overlap_shared_genes.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene_Coordinates", "BIO15", "BIO18", "fog", "N_Variables", "In_Inversion_Region"])
        for g in all_shared:
            flags = [v for v in VARS if g in gene_sets[v]]
            w.writerow([g, "Y" if "BIO15" in flags else "",
                        "Y" if "BIO18" in flags else "",
                        "Y" if "fog" in flags else "",
                        len(flags), "Y" if in_inversion(g) else ""])
    print(f"Wrote {GEMMA_DIR}/gene_overlap_shared_genes.tsv")

    # --- 3-set Venn diagram (drawn manually; matplotlib-venn not installed) ---
    fig, ax = plt.subplots(figsize=(7, 6))
    centers = {"BIO15": (-0.7, 0.4), "BIO18": (0.7, 0.4), "fog": (0, -0.5)}
    colors = {"BIO15": "#4D4D4D", "BIO18": "#969696", "fog": "#1A1A1A"}
    radius = 1.5
    centroid = (sum(c[0] for c in centers.values()) / 3, sum(c[1] for c in centers.values()) / 3)
    # names outside the circles; single-variable counts at 2*center - centroid
    name_pos = {"BIO15": (-1.5, 2.2), "BIO18": (1.5, 2.2), "fog": (0, -2.3)}
    for v in VARS:
        ax.add_patch(Circle(centers[v], radius, alpha=0.45, facecolor=colors[v],
                             edgecolor="black", linewidth=1.2, label=v))
        cx, cy = centers[v]
        ox, oy = 2 * cx - centroid[0], 2 * cy - centroid[1]
        ax.text(*name_pos[v], f"{v}\n(n={len(gene_sets[v])})", ha="center", va="center",
                fontsize=12, fontweight="bold", color="black")
        only_n = len(gene_sets[v] - set().union(*[gene_sets[o] for o in VARS if o != v]))
        ax.text(ox, oy, str(only_n), ha="center", va="center", fontsize=11)

    import math

    def pairwise_only_label_pos(a, b, c, r):
        """Label position inside the A/B-only region (in A and B, outside C):
        part way from the A/B center-midpoint toward the A/B circle
        intersection point farther from C."""
        ax_, ay_ = a
        bx_, by_ = b
        cx_, cy_ = c
        dx, dy = bx_ - ax_, by_ - ay_
        d = math.hypot(dx, dy)
        mx, my = (ax_ + bx_) / 2, (ay_ + by_) / 2
        h = math.sqrt(max(r**2 - (d / 2) ** 2, 0))
        ux, uy = -dy / d, dx / d  # unit vector perpendicular to A->B
        p1 = (mx + h * ux, my + h * uy)
        p2 = (mx - h * ux, my - h * uy)
        far = max([p1, p2], key=lambda p: math.hypot(p[0] - cx_, p[1] - cy_))
        # t = 0.5 sits too close to C's boundary
        t = 0.65
        return (mx + t * (far[0] - mx), my + t * (far[1] - my))

    b15_b18 = len(gene_sets["BIO15"] & gene_sets["BIO18"] - gene_sets["fog"])
    b15_fog = len(gene_sets["BIO15"] & gene_sets["fog"] - gene_sets["BIO18"])
    b18_fog = len(gene_sets["BIO18"] & gene_sets["fog"] - gene_sets["BIO15"])
    ax.text(*pairwise_only_label_pos(centers["BIO15"], centers["BIO18"], centers["fog"], radius),
            str(b15_b18), ha="center", va="center", fontsize=11)
    ax.text(*pairwise_only_label_pos(centers["BIO15"], centers["fog"], centers["BIO18"], radius),
            str(b15_fog), ha="center", va="center", fontsize=11)
    ax.text(*pairwise_only_label_pos(centers["BIO18"], centers["fog"], centers["BIO15"], radius),
            str(b18_fog), ha="center", va="center", fontsize=11)
    ax.text(*centroid, str(len(triple)), ha="center", va="center", fontsize=13, fontweight="bold")

    ax.set_xlim(-3, 3)
    ax.set_ylim(-2.9, 2.9)
    ax.set_aspect("equal")
    ax.axis("off")
    print(f"Venn context for the caption: max pairwise environmental Pearson r = 0.54; "
          f"gene universe n={gene_universe_n}")
    # crop the saved image to its non-white content plus a margin
    from PIL import Image, ImageChops
    png = f"{GEMMA_DIR}/gene_overlap_venn.png"
    fig.savefig(png, dpi=200)
    img = Image.open(png).convert("RGB")
    box = ImageChops.difference(img, Image.new("RGB", img.size, "white")).getbbox()
    pad = 12
    img = img.crop((max(box[0] - pad, 0), max(box[1] - pad, 0),
                    min(box[2] + pad, img.width), min(box[3] + pad, img.height)))
    img.save(png)
    img.save(f"{GEMMA_DIR}/gene_overlap_venn.jpg", quality=95)
    print(f"Wrote {png} and .jpg ({img.width} x {img.height} px)")


if __name__ == "__main__":
    main()
