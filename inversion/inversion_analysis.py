#!/usr/bin/env python3
"""
Combined gene-density / repeat-content / boundary analysis for the candidate
inversion region on DEINVI_00006, vs. the rest of the genome.

This merges what used to be three separate scripts (gene_density_analysis.py,
repeat_analysis.py, boundary_analysis.py) into one, so that the genome FASTA
and the (very large) RepeatMasker .out file are each parsed exactly once and
shared across sections instead of being re-read per script.

DEINVI_00006 carries TWO adjacent putative inversions:
  - INV2 (not the region of interest): ~13.0-43.0 Mb
  - INV1 (the region of interest):      ~43.0-49.3 Mb
INV2 is excluded from the "background" comparison group below (rather than
silently folded into "rest of DEINVI_00006"), since it is presumably also
non-recombining and would otherwise contaminate the control group.

Sections:
  1. Shared parameters & parsing (genome FASTA, GFF genes, RepeatMasker .out)
  2. Gene density analysis
  3. Repeat / TE content analysis
  4. Fine-scale boundary analysis
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# ════════════════════════════════════════════════════════════════════════════
# 1. Shared parameters & parsing
# ════════════════════════════════════════════════════════════════════════════

GFF        = "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
FASTA      = "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked"
RM_OUT     = "/private/groups/meyerlab/esteiger/inversions/gene_density/DEINVI_genome.fasta.LTRfix.out"
CHROM      = "DEINVI_00006"
LD_START   = 43_000_000   # bp -- INV1 (region of interest)
LD_END     = 49_300_000   # bp
INV2_START = 13_000_000   # bp -- INV2 (second, adjacent putative inversion; not of interest)
INV2_END   = 43_000_000   # bp
WIN_SIZE   = 500_000      # 500 kb windows for the genome-wide sliding-window plots
GAP_MB     = 0.2          # tiny visual gap drawn between the INV2 and INV1 grey boxes
OUT_DIR    = "/private/groups/meyerlab/esteiger/inversions/gene_density"

REGION_COLOR = {"INV1": "#e07b54", "INV2": "#9b59b6", "background": "#4a8db5"}


def classify_region(row):
    if row.start >= LD_START and row.end <= LD_END:
        return "INV1"
    if row.start >= INV2_START and row.end <= INV2_END:
        return "INV2"
    return "background"


# ── 1a. Parse genome sizes from FASTA (shared by gene density & repeat sections) ──
print("Parsing genome sizes...", flush=True)
chrom_sizes = {}
current = None
length = 0
with open(FASTA) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if current is not None:
                chrom_sizes[current] = length
            current = line[1:].split()[0]
            length = 0
        else:
            length += len(line)
    if current is not None:
        chrom_sizes[current] = length
print(f"  {len(chrom_sizes)} chromosomes/scaffolds found")
print(f"  {CHROM}: {chrom_sizes[CHROM]:,} bp")
chrom6_len = chrom_sizes[CHROM]

# ── 1b. Parse gene coordinates from GFF (gene density section only) ──────────
print("Parsing gene features from GFF...", flush=True)
genes = []
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 9:
            continue
        if parts[2] != "gene":
            continue
        genes.append((parts[0], int(parts[3]), int(parts[4])))

genes = pd.DataFrame(genes, columns=["chrom", "start", "end"])
genes["length"] = genes["end"] - genes["start"] + 1
print(f"  {len(genes):,} genes parsed")


# ── 1c. Parse RepeatMasker .out ONCE (shared by repeat & boundary sections) ──
def broad_class(cls):
    c = cls.split("/")[0] if "/" in cls else cls
    if c == "LTR":
        sub = cls.split("/")[1] if "/" in cls else "Unknown"
        if "Gypsy" in sub:
            return "LTR/Gypsy"
        if "Copia" in sub:
            return "LTR/Copia"
        return "LTR/Other"
    if c in ("Simple_repeat", "Low_complexity", "Satellite"):
        return c
    if c == "LINE":
        return "LINE"
    if c == "SINE":
        return "SINE"
    if c == "DNA":
        return "DNA"
    if c == "RC":
        return "DNA"       # rolling-circle (Helitron) -> lump with DNA
    if c == "PLE":
        return "DNA"
    if c == "rRNA":
        return "Other"
    if c == "Unknown":
        return "Unknown"
    return "Other"


CLASS_ORDER = ["LTR/Gypsy", "LTR/Copia", "LTR/Other", "LINE", "SINE",
               "DNA", "Satellite", "Simple_repeat", "Low_complexity", "Unknown", "Other"]

CLASS_COLORS = {
    "LTR/Gypsy":      "#c0392b",
    "LTR/Copia":      "#e67e22",
    "LTR/Other":      "#f1c40f",
    "LINE":           "#27ae60",
    "SINE":           "#2ecc71",
    "DNA":            "#2980b9",
    "Satellite":      "#8e44ad",
    "Simple_repeat":  "#bdc3c7",
    "Low_complexity": "#95a5a6",
    "Unknown":        "#7f8c8d",
    "Other":          "#ecf0f1",
}

print("Parsing RepeatMasker .out...", flush=True)
records = []
with open(RM_OUT) as fh:
    for i, line in enumerate(fh):
        if i < 3:
            continue
        line = line.rstrip()
        if not line:
            continue
        # Strip trailing asterisk (overlapping lower-score hit) -- exclude from coverage
        if line.endswith("*"):
            continue
        parts = line.split()
        if len(parts) < 11:
            continue
        records.append((parts[4], int(parts[5]), int(parts[6]), float(parts[1]), parts[9], parts[10]))

rm_all = pd.DataFrame(records, columns=["chrom", "start", "end", "div", "name", "class"])
rm_all["length"] = rm_all["end"] - rm_all["start"] + 1
rm_all["broad"] = rm_all["class"].apply(broad_class)
print(f"  {len(rm_all):,} non-overlapping repeat elements parsed")

chrom6_rm = rm_all[rm_all["chrom"] == CHROM].copy()   # used by both repeat & boundary sections
other_rm = rm_all[rm_all["chrom"] != CHROM]


# ════════════════════════════════════════════════════════════════════════════
# 2. Gene density analysis
# ════════════════════════════════════════════════════════════════════════════

def density(n, bp):
    return n / (bp / 1e6)


chrom6_genes = genes[genes["chrom"] == CHROM].copy()

in_ld = (chrom6_genes["start"] >= LD_START) & (chrom6_genes["end"] <= LD_END)
in_inv2 = (chrom6_genes["start"] >= INV2_START) & (chrom6_genes["end"] <= INV2_END)
in_bg = ~in_ld & ~in_inv2

ld_genes = chrom6_genes[in_ld]
inv2_genes = chrom6_genes[in_inv2]
bg_genes = chrom6_genes[in_bg]          # true background: chr6 excluding BOTH inversions
non_ld_genes = chrom6_genes[~in_ld]     # kept for backwards-compat summary row only
other_genes = genes[genes["chrom"] != CHROM]

ld_size = LD_END - LD_START
inv2_size = INV2_END - INV2_START
bg_size = chrom6_len - ld_size - inv2_size
non_ld_size = chrom6_len - ld_size
other_size = sum(v for k, v in chrom_sizes.items() if k != CHROM)

ld_dens = density(len(ld_genes), ld_size)
inv2_dens = density(len(inv2_genes), inv2_size)
bg_dens = density(len(bg_genes), bg_size)
non_ld_dens = density(len(non_ld_genes), non_ld_size)
other_dens = density(len(other_genes), other_size)
chrom6_dens = density(len(chrom6_genes), chrom6_len)

print(f"\nGene density summary (genes / Mb):")
print(f"  INV1 / region of interest ({LD_START/1e6:.1f}-{LD_END/1e6:.1f} Mb on {CHROM}):  {ld_dens:.2f}  ({len(ld_genes)} genes in {ld_size/1e6:.1f} Mb)")
print(f"  INV2 / second putative inversion ({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f} Mb):  {inv2_dens:.2f}  ({len(inv2_genes)} genes in {inv2_size/1e6:.1f} Mb)")
print(f"  True background on {CHROM} (excl. both inversions): {bg_dens:.2f}  ({len(bg_genes)} genes in {bg_size/1e6:.1f} Mb)")
print(f"  [Old definition] Rest of {CHROM} (incl. INV2):        {non_ld_dens:.2f}  ({len(non_ld_genes)} genes in {non_ld_size/1e6:.1f} Mb)")
print(f"  Whole {CHROM}:                    {chrom6_dens:.2f}")
print(f"  Rest of genome (other chroms):   {other_dens:.2f}  ({len(other_genes)} genes in {other_size/1e6:.0f} Mb)")

# ── Sliding-window gene density along DEINVI_00006 ────────────────────────────
print("\nComputing sliding-window gene density on DEINVI_00006...", flush=True)
windows = []
for start in range(0, chrom6_len, WIN_SIZE):
    end = min(start + WIN_SIZE, chrom6_len)
    n = ((chrom6_genes["start"] >= start) & (chrom6_genes["end"] <= end)).sum()
    windows.append({"start": start, "end": end, "n_genes": n,
                     "density": density(n, end - start)})
gd_win_df = pd.DataFrame(windows)
gd_win_df["mid"] = (gd_win_df["start"] + gd_win_df["end"]) / 2
gd_win_df["region"] = gd_win_df.apply(classify_region, axis=1)

# ── Plot 1: sliding-window gene density along DEINVI_00006 ───────────────────
print("Generating plot: gene density along DEINVI_00006...", flush=True)

fig, ax = plt.subplots(figsize=(14, 4))

ax.axvspan(INV2_START / 1e6, LD_START / 1e6 - GAP_MB, color="#808080", alpha=0.15, zorder=0)
ax.axvspan(LD_START / 1e6 + GAP_MB, LD_END / 1e6, color="#808080", alpha=0.15, zorder=0)

ax.bar(gd_win_df["mid"] / 1e6, gd_win_df["density"], width=WIN_SIZE / 1e6 * 0.9,
       color="black", edgecolor="none", alpha=0.85)

ax.axhline(other_dens, color="black", linestyle="--", linewidth=1.2,
           label="Genome-wide density")

ax.set_xlabel(f"Position on {CHROM} (Mb)", fontsize=16)
ax.set_ylabel("Gene density (genes / Mb)", fontsize=16)
ax.set_xlim(0, chrom6_len / 1e6)
ax.tick_params(axis="both", labelsize=13)
ax.legend(fontsize=13, loc="upper right")

plt.tight_layout()
out1 = f"{OUT_DIR}/gene_density_DEINVI00006_sliding_window.jpg"
plt.savefig(out1, dpi=150)
plt.close()
print(f"  Saved: {out1}")

# ── Plot 2: bar chart comparing regional gene densities ──────────────────────
print("Generating plot: regional gene density comparison...", flush=True)

regions = [
    (f"INV1\n({LD_START/1e6:.1f}-{LD_END/1e6:.1f} Mb)", ld_dens, len(ld_genes), ld_size, "#e07b54"),
    (f"INV2\n({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f} Mb)", inv2_dens, len(inv2_genes), inv2_size, "#9b59b6"),
    (f"Background\n(rest of {CHROM})", bg_dens, len(bg_genes), bg_size, "#4a8db5"),
    ("Rest of\ngenome", other_dens, len(other_genes), other_size, "#6aab6a"),
]

fig, ax = plt.subplots(figsize=(7, 5))
xpos = np.arange(len(regions))
bars = ax.bar(xpos, [r[1] for r in regions], color=[r[4] for r in regions],
              width=0.55, edgecolor="white", linewidth=0.5)

for bar, (label, dens, n, bp, col) in zip(bars, regions):
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.05,
            f"{dens:.1f}\n({n} genes)", ha="center", va="bottom", fontsize=9)

ax.set_xticks(xpos)
ax.set_xticklabels([r[0] for r in regions], fontsize=10)
ax.set_ylabel("Gene density (genes / Mb)", fontsize=12)
ax.set_title(f"Gene density: INV1 vs INV2 vs true background vs genome", fontsize=12)
ax.set_ylim(0, max(r[1] for r in regions) * 1.25)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
out2 = f"{OUT_DIR}/gene_density_regional_comparison.jpg"
plt.savefig(out2, dpi=150)
plt.close()
print(f"  Saved: {out2}")

# ── Statistical tests: window-level, Mann-Whitney U ──────────────────────────
print("\nStatistical comparison (window-level, Mann-Whitney U):")
ld_win = gd_win_df[gd_win_df["region"] == "INV1"]["density"]
inv2_win = gd_win_df[gd_win_df["region"] == "INV2"]["density"]
bg_win = gd_win_df[gd_win_df["region"] == "background"]["density"]
non_ld_win = gd_win_df[gd_win_df["region"] != "INV1"]["density"]   # old (contaminated) comparison, kept for reference

u_stat, p_val = stats.mannwhitneyu(ld_win, bg_win, alternative="two-sided")
print(f"  INV1 windows:        mean = {ld_win.mean():.2f} genes/Mb, n = {len(ld_win)}")
print(f"  True background windows (excl. INV2): mean = {bg_win.mean():.2f} genes/Mb, n = {len(bg_win)}")
print(f"  Mann-Whitney U = {u_stat:.0f}, p = {p_val:.4g}   <-- corrected primary comparison")

u_stat_old, p_val_old = stats.mannwhitneyu(ld_win, non_ld_win, alternative="two-sided")
print(f"  [Old/contaminated] Non-INV1 windows (incl. INV2): mean = {non_ld_win.mean():.2f} genes/Mb, n = {len(non_ld_win)}")
print(f"  [Old/contaminated] Mann-Whitney U = {u_stat_old:.0f}, p = {p_val_old:.4g}")

u_inv2, p_inv2 = stats.mannwhitneyu(inv2_win, bg_win, alternative="two-sided")
print(f"  INV2 windows:         mean = {inv2_win.mean():.2f} genes/Mb, n = {len(inv2_win)}")
print(f"  INV2 vs background: Mann-Whitney U = {u_inv2:.0f}, p = {p_inv2:.4g}")

gene_summary = {
    "Region": [
        f"INV1 / region of interest ({LD_START/1e6:.1f}-{LD_END/1e6:.1f}Mb, DEINVI_00006)",
        f"INV2 / second putative inversion ({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f}Mb, DEINVI_00006)",
        f"True background (DEINVI_00006, excl. both inversions)",
        "Rest of genome",
    ],
    "Size_Mb": [ld_size / 1e6, inv2_size / 1e6, bg_size / 1e6, other_size / 1e6],
    "N_genes": [len(ld_genes), len(inv2_genes), len(bg_genes), len(other_genes)],
    "Genes_per_Mb": [ld_dens, inv2_dens, bg_dens, other_dens],
}
pd.DataFrame(gene_summary).to_csv(f"{OUT_DIR}/gene_density_summary.tsv", sep="\t", index=False)
print(f"\nSummary saved to {OUT_DIR}/gene_density_summary.tsv")

# ── Plot 3: box-whisker of window-level density, INV1 vs INV2 vs background ──
print("\nGenerating plot: box-whisker of window-level gene density...", flush=True)

rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(5.5, 5))

box_data = [ld_win.values, inv2_win.values, bg_win.values]
box_color = "#808080"
bp = ax.boxplot(box_data, positions=[0, 1, 2], widths=0.5, showfliers=False,
                 patch_artist=True, medianprops=dict(color="black", linewidth=1.5))
for patch in bp["boxes"]:
    patch.set_facecolor(box_color)
    patch.set_alpha(0.35)

for i, vals in enumerate(box_data):
    jitter = rng.uniform(-0.12, 0.12, size=len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals, s=18,
               color=box_color, alpha=0.7, edgecolor="white", linewidth=0.4)

ax.axhline(other_dens, color="black", linestyle="--", linewidth=1.2)
xlim = ax.get_xlim()
ax.text(xlim[0] + 0.05, other_dens, "Genome-wide density", ha="left", va="bottom", fontsize=9)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels([f"43-49.3Mb region\nn={len(ld_win)} windows",
                     f"13-43Mb region\nn={len(inv2_win)} windows",
                     f"Background\n(rest of {CHROM})\nn={len(bg_win)} windows"], fontsize=9)
ax.set_ylabel("Gene density per Mb, calculated per 500kb window", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
out3 = f"{OUT_DIR}/gene_density_boxplot.jpg"
plt.savefig(out3, dpi=150)
plt.close()
print(f"  Saved: {out3}")


# ════════════════════════════════════════════════════════════════════════════
# 3. Repeat / TE content analysis
# ════════════════════════════════════════════════════════════════════════════

in_ld = (chrom6_rm["start"] >= LD_START) & (chrom6_rm["end"] <= LD_END)
in_inv2 = (chrom6_rm["start"] >= INV2_START) & (chrom6_rm["end"] <= INV2_END)
in_bg = ~in_ld & ~in_inv2

ld_rm = chrom6_rm[in_ld]
inv2_rm = chrom6_rm[in_inv2]
bg_rm = chrom6_rm[in_bg]        # true background: chr6 excluding BOTH inversions
non_ld_rm = chrom6_rm[~in_ld]   # kept for backwards-compat reference only


def repeat_fraction(df, size_bp):
    """% of region covered by each broad class (summed bp / region size)."""
    covered = df.groupby("broad")["length"].sum()
    return (covered / size_bp * 100).reindex(CLASS_ORDER, fill_value=0)


ld_frac = repeat_fraction(ld_rm, ld_size)
inv2_frac = repeat_fraction(inv2_rm, inv2_size)
bg_frac = repeat_fraction(bg_rm, bg_size)
non_ld_frac = repeat_fraction(non_ld_rm, non_ld_size)
other_frac = repeat_fraction(other_rm, other_size)

total_ld = ld_frac.sum()
total_inv2 = inv2_frac.sum()
total_bg = bg_frac.sum()
total_non_ld = non_ld_frac.sum()
total_other = other_frac.sum()

print(f"\nTotal repeat content (% bases covered):")
print(f"  INV1 / region of interest ({LD_START/1e6:.1f}-{LD_END/1e6:.1f} Mb):  {total_ld:.1f}%")
print(f"  INV2 / second putative inversion ({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f} Mb):  {total_inv2:.1f}%")
print(f"  True background on {CHROM} (excl. both inversions): {total_bg:.1f}%")
print(f"  [Old definition] Rest of {CHROM} (incl. INV2):        {total_non_ld:.1f}%")
print(f"  Rest of genome:         {total_other:.1f}%")
print(f"\nBy class (%):")
repeat_class_summary = pd.DataFrame({"INV1": ld_frac, "INV2": inv2_frac, "Background_chr6": bg_frac, "Rest_genome": other_frac})
print(repeat_class_summary.to_string(float_format="%.2f"))

# ── Sliding-window repeat content along DEINVI_00006 ──────────────────────────
print("\nComputing sliding windows on DEINVI_00006...", flush=True)
windows = []
for wstart in range(0, chrom6_len, WIN_SIZE):
    wend = min(wstart + WIN_SIZE, chrom6_len)
    wsize = wend - wstart
    mask = (chrom6_rm["start"] >= wstart) & (chrom6_rm["end"] <= wend)
    sub = chrom6_rm[mask]
    row = {"start": wstart, "end": wend, "mid": (wstart + wend) / 2}
    for cls in CLASS_ORDER:
        bp = sub.loc[sub["broad"] == cls, "length"].sum()
        row[cls] = bp / wsize * 100
    row["total"] = sub["length"].sum() / wsize * 100
    windows.append(row)
rep_win_df = pd.DataFrame(windows)
rep_win_df["region"] = rep_win_df.apply(classify_region, axis=1)
rep_x = rep_win_df["mid"].values / 1e6

# ── Plot: stacked area -- repeat class composition along DEINVI_00006 ────────
print("Generating plot: stacked repeat composition along DEINVI_00006...")
fig, ax = plt.subplots(figsize=(14, 5))
ax.axvspan(LD_START / 1e6, LD_END / 1e6, alpha=0.15, color="#e07b54", zorder=0)
ax.axvspan(INV2_START / 1e6, INV2_END / 1e6, alpha=0.15, color="#9b59b6", zorder=0)

bottom = np.zeros(len(rep_win_df))
for cls in CLASS_ORDER:
    y = rep_win_df[cls].values
    ax.bar(rep_x, y, bottom=bottom, width=WIN_SIZE / 1e6 * 0.92,
           color=CLASS_COLORS[cls], label=cls, edgecolor="none", alpha=0.9)
    bottom += y

ax.set_xlabel(f"Position on {CHROM} (Mb)", fontsize=12)
ax.set_ylabel("Repeat content (% bases)", fontsize=12)
ax.set_title(f"Repeat class composition along {CHROM} ({WIN_SIZE//1000} kb windows)", fontsize=13)
ax.set_xlim(0, chrom6_len / 1e6)

inv1_label = f"INV1 ({LD_START/1e6:.1f}-{LD_END/1e6:.1f} Mb)"
inv2_label = f"INV2 ({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f} Mb)"
inv1_patch = mpatches.Patch(color="#e07b54", alpha=0.4, label=inv1_label)
inv2_patch = mpatches.Patch(color="#9b59b6", alpha=0.4, label=inv2_label)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=[inv1_patch, inv2_patch] + handles, labels=[inv1_label, inv2_label] + labels,
          fontsize=8, loc="upper left", ncol=2, framealpha=0.8)

plt.tight_layout()
out1 = f"{OUT_DIR}/repeat_composition_DEINVI00006_sliding_window.jpg"
plt.savefig(out1, dpi=150)
plt.close()
print(f"  Saved: {out1}")

# ── Plot: total repeat content sliding window ─────────────────────────────────
print("Generating plot: total repeat content along DEINVI_00006...")
fig, ax = plt.subplots(figsize=(14, 4))
ax.axvspan(INV2_START / 1e6, LD_START / 1e6 - GAP_MB, color="#808080", alpha=0.15, zorder=0)
ax.axvspan(LD_START / 1e6 + GAP_MB, LD_END / 1e6, color="#808080", alpha=0.15, zorder=0)

ax.bar(rep_x, rep_win_df["total"], width=WIN_SIZE / 1e6 * 0.9,
       color="black", edgecolor="none", alpha=0.85)

ax.axhline(total_other, color="black", linestyle="--", linewidth=1.2,
           label="Genome-wide density")

ax.set_xlabel(f"Position on {CHROM} (Mb)", fontsize=16)
ax.set_ylabel("Total repeat content (% bases)", fontsize=16)
ax.set_xlim(0, chrom6_len / 1e6)
ax.tick_params(axis="both", labelsize=13)
ax.legend(fontsize=13, loc="upper right")
plt.tight_layout()
out2 = f"{OUT_DIR}/repeat_total_DEINVI00006_sliding_window.jpg"
plt.savefig(out2, dpi=150)
plt.close()
print(f"  Saved: {out2}")

# ── Plot: regional comparison grouped bar chart ───────────────────────────────
print("Generating plot: regional repeat class comparison...")
region_labels = [f"INV1\n({LD_START/1e6:.1f}-{LD_END/1e6:.1f} Mb)",
                  f"INV2\n({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f} Mb)",
                  f"Background\n(rest of {CHROM})", "Rest of\ngenome"]
fracs = [ld_frac, inv2_frac, bg_frac, other_frac]

shown = [c for c in CLASS_ORDER if max(f[c] for f in fracs) > 0.5]
x_pos = np.arange(len(shown))
width = 0.2

fig, ax = plt.subplots(figsize=(14, 5))
region_bar_colors = ["#e07b54", "#9b59b6", "#4a8db5", "#6aab6a"]
for i, (label, frac) in enumerate(zip(region_labels, fracs)):
    vals = [frac[c] for c in shown]
    ax.bar(x_pos + i * width, vals, width, label=label,
           color=region_bar_colors[i], alpha=0.85, edgecolor="white")

ax.set_xticks(x_pos + 1.5 * width)
ax.set_xticklabels(shown, rotation=30, ha="right", fontsize=10)
ax.set_ylabel("% bases covered", fontsize=12)
ax.set_title("Repeat class composition: INV1 vs INV2 vs background vs genome", fontsize=13)
ax.legend(fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
out3 = f"{OUT_DIR}/repeat_regional_comparison.jpg"
plt.savefig(out3, dpi=150)
plt.close()
print(f"  Saved: {out3}")

# ── Plot: LTR divergence distributions (element age) ─────────────────────────
print("Generating plot: LTR element divergence distributions...")
ltr_classes = ["LTR/Gypsy", "LTR/Copia", "LTR/Other"]
ltr_ld = chrom6_rm[in_ld & chrom6_rm["broad"].isin(ltr_classes)]
ltr_inv2 = chrom6_rm[in_inv2 & chrom6_rm["broad"].isin(ltr_classes)]
ltr_bg = chrom6_rm[in_bg & chrom6_rm["broad"].isin(ltr_classes)]
ltr_other = rm_all[(rm_all["chrom"] != CHROM) & rm_all["broad"].isin(ltr_classes)]

fig, axes = plt.subplots(1, 4, figsize=(18, 4), sharey=False)
bins = np.arange(0, 55, 2)

for ax, (data, label, color) in zip(axes, [
    (ltr_ld, f"INV1 ({LD_START/1e6:.1f}-{LD_END/1e6:.1f} Mb)", "#e07b54"),
    (ltr_inv2, f"INV2 ({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f} Mb)", "#9b59b6"),
    (ltr_bg, f"Background ({CHROM})", "#4a8db5"),
    (ltr_other, "Rest of genome", "#6aab6a"),
]):
    for cls, ls in zip(ltr_classes, ["-", "--", ":"]):
        subset = data[data["broad"] == cls]["div"]
        if len(subset) > 10:
            ax.hist(subset, bins=bins, density=True, histtype="step",
                    linewidth=1.5, linestyle=ls, color=CLASS_COLORS[cls], label=cls)
    med = data["div"].median() if len(data) > 0 else 0
    ax.axvline(med, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_title(f"{label}\n(n={len(data):,}, median div={med:.1f}%)", fontsize=9)
    ax.set_xlabel("% divergence from consensus", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.legend(fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)

plt.suptitle("LTR element age (% divergence) by region", fontsize=12, y=1.02)
plt.tight_layout()
out4 = f"{OUT_DIR}/LTR_divergence_by_region.jpg"
plt.savefig(out4, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {out4}")

# ── Statistical tests ─────────────────────────────────────────────────────────
print("\nStatistical comparisons (Mann-Whitney U, window-level total repeat %):")
ld_win = rep_win_df[rep_win_df["region"] == "INV1"]["total"]
inv2_win = rep_win_df[rep_win_df["region"] == "INV2"]["total"]
bg_win = rep_win_df[rep_win_df["region"] == "background"]["total"]
non_ld_win = rep_win_df[rep_win_df["region"] != "INV1"]["total"]   # old (contaminated) comparison, kept for reference

u, p = stats.mannwhitneyu(ld_win, bg_win, alternative="two-sided")
print(f"  INV1:                mean={ld_win.mean():.1f}%, n={len(ld_win)}")
print(f"  True background:     mean={bg_win.mean():.1f}%, n={len(bg_win)}")
print(f"  U={u:.0f}, p={p:.3e}   <-- corrected primary comparison")

u_old, p_old = stats.mannwhitneyu(ld_win, non_ld_win, alternative="two-sided")
print(f"  [Old/contaminated] Non-INV1 (incl. INV2): mean={non_ld_win.mean():.1f}%, n={len(non_ld_win)}")
print(f"  [Old/contaminated] U={u_old:.0f}, p={p_old:.3e}")

u_inv2, p_inv2 = stats.mannwhitneyu(inv2_win, bg_win, alternative="two-sided")
print(f"  INV2:                 mean={inv2_win.mean():.1f}%, n={len(inv2_win)}")
print(f"  INV2 vs background: U={u_inv2:.0f}, p={p_inv2:.3e}")

print("\nLTR divergence (INV1 vs true background):")
u2, p2 = stats.mannwhitneyu(ltr_ld["div"], ltr_bg["div"], alternative="two-sided")
print(f"  INV1 LTR median div:       {ltr_ld['div'].median():.1f}%  (n={len(ltr_ld):,})")
print(f"  Background LTR median div: {ltr_bg['div'].median():.1f}%  (n={len(ltr_bg):,})")
print(f"  U={u2:.0f}, p={p2:.3e}")

repeat_class_summary.to_csv(f"{OUT_DIR}/repeat_class_summary.tsv", sep="\t")
print(f"\nClass summary saved to {OUT_DIR}/repeat_class_summary.tsv")

# ── Box-whisker: window-level total repeat %, INV1 vs INV2 vs background ─────
print("\nGenerating box-whisker: window-level total repeat content...", flush=True)

rng = np.random.default_rng(0)
box_color = "#808080"


def box_whisker_plot(vals_lst, box_labels, genome_ref, ylabel, out_path, u_stat, p_stat, show_legend=True):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    positions = list(range(len(vals_lst)))
    bp = ax.boxplot([v.values for v in vals_lst], positions=positions, widths=0.5, showfliers=False,
                     patch_artist=True, medianprops=dict(color="black", linewidth=1.5))
    for patch in bp["boxes"]:
        patch.set_facecolor(box_color)
        patch.set_alpha(0.35)
    for i, vals in enumerate(vals_lst):
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals.values, s=18,
                   color=box_color, alpha=0.7, edgecolor="white", linewidth=0.4)
    ax.axhline(genome_ref[1], color="black", linestyle="--", linewidth=1.2, label=genome_ref[0])
    ax.set_xticks(positions)
    ax.set_xticklabels(box_labels, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=12)
    if show_legend:
        ax.legend(fontsize=9, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


box_whisker_plot(
    [ld_win, inv2_win, bg_win],
    [f"43-49.3Mb region\nn={len(ld_win)} windows",
     f"13-43Mb region\nn={len(inv2_win)} windows",
     f"Background\n(rest of {CHROM})\nn={len(bg_win)} windows"],
    (f"Genome-wide avg (excl. {CHROM}): {total_other:.1f}%", total_other),
    "Total repeat content (% bases), per 500 kb window",
    f"{OUT_DIR}/repeat_total_boxplot.jpg",
    u, p,
)

# ── LINE-specific analysis -- LINE content is ~3.5x enriched in INV1
#    (3.2% vs 0.9% background / 1.1% rest-of-genome), a much larger
#    fold-difference than any other repeat class ─────────────────────────────
print("\nGenerating LINE-specific plots...", flush=True)

fig, ax = plt.subplots(figsize=(14, 4))
ax.axvspan(INV2_START / 1e6, LD_START / 1e6 - GAP_MB, color="#808080", alpha=0.15, zorder=0)
ax.axvspan(LD_START / 1e6 + GAP_MB, LD_END / 1e6, color="#808080", alpha=0.15, zorder=0)
ax.bar(rep_x, rep_win_df["LINE"], width=WIN_SIZE / 1e6 * 0.9,
       color="black", edgecolor="none", alpha=0.85)
line_other = other_frac["LINE"]
ax.axhline(line_other, color="black", linestyle="--", linewidth=1.2,
           label="Genome-wide density")

ymax = rep_win_df["LINE"].max()

ax.set_xlabel(f"Position on {CHROM} (Mb)", fontsize=16)
ax.set_ylabel("LINE content (% bases)", fontsize=16)
ax.set_xlim(0, chrom6_len / 1e6)
ax.set_ylim(top=ymax * 1.15)
ax.tick_params(axis="both", labelsize=13)
ax.legend(fontsize=13, loc="upper right")
plt.tight_layout()
out_line1 = f"{OUT_DIR}/LINE_content_DEINVI00006_sliding_window.jpg"
plt.savefig(out_line1, dpi=150)
plt.close()
print(f"  Saved: {out_line1}")

# Box-whisker: window-level LINE %, INV1 vs INV2 vs background
ld_line_win = rep_win_df[rep_win_df["region"] == "INV1"]["LINE"]
inv2_line_win = rep_win_df[rep_win_df["region"] == "INV2"]["LINE"]
bg_line_win = rep_win_df[rep_win_df["region"] == "background"]["LINE"]
u_line, p_line = stats.mannwhitneyu(ld_line_win, bg_line_win, alternative="two-sided")
print(f"  INV1 LINE:       mean={ld_line_win.mean():.2f}%, n={len(ld_line_win)}")
print(f"  Background LINE: mean={bg_line_win.mean():.2f}%, n={len(bg_line_win)}")
print(f"  U={u_line:.0f}, p={p_line:.3e}")

box_whisker_plot(
    [ld_line_win, inv2_line_win, bg_line_win],
    [f"43-49.3Mb region\nn={len(ld_line_win)} windows",
     f"13-43Mb region\nn={len(inv2_line_win)} windows",
     f"Background\n(rest of {CHROM})\nn={len(bg_line_win)} windows"],
    (f"Genome-wide avg (excl. {CHROM}): {line_other:.2f}%", line_other),
    "% LINE content per 500kb window",
    f"{OUT_DIR}/LINE_content_boxplot.jpg",
    u_line, p_line,
    show_legend=False,
)

# ── Satellite-specific analysis -- checked because gene-poor, TE-enriched
#    regions are sometimes (peri)centromeric; satellite content here is
#    actually LOWER than background, arguing against that ────────────────────
print("\nGenerating satellite-specific plots...", flush=True)

fig, ax = plt.subplots(figsize=(14, 4))
ax.axvspan(INV2_START / 1e6, LD_START / 1e6 - GAP_MB, color="#808080", alpha=0.15, zorder=0)
ax.axvspan(LD_START / 1e6 + GAP_MB, LD_END / 1e6, color="#808080", alpha=0.15, zorder=0)
ax.bar(rep_x, rep_win_df["Satellite"], width=WIN_SIZE / 1e6 * 0.9,
       color="black", edgecolor="none", alpha=0.85)
sat_other = other_frac["Satellite"]
ax.axhline(sat_other, color="black", linestyle="--", linewidth=1.2,
           label="Genome-wide density")

ax.set_xlabel(f"Position on {CHROM} (Mb)", fontsize=16)
ax.set_ylabel("Satellite content (% bases)", fontsize=16)
ax.set_xlim(0, chrom6_len / 1e6)
ax.tick_params(axis="both", labelsize=13)
ax.legend(fontsize=13, loc="upper right")
plt.tight_layout()
out_sat1 = f"{OUT_DIR}/satellite_content_DEINVI00006_sliding_window.jpg"
plt.savefig(out_sat1, dpi=150)
plt.close()
print(f"  Saved: {out_sat1}")

# Box-whisker: window-level satellite %, INV1 vs INV2 vs background
ld_sat_win = rep_win_df[rep_win_df["region"] == "INV1"]["Satellite"]
inv2_sat_win = rep_win_df[rep_win_df["region"] == "INV2"]["Satellite"]
bg_sat_win = rep_win_df[rep_win_df["region"] == "background"]["Satellite"]
u_sat, p_sat = stats.mannwhitneyu(ld_sat_win, bg_sat_win, alternative="two-sided")
print(f"  INV1 Satellite:       mean={ld_sat_win.mean():.3f}%, n={len(ld_sat_win)}")
print(f"  Background Satellite: mean={bg_sat_win.mean():.3f}%, n={len(bg_sat_win)}")
print(f"  U={u_sat:.0f}, p={p_sat:.3e}")

box_whisker_plot(
    [ld_sat_win, inv2_sat_win, bg_sat_win],
    [f"43-49.3Mb region\nn={len(ld_sat_win)} windows",
     f"13-43Mb region\nn={len(inv2_sat_win)} windows",
     f"Background\n(rest of {CHROM})\nn={len(bg_sat_win)} windows"],
    (f"Genome-wide avg (excl. {CHROM}): {sat_other:.3f}%", sat_other),
    "% Satellite content per 500kb window",
    f"{OUT_DIR}/satellite_content_boxplot.jpg",
    u_sat, p_sat,
)

# ── Window tables for the stacked DEINVI_00006 figure (DEINVI_0006_stacked.py) ──
gd_win_df.to_csv(f"{OUT_DIR}/gene_density_windows_{CHROM}.tsv", sep="\t", index=False)
rep_win_df.to_csv(f"{OUT_DIR}/repeat_windows_{CHROM}.tsv", sep="\t", index=False)
pd.DataFrame({"measure": ["gene_density", "total_repeat"] + [f"repeat_{c}" for c in CLASS_ORDER],
              "genome_wide": [other_dens, total_other] + [other_frac[c] for c in CLASS_ORDER]}
             ).to_csv(f"{OUT_DIR}/genome_wide_refs.tsv", sep="\t", index=False)
print(f"\nWindow tables saved to {OUT_DIR}/{{gene_density,repeat}}_windows_{CHROM}.tsv")


# ════════════════════════════════════════════════════════════════════════════
# 4. Fine-scale boundary analysis
# ════════════════════════════════════════════════════════════════════════════
# Looks for sharp compositional transitions and shared TE families at both
# breakpoints (~43 Mb and ~49.3 Mb).
#
# NOTE: DEINVI_00006 carries a SECOND, adjacent putative inversion (INV2,
# ~13.0-43.0 Mb) immediately upstream of the region of interest (INV1). The
# left boundary analyzed here (~43 Mb) is therefore the shared breakpoint
# between INV2 and INV1, not a boundary between INV1 and neutral background.
# The "outside"/flanking reference used for the left-boundary TE family
# enrichment calculation below is INV2 interior sequence, not background --
# see the caveat printed at that step.

BOUNDARY_WIN_SIZE = 50_000    # 50 kb fine-scale windows
FLANK = 3_000_000             # +/-3 Mb around each boundary

CLASS_ORDER_BOUNDARY = list(CLASS_COLORS.keys())


def windowed_composition(df, region_start, region_end, win_size):
    rows = []
    for ws in range(region_start, region_end, win_size):
        we = min(ws + win_size, region_end)
        wsz = we - ws
        sub = df[(df["start"] >= ws) & (df["end"] <= we)]
        row = {"mid": (ws + we) / 2}
        for cls in CLASS_ORDER_BOUNDARY:
            row[cls] = sub.loc[sub["broad"] == cls, "length"].sum() / wsz * 100
        row["total"] = sub["length"].sum() / wsz * 100
        rows.append(row)
    return pd.DataFrame(rows)


b1_start, b1_end = LD_START - FLANK, LD_START + FLANK
b2_start, b2_end = LD_END - FLANK, LD_END + FLANK

print(f"Computing fine-scale windows around boundary 1 (~{LD_START/1e6:.1f} Mb)...", flush=True)
win_b1 = windowed_composition(chrom6_rm, b1_start, b1_end, BOUNDARY_WIN_SIZE)
print(f"Computing fine-scale windows around boundary 2 (~{LD_END/1e6:.1f} Mb)...", flush=True)
win_b2 = windowed_composition(chrom6_rm, b2_start, b2_end, BOUNDARY_WIN_SIZE)


# ── Plot: stacked composition + total repeat at each boundary ────────────────
def plot_boundary(win_df, boundary_pos, boundary_name, filename):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1]})
    x = win_df["mid"].values / 1e6
    xlo, xhi = x.min(), x.max()

    inv2_lo = max(INV2_START / 1e6, xlo)
    inv2_hi = min(INV2_END / 1e6, xhi)
    if inv2_lo < inv2_hi:
        for ax in (ax1, ax2):
            ax.axvspan(inv2_lo, inv2_hi, color="#9b59b6", alpha=0.10, zorder=0,
                       label="INV2 (second putative inversion)" if ax is ax1 else None)

    bottom = np.zeros(len(win_df))
    for cls in CLASS_ORDER_BOUNDARY:
        y = win_df[cls].values
        ax1.bar(x, y, bottom=bottom, width=BOUNDARY_WIN_SIZE / 1e6 * 0.9,
                color=CLASS_COLORS[cls], label=cls, edgecolor="none", alpha=0.9)
        bottom += y

    ax1.axvline(boundary_pos / 1e6, color="black", linewidth=1.5,
                linestyle="--", label=f"LD boundary ({boundary_pos/1e6:.1f} Mb)")
    ax1.set_ylabel("Repeat content (% bases)", fontsize=11)
    ax1.set_title(f"Repeat composition at LD block boundary: {boundary_name}", fontsize=12)
    ax1.legend(fontsize=7, ncol=4, loc="upper right")

    ax2.bar(x, win_df["total"], width=BOUNDARY_WIN_SIZE / 1e6 * 0.9,
            color="#888888", edgecolor="none", alpha=0.7)
    ax2.axvline(boundary_pos / 1e6, color="black", linewidth=1.5, linestyle="--")
    ax2.set_ylabel("Total\nrepeat %", fontsize=10)
    ax2.set_xlabel("Position on DEINVI_00006 (Mb)", fontsize=11)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"  Saved: {filename}")


print("Generating boundary plots...")
plot_boundary(win_b1, LD_START, f"~{LD_START/1e6:.1f} Mb (left boundary)",
              f"{OUT_DIR}/boundary_left_{LD_START/1e6:.1f}Mb.jpg")
plot_boundary(win_b2, LD_END, f"~{LD_END/1e6:.1f} Mb (right boundary)",
              f"{OUT_DIR}/boundary_right_{LD_END/1e6:.1f}Mb.jpg")

# ── TE family enrichment at each boundary ─────────────────────────────────────
# Compare TE families in a narrow window AT the boundary vs. flanking regions
# "at boundary": +/-500 kb of the breakpoint
# "inside LD":   500 kb-2 Mb inside the block
# "outside LD":  500 kb-2 Mb outside the block
NARROW = 500_000
WIDE = 2_000_000


def boundary_families(df, bp, label):
    at_bp = df[(df["start"] >= bp - NARROW) & (df["end"] <= bp + NARROW)]
    inside = df[(df["start"] >= bp + NARROW) & (df["end"] <= bp + WIDE)] if bp == LD_START else \
        df[(df["start"] >= bp - WIDE) & (df["end"] <= bp - NARROW)]
    outside = df[(df["start"] >= bp - WIDE) & (df["end"] <= bp - NARROW)] if bp == LD_START else \
        df[(df["start"] >= bp + NARROW) & (df["end"] <= bp + WIDE)]

    def family_density(sub, size_bp):
        return (sub.groupby("name")["length"].sum() / size_bp * 100).sort_values(ascending=False)

    at_dens = family_density(at_bp, 2 * NARROW)
    in_dens = family_density(inside, WIDE - NARROW)
    out_dens = family_density(outside, WIDE - NARROW)

    background = (in_dens.add(out_dens, fill_value=0)) / 2
    enrichment = (at_dens / background.reindex(at_dens.index, fill_value=0.01)).fillna(0)
    enrichment = enrichment[enrichment > 1].sort_values(ascending=False).head(30)

    result = pd.DataFrame({
        "at_boundary_%": at_dens.reindex(enrichment.index, fill_value=0).round(3),
        "inside_%":      in_dens.reindex(enrichment.index, fill_value=0).round(3),
        "outside_%":     out_dens.reindex(enrichment.index, fill_value=0).round(3),
        "enrichment":    enrichment.round(2),
        "class":         [df.loc[df["name"] == n, "class"].iloc[0] if (df["name"] == n).any() else "?"
                          for n in enrichment.index],
    })
    print(f"\nTop enriched TE families at {label} (+/-{NARROW//1000} kb vs flanks):")
    if bp == LD_START:
        print(f"  CAVEAT: 'outside'/flanking reference here ({(bp-WIDE)/1e6:.1f}-{(bp-NARROW)/1e6:.1f} Mb) "
              f"falls inside INV2 ({INV2_START/1e6:.1f}-{INV2_END/1e6:.1f} Mb), not neutral background. "
              f"Enrichment values below are relative to INV2 interior sequence, not genome background.")
    print(result.to_string())
    return result


print("\n-- TE family enrichment at boundaries --")
fam_b1 = boundary_families(chrom6_rm, LD_START, f"left boundary (~{LD_START/1e6:.1f} Mb)")
fam_b2 = boundary_families(chrom6_rm, LD_END, f"right boundary (~{LD_END/1e6:.1f} Mb)")

# ── Shared enriched families at BOTH boundaries ──────────────────────────────
shared = set(fam_b1.index) & set(fam_b2.index)
print(f"\n-- Families enriched at BOTH boundaries ({len(shared)} found) --")
if shared:
    shared_df = pd.DataFrame({
        "class":       [fam_b1.loc[f, "class"] for f in shared],
        "enrich_b1":   [fam_b1.loc[f, "enrichment"] for f in shared],
        "enrich_b2":   [fam_b2.loc[f, "enrichment"] for f in shared],
        "at_b1_%":     [fam_b1.loc[f, "at_boundary_%"] for f in shared],
        "at_b2_%":     [fam_b2.loc[f, "at_boundary_%"] for f in shared],
    }, index=list(shared)).sort_values("enrich_b1", ascending=False)
    print(shared_df.to_string())
    shared_df.to_csv(f"{OUT_DIR}/shared_boundary_families.tsv", sep="\t")
else:
    print("  None found.")

# ── Side-by-side LTR Copia composition at both boundaries ────────────────────
print("\nGenerating side-by-side LTR/Copia plot at both boundaries...")
fig, axes = plt.subplots(1, 2, figsize=(14, 4), sharey=True)
for ax, (win_df, bp, label) in zip(axes, [
    (win_b1, LD_START, f"Left boundary (~{LD_START/1e6:.1f} Mb)"),
    (win_b2, LD_END, f"Right boundary (~{LD_END/1e6:.1f} Mb)"),
]):
    x = win_df["mid"].values / 1e6
    for cls in ["LTR/Gypsy", "LTR/Copia", "LTR/Other", "Unknown"]:
        ax.plot(x, win_df[cls], label=cls, color=CLASS_COLORS[cls], linewidth=1.5)
    ax.axvline(bp / 1e6, color="black", linewidth=1.5, linestyle="--",
               label=f"Boundary ({bp/1e6:.0f} Mb)")
    ax.set_title(label, fontsize=11)
    ax.set_xlabel("Position on DEINVI_00006 (Mb)", fontsize=10)
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

axes[0].set_ylabel("% bases covered", fontsize=11)
plt.suptitle("LTR class composition at LD block boundaries (50 kb windows)", fontsize=12)
plt.tight_layout()
out5 = f"{OUT_DIR}/boundary_LTR_comparison.jpg"
plt.savefig(out5, dpi=150)
plt.close()
print(f"  Saved: {out5}")

# ── Zoomed-in, family-resolved view of the boundary enrichment signal ────────
# The enrichment table above is real (20-200x), but it is a LOCAL, narrow
# (+/-500kb), FAMILY-specific signal -- diluted away in the boundary plot
# because that plot averages whole broad classes (Gypsy/Copia/Other all
# lumped) over a wide 50kb x +/-3Mb view, where genome-typical background
# from those same classes elsewhere in the flank swamps the local spike.
# This plot narrows the window (10 kb) and flank (+/-1 Mb) and adds a rug of
# the exact element positions for the top enriched families, so the local
# cluster is visible directly rather than inferred from an averaged
# composition curve.
print("\nGenerating zoomed-in boundary enrichment plot (10 kb windows, +/-1 Mb, top families)...")

ZOOM_WIN = 10_000
ZOOM_FLANK = 1_000_000
TOP_N_FAM = 5


def zoomed_boundary_plot(df, bp, fam_table, label, filename):
    z_start, z_end = bp - ZOOM_FLANK, bp + ZOOM_FLANK
    zwin = windowed_composition(df, z_start, z_end, ZOOM_WIN)
    ltr_total = zwin[["LTR/Gypsy", "LTR/Copia", "LTR/Other"]].sum(axis=1)

    top_fams = fam_table.head(TOP_N_FAM).index.tolist()
    fam_colors = plt.cm.tab10(np.linspace(0, 1, TOP_N_FAM))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                    gridspec_kw={"height_ratios": [2, 1.3]})

    xlo, xhi = z_start / 1e6, z_end / 1e6
    inv2_lo = max(INV2_START / 1e6, xlo)
    inv2_hi = min(INV2_END / 1e6, xhi)
    if inv2_lo < inv2_hi:
        ax1.axvspan(inv2_lo, inv2_hi, color="#9b59b6", alpha=0.10, zorder=0, label="INV2")
        ax2.axvspan(inv2_lo, inv2_hi, color="#9b59b6", alpha=0.10, zorder=0)

    ax1.plot(zwin["mid"].values / 1e6, ltr_total.values, color="#c0392b",
              linewidth=1.8, label="Total LTR (Gypsy+Copia+Other)")
    ax1.axvline(bp / 1e6, color="black", linewidth=1.5, linestyle="--",
                label=f"Boundary ({bp/1e6:.2f} Mb)")
    ax1.set_ylabel("LTR content\n(% bases, 10 kb windows)", fontsize=10)
    ax1.set_title(f"Zoomed boundary view: {label}", fontsize=12)
    ax1.legend(fontsize=8, loc="upper right")
    ax1.spines[["top", "right"]].set_visible(False)

    zoom_df = df[(df["start"] >= z_start) & (df["end"] <= z_end)]
    for i, fam in enumerate(top_fams):
        fam_hits = zoom_df[zoom_df["name"] == fam]
        mid = (fam_hits["start"] + fam_hits["end"]) / 2 / 1e6
        ax2.scatter(mid, np.full(len(mid), i), s=40, color=fam_colors[i],
                    label=f"{fam} ({fam_table.loc[fam, 'class']}, {fam_table.loc[fam,'enrichment']:.0f}x)")
    ax2.axvline(bp / 1e6, color="black", linewidth=1.5, linestyle="--")
    ax2.set_yticks(range(TOP_N_FAM))
    ax2.set_yticklabels([])
    ax2.set_ylabel(f"Top {TOP_N_FAM}\nenriched families", fontsize=10)
    ax2.set_xlabel(f"Position on {CHROM} (Mb)", fontsize=11)
    ax2.legend(fontsize=7, loc="upper right", ncol=1)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.set_ylim(-0.7, TOP_N_FAM - 0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"  Saved: {filename}")


zoomed_boundary_plot(chrom6_rm, LD_START, fam_b1, f"~{LD_START/1e6:.1f} Mb (left boundary)",
                      f"{OUT_DIR}/boundary_left_{LD_START/1e6:.1f}Mb_zoomed.jpg")
zoomed_boundary_plot(chrom6_rm, LD_END, fam_b2, f"~{LD_END/1e6:.1f} Mb (right boundary)",
                      f"{OUT_DIR}/boundary_right_{LD_END/1e6:.1f}Mb_zoomed.jpg")

print("Done.")
