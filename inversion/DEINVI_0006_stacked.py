#!/usr/bin/env python3
"""
Stacked, x-aligned summary figure for putative inversions on one scaffold
(originally DEINVI_00006; per-scaffold settings live in SCAFFOLDS below).

Usage: DEINVI_0006_stacked.py [scaffold]     (default DEINVI_00006)

One column of panels sharing a single x-axis (position on the scaffold):
  1. fog GEMMA Manhattan (SNPs, Bonferroni threshold, gene labels exported by
     gvtp_popgen/selection/gemma/annotated_manhattan.R, so they match its
     per-scaffold Manhattan plot)
  2. Gene density, 500 kb windows    } window tables written by
  3. LINE content, 500 kb windows    } scaffold_windows.py
  4. % variance explained by PC1 per window, from the same run as panel 5
  5. winPCA PC1, one line per sample
  6. winPCA SNP heterozygosity -- GT-based (winpca only computes hetp from
     called genotypes), no MAF filter, >=8x samples only
Panels 4-6 are left blank for a scaffold with no winPCA run yet.

For DEINVI_00006: PC1 is the all-locality 5% MAF GL run
(winPCA_maf5_DEINVI_00006.sh, its "_flipped" copy with auto-polarization sign
flips corrected); heterozygosity is the all-locality >=8x GT run (_all8x).
For DEINVI_00003 / 00008 / 00017: all-locality 5% MAF runs from
winPCA_maf5.sh (+ the all-locality GT run for 00017's heterozygosity).

winPCA lines are colored by fog with winpca's own numeric color scale (see
small_chromplot.py). metadata.tsv stores fog NEGATED relative to
env_var.txt (exactly -1x; presumably so winpca's navy->yellow scale puts
high fog at navy), so it is flipped back here to real fog (days per month)
and the color scale reversed to match -- each sample keeps the same color
as in the standalone winPCA plots, but the colorbar reads true fog values.

Lines are drawn in descending fog, so the lowest-fog samples sit on top.
Strauss is 277 of the 339 all-locality samples, so for plotting only (the
PCA itself used everyone) each locality is capped at PER_SITE_CAP randomly
chosen individuals (fixed seed). Samples present in both the PC1 and
heterozygosity runs are picked first, so the two panels show the same
individuals wherever possible; a site with fewer such samples than the cap
is topped up from the PC1-only ones.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
import temporal_null as tn      # centromere-matched random-SNP null (shared with fog_peaks_temporal.py)

WINPCA_DIR = "/private/groups/meyerlab/esteiger/inversions/all_localities"
OUT_DIR    = "/private/groups/meyerlab/esteiger/inversions"
# Per-scaffold settings. shade: grey boxes (Mb; None = the scaffold's centromere
# call from call_centromeres.py, CENTROMERES); lines: draw the per-sample winPCA
# PC1 and heterozygosity panels (False = only the PC1 % variance track; for
# chr 6 those panels are in the supplementary figure, chr6_synteny_panels.py);
# winpca / hetp: winPCA output
# prefixes (None = panel left blank); labels: genes to annotate on the
# Manhattan panel (None = every exported label of <= MAX_LABEL_CHARS
# characters, labels at the same SNP joined with "/").
SCAFFOLDS = {
    "DEINVI_00006": dict(
        shade=None,
        lines=False,
        winpca=f"{WINPCA_DIR}/DEINVI_00006_GL_maf5_flipped",
        hetp=f"{WINPCA_DIR}/DEINVI_00006_GT_all8x",
        labels=["GNS1"],
        label_text={"GNS1": "GNS1, SNF1"},   # display text for a label (the peak spans both genes)
        synteny=False,             # top panel: ancestral segments (see SYNTENY_* below); off --
                                   # synteny is in the supplementary figure (chr6_synteny_panels.py)
        bottom_row=True,           # bottom panels (chr 6 fog peak data; see BOTTOM_* below)
        out=f"{OUT_DIR}/DEINVI_0006.jpg"),
    # 5% MAF runs from winPCA_maf5.sh (submitted 2026-10-03); switch to the
    # "_flipped" copies once sign flips have been checked/fixed
    "DEINVI_00003": dict(shade=None, winpca=f"{WINPCA_DIR}/DEINVI_00003_GL_maf5", hetp=None,
                         labels=None, out=f"{OUT_DIR}/DEINVI_0003.jpg"),
    "DEINVI_00008": dict(shade=None, winpca=f"{WINPCA_DIR}/DEINVI_00008_GL_maf5", hetp=None,
                         labels=None, out=f"{OUT_DIR}/DEINVI_0008.jpg"),
    "DEINVI_00017": dict(shade=None, winpca=f"{WINPCA_DIR}/DEINVI_00017_GL_maf5",
                         hetp=f"{WINPCA_DIR}/DEINVI_00017_GT_all8x",
                         labels=None, out=f"{OUT_DIR}/DEINVI_0017.jpg"),
}
CHROM = sys.argv[1] if len(sys.argv) > 1 else "DEINVI_00006"
CFG = SCAFFOLDS[CHROM]
CHROM_NUM = int(CHROM.split("_")[1])
CENTROMERES = "/private/groups/meyerlab/esteiger/inversions/centromere/centromeres.tsv"
if CFG["shade"] is None:
    _cen = pd.read_csv(CENTROMERES, sep="\t").query("confident").set_index("chr")
    CFG["shade"] = [(_cen.loc[CHROM, "cen_start"] / 1e6, _cen.loc[CHROM, "cen_end"] / 1e6)] if CHROM in _cen.index else []
lines = CFG.get("lines", True)
SHADE_ALPHA   = 0.25           # grey box opacity
MAX_LABEL_CHARS = 12

GENE_DENSITY_DIR = "/private/groups/meyerlab/esteiger/inversions/gene_density"
GEMMA_DIR        = "/private/groups/meyerlab/esteiger/selection/gemma"
MANHATTAN_VAR    = "fog"
WINPCA_META      = "/private/groups/meyerlab/esteiger/inversions/metadata.tsv"
WINPCA_COLOR_BY  = "fog"

# winpca's numeric color scale, reversed because fog is un-negated below
WINPCA_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "winpca_r", list(zip([0, 0.25, 0.5, 0.75, 1], ["#FFFF00", "#FF0000", "#800080", "#0000FF", "#000080"])))

LABEL_SIZE       = 24         # every axis title, the Manhattan gene label and the colorbar label
TICK_SIZE        = 19         # every tick label and legend entry
PANEL_LETTERS    = "CDEFGHIJKLM"  # panel labels: column panels top to bottom, then the bottom row left to right
PANEL_LETTER_SIZE = 28

# Bottom row (panels G-I), chr 6 only; all Tranquillon Mtn (Strauss) data, drawn
# in its locality colour (TRANQ_COLOR, as in the Tajima's D figures):
#   G: unit long-term fog vs favorable-allele frequency, every fog peak (chr 6 =
#      diamonds, area ~ plants per unit, dashed line; other peaks = lines only)
#   H: germination-year (June) fog by germination year (one value per plant)
#   I: germination-year fog vs unit-adjusted favorable-allele frequency per
#      cohort, every fog peak (chr 6 = diamonds +-95% CI, labelled by year,
#      dashed within-unit line; other peaks = lines only)
# Data exported by gvtp_popgen/selection/gemma/fog_peak_figures.py (G, I; same
# plants, polarization and spatial units as fog_hypothesis_tests.py) and
# fog_deinvi6_figure.R (H).
BOTTOM_BY_YEAR = "/private/groups/meyerlab/esteiger/selection/gemma/fog/fog_deinvi6_by_year.tsv"
BOTTOM_UNITS   = "/private/groups/meyerlab/esteiger/selection/gemma/fog/fog_peaks_by_unit.tsv"
BOTTOM_PEAK_YR = "/private/groups/meyerlab/esteiger/selection/gemma/fog/fog_peaks_by_year.tsv"
TRANQ_COLOR    = "#009e73"
OTHER_ALPHA    = 0.45         # non-chr 6 peak lines in G and I
BOTTOM_WIDTHS  = [1.45, 0.65, 0.85]   # relative widths of panels G, H, I
# panel I year labels: (x, y offset in points, alignment), placed clear of the CIs and line
YEAR_LABEL_OFFSET = {2020: (12, -16, "left"), 2021: (-14, 14, "right"), 2022: (14, -16, "left")}

# Synteny panel (top; chr 6 only): black = stretches syntenic with sunflower
# (Helianthus annuus) and lettuce (Lactuca sativa), i.e. runs of GENESPACE
# syntenic gene pairs (blkID assigned) split wherever consecutive hits are
# > SYNTENY_GAP Mb apart. Hatched = SYNTENY_REVERSED, a run of genes whose
# order is reversed relative to sunflower, H. anomalus and lettuce (from the
# gene order of the syntenic hits; minimal extent, outermost reversed genes).
SYNTENY_DIR      = "/private/groups/meyerlab/esteiger/inversions/synteny/genespace/syntenicHits"
SYNTENY_SPECIES  = ["Helianthus_annuus", "Lactuca_sativa"]
SYNTENY_GAP      = 5.0
SYNTENY_REVERSED = [(19.34, 20.83)]

# Temporal checks (printed and written to TEMPORAL_OUT): is the favorable
# haplotype's frequency in the foggiest germination cohort unusual?
# Statistic (fixed in advance): allele frequency in FOGGY_YEAR minus the
# pooled other years, from mean dosage / 2.
#   1. Random-SNP null: the same statistic for random SNPs elsewhere in the
#      genome (one per random 2 kb region, not on DEINVI_00006), each polarized
#      exactly like the outliers (favorable = allele whose coverage-adjusted
#      dosage increases with long-term fog across ALL samples, as in
#      germ_year_polygenic.R / fog_deinvi6_figure.R), with favorable-allele
#      frequency in Strauss matched to the peak haplotype (+-NULL_FREQ_TOL),
#      same individuals and years. Imputed dosages (DS) from the GWAS VCF for
#      all samples; cached in NULL_CACHE (shared with fog_peaks_temporal.py).
#   2. Within-stand permutation: years shuffled within each stand (OHANA_order,
#      field_data.txt), so stand composition per year is held fixed.
#   3. Centromere-matched null: as 1., but the random SNPs lie at a similar
#      distance from their own scaffold's centromere as the peak does from the
#      chr 6 centromere (temporal_null.py; same cache and band as
#      fog_peaks_temporal.py), since low recombination near centromeres can
#      make allele frequencies more variable.
FOGGY_YEAR     = 2021
VCF_DS         = "/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog_MAF_auto_GLimp_clean.vcf.gz"
FAI            = "/private/groups/meyerlab/esteiger/inversions/filtered_scaffolds.fai"
FIELD_DATA     = "/private/groups/meyerlab/esteiger/field_data.txt"
ENV_VAR        = "/private/groups/meyerlab/esteiger/env_var.txt"
NULL_CACHE     = "/private/groups/meyerlab/esteiger/selection/gemma/fog/temporal_null_random_snps_allsamples.tsv.gz"
TEMPORAL_OUT   = "/private/groups/meyerlab/esteiger/selection/gemma/fog/temporal_checks.tsv"
PEAK_SPAN      = ("DEINVI_00006", 49031257, 49146391)   # the chr 6 fog peak (19 retained fog SNPs)
FOG_OUTLIERS   = "/private/groups/meyerlab/esteiger/selection/gemma/fog/outliers_for_annot_40000_3_fog.txt"
N_NULL_REGIONS, NULL_FREQ_TOL, N_PERM = 8000, 0.04, 10000
PER_SITE_CAP     = 15         # max individuals drawn per locality in the winPCA panels (None = all)
PC1_YEARS        = None       # sampling years (metadata Year.x) drawn in the PC1 panel (None = all);
                              # plotting only -- the PCA itself used every sample
SUBSAMPLE_SEED   = 1
LINE_WIDTH       = 1.8        # winPCA panels' per-sample lines
LINE_ALPHA       = 0.45
CBAR_ALPHA       = 0.6        # colorbar opacity (a bit above LINE_ALPHA: overlapping lines read darker)
Y_PAD            = 0.01       # y-axis headroom beyond the most extreme value (fraction of data range)

# ── Load data ────────────────────────────────────────────────────────────────
man_df = pd.read_csv(f"{GEMMA_DIR}/{MANHATTAN_VAR}/manhattan_{CHROM}_{MANHATTAN_VAR}.tsv.gz", sep="\t")
man_labels = pd.read_csv(f"{GEMMA_DIR}/{MANHATTAN_VAR}/manhattan_labels_{CHROM}_{MANHATTAN_VAR}.tsv", sep="\t")
if CFG["labels"] is not None:
    man_labels = man_labels[man_labels["Label"].isin(CFG["labels"])]
elif len(man_labels):
    man_labels = (man_labels[man_labels["Label"].str.len() <= MAX_LABEL_CHARS]
                  .groupby(["Position", "p_final"], as_index=False)["Label"].agg("/".join))
man_df["logp"] = -np.log10(man_df["p_final"])
man_labels["logp"] = -np.log10(man_labels["p_final"])
thresh_logp = -np.log10(man_df["thresh_bonf"].iloc[0])

gd_win_df = pd.read_csv(f"{GENE_DENSITY_DIR}/gene_density_windows_{CHROM}.tsv", sep="\t")
rep_win_df = pd.read_csv(f"{GENE_DENSITY_DIR}/repeat_windows_{CHROM}.tsv", sep="\t")
refs = pd.read_csv(f"{GENE_DENSITY_DIR}/genome_wide_refs_{CHROM}.tsv", sep="\t", index_col=0)["genome_wide"]
CHROM_LEN = man_df["scaf_len"].iloc[0]


def read_winpca(prefix, suffix, what):
    """winPCA output table, or None (panel left blank) if not available yet."""
    path = f"{prefix}.{suffix}.tsv.gz" if prefix else None
    if path is None or not os.path.exists(path):
        print(f"NOTE: no {what} for {CHROM} ({path}) -- panel left blank")
        return None
    return pd.read_csv(path, sep="\t", index_col=0)


pc1 = read_winpca(CFG["winpca"], "pc_1", "winPCA PC1")
pca_stat = read_winpca(CFG["winpca"], "stat", "winPCA stats")
# GT-based heterozygosity comes from a separate (slow) winPCA_GT.sh run
hetp = read_winpca(CFG["hetp"], "hetp", "heterozygosity")
pca_meta = pd.read_csv(WINPCA_META, sep="\t", index_col=0)
pca_samples = pd.Index([])
for df in (pc1, hetp):
    if df is not None:
        pca_samples = pca_samples.union(df.columns)
fog_val = -pd.to_numeric(pca_meta[WINPCA_COLOR_BY], errors="coerce").reindex(pca_samples)
site = pca_meta["Clusters"].reindex(pca_samples)
year = pca_meta["Year.x"].reindex(pca_samples)
rng = np.random.default_rng(SUBSAMPLE_SEED)
present = [set(df.columns) for df in (pc1, hetp) if df is not None]
in_both = set.intersection(*present) if present else set()
drawn_samples = set()
for loc, members in site.groupby(site):
    both = sorted(s for s in members.index if s in in_both)
    rest = sorted(s for s in members.index if s not in in_both)
    if PER_SITE_CAP is None:
        picked = both + rest
    else:
        picked = list(rng.permutation(both))[:PER_SITE_CAP]
        if len(picked) < PER_SITE_CAP:
            picked += list(rng.permutation(rest))[:PER_SITE_CAP - len(picked)]
    drawn_samples.update(picked)
    print(f"  {loc}: {len(picked)} of {len(members)} drawn")

# ── Temporal checks ───────────────────────────────────────────────────────────
def foggy_shift(dosage, years):
    """Allele frequency in FOGGY_YEAR minus pooled other years (dosage / 2)."""
    f = years == FOGGY_YEAR
    return dosage[f].mean() / 2 - dosage[~f].mean() / 2


def random_snp_dosages(samples):
    """Imputed DS dosages for one random SNP per random 2 kb region (not chr 6), cached."""
    if os.path.exists(NULL_CACHE):
        return pd.read_csv(NULL_CACHE, sep="\t", index_col=[0, 1])
    import subprocess, tempfile
    fai = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"])
    fai = fai[fai.chrom != "DEINVI_00006"]
    rng = np.random.default_rng(7)
    chroms = rng.choice(fai.chrom.values, N_NULL_REGIONS, p=fai.len.values / fai.len.sum())
    lens = fai.set_index("chrom")["len"]
    starts = [int(rng.integers(1, lens[c] - 2000)) for c in chroms]
    reg = pd.DataFrame({"chrom": chroms, "start": starts, "end": [x + 2000 for x in starts]}).sort_values(["chrom", "start"])
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        reg.to_csv(fh, sep="\t", header=False, index=False); regfile = fh.name
    out = subprocess.run(["bcftools", "query", "-R", regfile, "-s", ",".join(samples),
                          "-f", "%CHROM\t%POS[\t%DS]\n", VCF_DS], capture_output=True, text=True).stdout.splitlines()
    os.remove(regfile)
    d = pd.DataFrame([l.split("\t") for l in out])
    d = d.drop_duplicates([0, 1]).set_index([0, 1]).astype(float)
    d.index.names = ["chrom", "pos"]; d.columns = samples
    d.to_csv(NULL_CACHE, sep="\t", compression="gzip")
    return d


def temporal_checks(by_year):
    rng = np.random.default_rng(11)
    by_year = by_year.dropna(subset=["pk_score"]).copy()
    by_year["SampleID"] = by_year["SampleID"].astype(str)
    years = by_year["Year"].values
    obs = foggy_shift(by_year["pk_score"].values, years)
    f_peak = by_year["pk_score"].mean() / 2

    # 1. random-SNP null, polarized by long-term fog across all samples
    env = pd.read_csv(ENV_VAR, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    env["Cov"] = pd.to_numeric(env["Cov"], errors="coerce")
    env = env.dropna(subset=["Cov", "fog"])
    env = env[~env.index.isin(tn.BAD_COORD_IDS)]   # untrustworthy-coordinate IDs excluded
    D = random_snp_dosages(list(env.index))
    D.columns = D.columns.astype(str)
    allS = [s for s in D.columns if s in env.index]
    X = np.column_stack([np.ones(len(allS)), env.loc[allS, "Cov"].values])
    H = X @ np.linalg.solve(X.T @ X, X.T)
    fog = env.loc[allS, "fog"].values.astype(float)
    M = D[allS].values.T                                    # samples x SNPs
    beta = (M - H @ M).T @ (fog - H @ fog)
    D = D[by_year["SampleID"].tolist()]
    D = D.where(pd.Series(beta > 0, index=D.index), 2 - D, axis=0)   # favorable dose
    f = D.mean(axis=1) / 2
    D = D[((f - f_peak).abs() <= NULL_FREQ_TOL) & (f > 0) & (f < 1)]
    # thin to one SNP per random 2 kb region (regions are 2 kb windows at random positions)
    region = [f"{c}:{int(p) // 2000}" for c, p in D.index]
    D = D.groupby(region, group_keys=False).apply(lambda x: x.sample(1, random_state=1))
    null = np.array([foggy_shift(row.values, years) for _, row in D.iterrows()])
    p_null = (np.sum(null >= obs) + 1) / (len(null) + 1)

    # 2. within-stand permutation
    field = pd.read_csv(FIELD_DATA, sep=None, engine="python", dtype={"SampleID": str})
    stand = by_year["SampleID"].map(field.drop_duplicates("SampleID").set_index("SampleID")["OHANA_order"]).values
    dos = by_year["pk_score"].values
    groups = [np.where(stand == st)[0] for st in pd.unique(stand)]
    perm = np.empty(N_PERM)
    for k in range(N_PERM):
        y = years.copy()
        for g in groups:
            y[g] = rng.permutation(y[g])
        perm[k] = foggy_shift(dos, y)
    p_stand = (np.sum(perm >= obs) + 1) / (N_PERM + 1)

    # 3. centromere-matched null
    chrom, a, b = PEAK_SPAN
    dist = tn.cen_dist(chrom, (a + b) / 2)
    band_id = f"{chrom}:{a}"
    o = pd.read_csv(FOG_OUTLIERS, sep=" ", header=None).iloc[:, :3]
    smap = pd.read_csv(f"{GEMMA_DIR}/scaffold_snp_map.txt", sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    o["chrom"] = o[0].map(dict(zip(smap[0], smap[2]))); o = o.sort_values(["chrom", 2])
    o["pk"] = ((o.chrom != o.chrom.shift()) | (o[2].diff() > 40_000)).cumsum()
    fog_peaks = o.groupby("pk").agg(chrom=("chrom", "first"), start=(2, "min"), end=(2, "max"))
    C = tn.cen_null_dosages(list(env.index), {band_id: tn.cen_band(dist)}, fog_peaks)
    C = C[C.index.get_level_values("band") == band_id]
    C = tn.polarize(C)[by_year["SampleID"].tolist()]
    cf = C.mean(axis=1) / 2
    C = C[((cf - f_peak).abs() <= NULL_FREQ_TOL) & (cf > 0) & (cf < 1)]
    creg = [f"{c}:{int(p) // 2000}" for c, p, _ in C.index]
    C = C.groupby(creg, group_keys=False).apply(lambda x: x.sample(1, random_state=1))
    cnull = np.array([foggy_shift(row.values, years) for _, row in C.iterrows()])
    p_cen = (np.sum(cnull >= obs) + 1) / (len(cnull) + 1)
    lo, hi = tn.cen_band(dist)

    res = pd.DataFrame([dict(statistic=f"freq({FOGGY_YEAR}) - freq(other years)", observed=obs,
                             peak_haplotype_freq=f_peak, n_random_snps=len(null),
                             random_snp_null_p=p_null, random_snp_null_95pct=np.percentile(null, 95),
                             within_stand_perm_p=p_stand, within_stand_perm_95pct=np.percentile(perm, 95),
                             n_individuals=len(by_year), n_stands=len(groups),
                             peak_cen_dist_Mb=dist / 1e6, cen_band_Mb=f"{lo/1e6:.2f}-{hi/1e6:.2f}",
                             n_random_snps_cen=len(cnull), cen_matched_null_p=p_cen,
                             cen_matched_null_95pct=np.percentile(cnull, 95))])
    res.to_csv(TEMPORAL_OUT, sep="\t", index=False)
    print(f"Temporal checks: favorable-haplotype frequency {FOGGY_YEAR} minus other years = {obs:+.4f} "
          f"(overall freq {f_peak:.3f}, n = {len(by_year)})")
    print(f"  1. vs {len(null)} fog-polarized, frequency-matched random SNPs: one-sided p = {p_null:.3f} "
          f"(null 95th percentile {np.percentile(null, 95):+.4f})")
    print(f"  2. years permuted within {len(groups)} stands ({N_PERM} permutations): one-sided p = {p_stand:.4f} "
          f"(95th percentile {np.percentile(perm, 95):+.4f})")
    print(f"  3. vs {len(cnull)} random SNPs {lo/1e6:.2f}-{hi/1e6:.2f} Mb from their centromere (peak: "
          f"{dist/1e6:.2f} Mb): one-sided p = {p_cen:.3f} (null 95th percentile {np.percentile(cnull, 95):+.4f})")
    print(f"  -> {TEMPORAL_OUT}")


# ── Figure ───────────────────────────────────────────────────────────────────
synteny = CFG.get("synteny", False)
height_ratios = ([0.65] if synteny else []) + [0.6, 0.5, 0.54, 0.42] + ([1.2, 0.8] if lines else [])
o = 1 if synteny else 0       # index offset of the panels below the synteny panel
bottom_row = CFG.get("bottom_row", False)
# bottom row: a spacer (F keeps its x-axis title) then panels G-I
grid_ratios = height_ratios + ([0.3, 1.0] if bottom_row else [])
fig = plt.figure(figsize=(14, 3.1 * sum(grid_ratios)))
gs = fig.add_gridspec(len(grid_ratios), 1, height_ratios=grid_ratios, hspace=0.04)
axes = [fig.add_subplot(gs[0])]
axes += [fig.add_subplot(gs[i], sharex=axes[0]) for i in range(1, len(height_ratios))]
for ax in axes[:-1]:
    ax.tick_params(labelbottom=False)
bottom_axes = []
if bottom_row:
    sub = gs[len(height_ratios) + 1].subgridspec(1, 3, wspace=0.75, width_ratios=BOTTOM_WIDTHS)
    bottom_axes = [fig.add_subplot(sub[0, k]) for k in range(3)]

for ax in axes:
    for lo, hi in CFG["shade"]:
        ax.axvspan(lo, hi, color="#808080", alpha=SHADE_ALPHA, zorder=0)
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.spines[["top", "right"]].set_visible(False)

# Synteny panel: black bars = ancestral segments, hatched = reversed run
if synteny:
    ax = axes[0]
    pos = []
    for sp in SYNTENY_SPECIES:
        h = pd.read_csv(f"{SYNTENY_DIR}/{sp}_vs_Deinandra.synHits.txt.gz", sep="\t", low_memory=False,
                        usecols=["chr2", "start2", "end2", "blkID"])
        h = h[(h.chr2 == CHROM) & h.blkID.notna()]
        pos += [h.start2.values / 1e6, h.end2.values / 1e6]
    pos = np.sort(np.concatenate(pos))
    breaks = np.where(np.diff(pos) > SYNTENY_GAP)[0]
    segs = list(zip(np.r_[pos[0], pos[breaks + 1]], np.r_[pos[breaks], pos[-1]]))
    print("Ancestral segments (Mb): " + ", ".join(f"{a:.2f}-{b:.2f}" for a, b in segs))
    ax.broken_barh([(a, b - a) for a, b in segs], (0.3, 0.4), facecolor="black", edgecolor="black")
    ax.broken_barh([(a, b - a) for a, b in SYNTENY_REVERSED], (0.3, 0.4), facecolor="white",
                   edgecolor="black", hatch="////", linewidth=1.2, zorder=3)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_ylabel("Ancestral\nsegments", fontsize=LABEL_SIZE)

# Fog Manhattan (peak-filtered points from annotated_manhattan.R).
# y-axis tight to the top SNP; gene labels sit off to the side of their peak.
ax = axes[o + 0]
ax.scatter(man_df["Position"] / 1e6, man_df["logp"], s=14, color="black", alpha=0.65,
           edgecolor="none", rasterized=True)
ax.axhline(thresh_logp, color="grey", linestyle="--", linewidth=1.5)
for _, row in man_labels.iterrows():
    ax.annotate(CFG.get("label_text", {}).get(row["Label"], row["Label"]), xy=(row["Position"] / 1e6, row["logp"]),
                xytext=(3, 6), textcoords="offset points", ha="left", va="center",
                fontsize=LABEL_SIZE, fontstyle="italic",
                arrowprops=dict(arrowstyle="-", color="grey", linewidth=0.8))
ax.set_ylim(0, man_df["logp"].max() * (1 + Y_PAD))
ax.set_ylabel("-log10\np-val", fontsize=LABEL_SIZE)

# Panels 2-4 share one style: a black line filled underneath with black.
# Gene density and LINE content are per 500 kb window, so they're drawn as
# steps (each window flat across its span); dashed grey line = genome-wide
# value, drawn on top so it shows over the fill.
def window_track(ax, win_df, col):
    edges = np.append(win_df["start"].values, win_df["end"].values[-1]) / 1e6
    ax.stairs(win_df[col].values, edges, fill=True, color="black")
    ax.set_ylim(0, win_df[col].max() * (1 + Y_PAD))


ax = axes[o + 1]
window_track(ax, gd_win_df, "density")
ax.axhline(refs["gene_density"], color="grey", linestyle="--", linewidth=1.5, zorder=3)
ax.set_ylabel("Genes\nper Mb", fontsize=LABEL_SIZE)

ax = axes[o + 2]
window_track(ax, rep_win_df, "LINE")
ax.axhline(refs["repeat_LINE"], color="grey", linestyle="--", linewidth=1.5, zorder=3)
ax.set_ylabel("% LINE\ncontent", fontsize=LABEL_SIZE)

# Panel 4: % variance explained by PC1 in each window
ax = axes[o + 3]
if pca_stat is not None:
    ve_x, ve = pca_stat.index.values / 1e6, pca_stat["pc_1_ve"].values
    ax.fill_between(ve_x, ve, 0, color="black", linewidth=0)
    ax.plot(ve_x, ve, color="black", linewidth=1.2)
    ax.set_ylim(0, np.nanmax(ve) * (1 + Y_PAD))
ax.set_ylabel("PC1\n% var.", fontsize=LABEL_SIZE)

# Panels 5-6: winPCA, one line per sample, drawn in descending fog so the
# lowest-fog samples are drawn last (on top); samples with no fog value first
all_fog = -pd.to_numeric(pca_meta[WINPCA_COLOR_BY], errors="coerce")
norm = matplotlib.colors.Normalize(vmin=fog_val.min() if len(fog_val) else all_fog.min(),
                                   vmax=fog_val.max() if len(fog_val) else all_fog.max())
keep = fog_val[fog_val.index.isin(drawn_samples)]
sample_order = list(keep.index[keep.isna()]) + list(keep.dropna().sort_values(ascending=False).index)


def winpca_panel(ax, df, ylabel, years=None):
    x = df.index.values / 1e6
    drawn = [s for s in sample_order if s in df.columns and (years is None or year[s] in years)]
    print(f"{ylabel} panel: {len(drawn)} samples drawn")
    for sample in drawn:
        val = fog_val[sample]
        color = "lightgrey" if pd.isna(val) else WINPCA_CMAP(norm(val))
        ax.plot(x, df[sample].values, color=color, linewidth=LINE_WIDTH, alpha=LINE_ALPHA, rasterized=True)
    ax.margins(y=Y_PAD)
    ax.set_ylabel(ylabel, fontsize=LABEL_SIZE)


if lines and pc1 is not None:
    winpca_panel(axes[o + 4], pc1, "PC1", years=PC1_YEARS)
elif lines:
    axes[o + 4].set_ylabel("PC1", fontsize=LABEL_SIZE)
if lines and hetp is not None:
    winpca_panel(axes[o + 5], hetp, "Proportion\nhet. sites")
elif lines:
    axes[o + 5].set_ylabel("Proportion\nhet. sites", fontsize=LABEL_SIZE)

axes[-1].set_xlabel(f"Position on chromosome {CHROM} (Mb)", fontsize=LABEL_SIZE)
axes[-1].set_xlim(0, CHROM_LEN / 1e6)

# Panels G-I (chr 6 only)
if bottom_row:
    by_year = pd.read_csv(BOTTOM_BY_YEAR, sep="\t")
    years = sorted(by_year["Year"].unique())
    temporal_checks(by_year)

    # G: unit long-term fog vs favorable-allele frequency, one line per peak
    ax = bottom_axes[0]
    un = pd.read_csv(BOTTOM_UNITS, sep="\t")
    xs = np.array([un.longterm_fog.min() - 0.05, un.longterm_fog.max() + 0.05])
    for pk, u in un.groupby("peak"):
        is6 = u.chr6.iloc[0]
        ax.plot(xs, u.intercept_longterm.iloc[0] + u.slope_longterm.iloc[0] * xs,
                **(dict(color="black", lw=2.5, ls="--", zorder=4) if is6 else
                   dict(color=TRANQ_COLOR, lw=1.2, alpha=OTHER_ALPHA, zorder=2)))
        if is6:
            ax.scatter(u.longterm_fog, u.freq, s=u.n * 1.5 + 20, marker="D", color=TRANQ_COLOR,
                       edgecolor="black", linewidth=1.3, zorder=5)
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
    ax.set_xlabel("Long-term fog\n(days/month)", fontsize=LABEL_SIZE)
    ax.set_ylabel("Favorable\nallele freq.", fontsize=LABEL_SIZE)

    # H: germination-year fog by germination year (Tranquillon Mtn plants)
    ax = bottom_axes[1]
    vals = [by_year.loc[by_year.Year == y, "germ_env"].dropna().values for y in years]
    ax.boxplot(vals, positions=range(len(years)), widths=0.6, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=TRANQ_COLOR, edgecolor="black", linewidth=1.5),
               medianprops=dict(color="black", linewidth=2.5),
               whiskerprops=dict(linewidth=1.5), capprops=dict(linewidth=0))
    ax.set_xticks(range(len(years)), [str(y) for y in years], rotation=45, ha="right", rotation_mode="anchor")
    ax.set_xlabel("Germination year", fontsize=LABEL_SIZE)
    ax.set_ylabel("Germination year\nfog (days/month)", fontsize=LABEL_SIZE)

    # I: germination-year fog vs unit-adjusted frequency per cohort, one line per peak
    ax = bottom_axes[2]
    py = pd.read_csv(BOTTOM_PEAK_YR, sep="\t")
    for pk, c in py.groupby("peak"):
        c = c.sort_values("year")
        is6 = c.chr6.iloc[0]
        jx = np.array([c.june_fog.min() - 0.3, c.june_fog.max() + 0.3])
        ax.plot(jx, c.mean_freq.iloc[0] + c.slope_june.iloc[0] * (jx - c.june_mean.iloc[0]),
                **(dict(color="black", lw=2.5, ls="--", zorder=4) if is6 else
                   dict(color=TRANQ_COLOR, lw=1.2, alpha=OTHER_ALPHA, zorder=2)))
        if is6:
            ax.errorbar(c.june_fog, c.unit_adjusted_freq, yerr=c.ci95_half, fmt="D", color=TRANQ_COLOR,
                        mec="black", mew=1.3, ecolor="black", markersize=13, elinewidth=2, capsize=8,
                        capthick=2, zorder=5)
            for y, xv, yv in zip(c.year, c.june_fog, c.unit_adjusted_freq):
                dx, dy, ha = YEAR_LABEL_OFFSET.get(y, (12, 4, "left"))
                ax.annotate(str(y), (xv, yv), xytext=(dx, dy), textcoords="offset points", ha=ha,
                            va="center", fontsize=TICK_SIZE - 2, zorder=6)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Germination year\nfog (days/month)", fontsize=LABEL_SIZE)
    ax.set_ylabel("Favorable\nallele freq.", fontsize=LABEL_SIZE)

    for ax in bottom_axes:
        ax.tick_params(axis="both", labelsize=TICK_SIZE)
        ax.spines[["top", "right"]].set_visible(False)

fig.align_ylabels(axes)
fig.subplots_adjust(left=0.1, right=0.9, top=0.99, bottom=0.04)

# Warn if neighbouring y-axis titles collide (a title may extend past its own
# panel, as long as it clears the next panel's title)
fig.canvas.draw()
label_boxes = [ax.yaxis.label.get_window_extent() for ax in axes]
for upper, lower, a_up, a_lo in zip(label_boxes, label_boxes[1:], axes, axes[1:]):
    gap = upper.y0 - lower.y1
    if gap < 15:
        print(f"WARNING: y-axis titles {a_up.get_ylabel()!r} and {a_lo.get_ylabel()!r} "
              f"collide ({gap:.0f} px apart)")

# Panel letters left-justified in one column, whose widest letter ends just
# left of the leftmost y-axis label
inv_fig = fig.transFigure.inverted()
label_x0 = min(inv_fig.transform(ax.yaxis.label.get_window_extent())[0, 0] for ax in axes)
letters = [fig.text(0, ax.get_position().y1, letter, fontsize=PANEL_LETTER_SIZE, ha="left", va="top")
           for ax, letter in zip(axes, PANEL_LETTERS)]
fig.canvas.draw()
widest = max(inv_fig.transform(t.get_window_extent())[1, 0] - inv_fig.transform(t.get_window_extent())[0, 0]
             for t in letters)
for t in letters:
    t.set_x(label_x0 - 0.003 - widest)
# G-I: the first sits in the C-F letter column; the others just left of
# their own y-axis title (and tick labels), level with the top of the panel
for k, (ax, letter) in enumerate(zip(bottom_axes, PANEL_LETTERS[len(axes):])):
    x0 = min(inv_fig.transform(ax.yaxis.label.get_window_extent())[0, 0],
             inv_fig.transform(ax.get_tightbbox(fig.canvas.get_renderer()))[0, 0])
    x = label_x0 - 0.003 - widest if k == 0 else x0 - 0.003 - widest
    fig.text(x, ax.get_position().y1 + 0.004, letter, fontsize=PANEL_LETTER_SIZE, ha="left", va="top")

# One colorbar spanning the winPCA panel(s), placed after layout is final
if lines:
    top_pos, bottom_pos = axes[o + 4].get_position(), axes[-1].get_position()
    cax = fig.add_axes([top_pos.x1 + 0.01, bottom_pos.y0 + 0.02 * bottom_pos.height,
                        0.012, top_pos.y1 - bottom_pos.y0 - 0.04 * bottom_pos.height])
    cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=WINPCA_CMAP), cax=cax, alpha=CBAR_ALPHA)
    cbar.set_label("Fog (days per month)", fontsize=LABEL_SIZE)
    cbar.ax.tick_params(labelsize=TICK_SIZE)

# bbox_inches="tight" trims the outer whitespace down to pad_inches
plt.savefig(CFG["out"], dpi=200, bbox_inches="tight", pad_inches=0.1,
            bbox_extra_artists=[ax.xaxis.label for ax in bottom_axes])
plt.close()
print(f"Saved: {CFG['out']}")
