# Methods

## Fog variability metrics

Fog exposure was quantified from the MODIS Monthly FLCC (Fog and Low Cloud Cover) raster dataset (Werner et al. 2022), which reports the number of days per month with detected fog/low cloud cover at 1 km resolution across coastal California and southern Oregon, for the fog season (June–September) of each year from 2000–2022. For each sample coordinate, the four within-season monthly rasters (June, July, August, September) were extracted and averaged to obtain a single fog-season mean value per site per year, giving a 23-year annual time series per sample.

From this per-sample annual time series, three summary metrics of long-term fog variability were computed:

- **`fog_SD`**: the standard deviation of the 23 annual fog-season means (previously computed; used here as a baseline).
- **`fog_CV`**: the coefficient of variation, `fog_SD / mean(annual values) × 100`, which normalizes interannual variability by each site's own long-term mean fog level. This was added because `fog_SD` alone is confounded with mean fog level (sites with more fog mechanically have more room for large absolute swings), whereas `fog_CV` is comparable across sites with different baseline fog exposure — consistent with the coefficient-of-variation convention already used for precipitation seasonality (BIO15) elsewhere in this pipeline.
- **`fog_detrend`**: the standard deviation of the residuals after fitting a linear model of the annual fog-season mean against year (`lm(annual_value ~ year)`) at each site (sites with fewer than 3 non-missing years were left as `NA`). This isolates year-to-year variability ("noise") from any long-term directional trend in fog (e.g., a steady multi-decadal decline would otherwise inflate raw SD without representing unpredictability).

These metrics were added to `env/extract_env.R`. Before modifying the production `env_var.txt` file, the extraction pipeline was rerun end-to-end and the recomputed `fog_SD` was checked against the existing `fog_SD` column; all 978 samples matched exactly (maximum absolute difference = 0), confirming the raster-extraction and row-ordering logic were reproduced correctly. The output was first written to a separate file and diffed column-by-column and row-by-row against the original `env_var.txt` to confirm that row order, row count, and all pre-existing columns were unchanged before the two new columns (`fog_CV`, `fog_detrend`) were added and the file was overwritten (with a backup retained).

## Population-level environmental and inbreeding covariates

For each of the five populations with an available Tajima's D track (CDFW/Gaviota, Strauss/Tranquillon Mtn, Point Conception, VAFB Lion's Head, VAFB Arguello — Dangermond was excluded, as no Tajima's D track exists for it), population-level covariates were computed as the mean across exactly the set of samples used to build that population's SAF/Tajima's D track — i.e., the same bamlist passed to ANGSD in `run_TajD.sh` for that population (16, 277, 16, 15, and 14 individuals for CDFW, Strauss, Point Conception, Lion's Head, and Arguello respectively). This ensured the environmental and inbreeding covariates describe the same individuals whose genomic data produced the Tajima's D estimates, rather than a possibly different or larger sampling pool. (Strauss's Tajima's D track was built from the full 277-sample bamlist, not the smaller 16-sample high-coverage subset used for an unrelated SAF run elsewhere in the project — this was initially confused and later corrected by cross-checking each population's ANGSD `.arg` log against the bamlist actually passed to `run_TajD.sh`.)

- **Inbreeding coefficient**: mean per-individual F from `ngsF`, averaged over each population's TajD sample set.
- **Fog metrics**: mean of `fog` (raw fog days/month), `fog_CV`, and `fog_detrend` (see above), averaged over the same sample sets.

## Windowed Tajima's D computation

Genome-wide, windowed Tajima's D was computed per population using a genotype-likelihood-based pipeline (`submit_TajD_all.sh`, `run_TajD.sh`, `Tajimas_SNPmidpoint.py`), run separately for each of the five populations with sufficient sample size (CDFW, Strauss, Point Conception, VAFB Lion's Head, VAFB Arguello; Dangermond was not run through this pipeline and has no Tajima's D track).

For each population, site allele frequency likelihoods (SAF) were estimated in ANGSD (`-doSaf 1`, `-GL 2`, the GATK genotype-likelihood model) from that population's BAM list, using the reference genome as both reference and ancestral sequence (no independent outgroup was available). Sites were restricted to a predefined genome-wide site list (`all_sites`, common across all populations) and a fixed set of scaffolds (`all_chrs.txt`, 21,073 scaffolds), with standard quality filters (`-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20`), a per-individual depth filter (`-setMinDepthInd 1 -setMaxDepthInd 25`), and a population-specific minimum-individual-coverage threshold (`-minInd`, scaled to each population's sample size: e.g. 50 of 277 for Strauss, 10–11 of 14–16 for the smaller populations). Because no ancestral-state information was available, the folded site frequency spectrum was estimated from the SAF file with `realSFS -fold 1`, and this per-population SFS was used as a prior to compute per-site posterior expectations of Watterson's θ and pairwise nucleotide diversity (π) via `realSFS saf2theta` followed by `thetaStat print`, which writes one log-scaled θ estimate per analyzed site (not only sites called as variant).

Non-overlapping windows of fixed SNP count (rather than fixed physical size) were then built from this per-site θ output by a custom script (`Tajimas_SNPmidpoint.py`, window size = 100 sites per window, restarting at each new scaffold). Window boundaries were placed at the genomic midpoint between the last site of one window and the first site of the next, so that adjacent windows partition the scaffold without gaps or overlap. Within each 100-site window, the per-site Watterson's θ and π estimates (stored as natural logs by `thetaStat`) were exponentiated and summed to give window-total θ_W and θ_π, and Tajima's D was calculated from these window totals using the standard Tajima (1989) estimator and its associated variance correction (constants a1, a2, b1, b2, c1, c2, e1, e2, computed from the population's sample size as 2N chromosomes):

D = (θ_π − θ_W) / √(e1·S + e2·S·(S−1))

where S = θ_W · a1 is the number of segregating sites implied by θ_W. This produced one Tajima's D (and θ_W, π) value per 100-site window per population, output as `<prefix>.Diana.out`, which is the genome-wide Tajima's D track used throughout the rest of this analysis (loaded into `tajd_by_pop` in `plot_TajD_extreme_combined.R`).

Note that because window size is defined in SNPs/sites rather than base pairs, windows are not of uniform physical length — they are narrower in high-diversity/high-density regions and wider in low-diversity/sparser regions of the genome.

## Tajima's D outlier windows and the putative inversion

For each environmental variable with a GEA (genotype–environment association, GEMMA) outlier scan, outlier SNPs were mapped to their nearest genomic window (nearest-match on window center, per scaffold) within each population's Tajima's D track. A (population, variable) combination was retained only if at least `MIN_WINDOWS = 5` outlier SNPs mapped to windows in that population, avoiding conclusions drawn from a handful of windows.

Two non-mutually-exclusive outlier categories were defined per population, pooled across all retained environmental variables (with windows flagged by more than one variable's scan counted only once — see below):

- **"Outliers"**: the full set of matched outlier windows, regardless of genomic location.
- **"Putative inversion outliers"**: the subset of matched outlier windows whose full window extent (start and end) falls entirely within the putative inversion region (scaffold `DEINVI_00006`, 43–50 Mb). Windows that only partially overlap the inversion boundary were excluded from this category (but remain in "Outliers"). By construction, "Putative inversion outliers" is a subset of "Outliers," not a disjoint category.

## Delta Tajima's D (population-level summary)

For each population, a genome-wide background Tajima's D was defined as the median Tajima's D across all windows in that population's track. Delta Tajima's D for a given category was then:

Δ*D* = median(Tajima's D of category windows) − genome-wide background median

This was computed separately for the "Outliers" and "Putative inversion outliers" categories per population (each requiring ≥ `MIN_WINDOWS` windows to be plotted), and separately for two variable-pooling schemes: pooled across all retained environmental variables, and restricted to outlier windows from the "fog" GEA scan only (`_fog`-suffixed outputs).

These Δ*D* values were plotted against each population's inbreeding coefficient and fog metrics, with one point per population per category: fill color indicating population (matching the palette used elsewhere in this project) and point shape indicating category (circle = "Outliers," triangle = "Putative inversion outliers," both drawn with a black outline and shape-fill mapped to population). The population-level background category (Δ*D* = 0 by definition) is represented only by the dashed horizontal reference line rather than as a plotted point.

## Mixed-effects test of the Delta Tajima's D–environment relationship

Because the population-level summary above reduces the data to five points per category, it cannot support a well-powered significance test and is sensitive to any single population. A linear mixed-effects model was therefore fit on the underlying window-level data rather than the population medians.

Unlike the population-level Δ*D* plots above (which pool outlier windows across all retained environmental variables), the mixed model was restricted to outlier windows from the "fog" GEA scan only, matching the `_fog`-suffixed descriptive plots. This restriction was necessary because the model tests whether Δ*D* in outlier windows scales with fog variability (or inbreeding): the windows entering that test need to already be fog-associated outliers, not GEA outliers for unrelated variables (e.g. precipitation seasonality, soil chemistry), which would dilute or confound a fog-specific signal. For each fog-associated outlier window contributing to either category, a window-level Δ*D* was computed as that window's Tajima's D minus its population's genome-wide background median. This produced 12,699 window-level observations across the two categories.

The model was:

```
delta_TajD ~ env_x * category + (1 + category | population)
```

fit by REML with `lme4`/`lmerTest`, where `env_x` is the population-level covariate (inbreeding coefficient or `fog_detrend`), `category` is "Outliers" vs. "Putative inversion outliers," and `population` is a random effect with both a random intercept and a random slope for `category`. The random slope is necessary because `category` varies *within* each population (every population contributes windows to both categories): a random-intercept-only model would let any population-specific idiosyncrasy in the "Outliers" vs. "Putative inversion outliers" contrast leak into the window-level residual term instead of the population-level random effect, artificially shrinking the standard error on the `env_x:category` interaction — the term that tests whether the inversion category responds more steeply to the environmental covariate than the rest of the outlier set. Because that comparison is fundamentally informed by only 5 independent populations, the model needs to reflect that in its uncertainty, not the ~12,700 window-level observations.

Fixed-effect contrasts of interest — the slope of Δ*D* on `env_x` within "Outliers," within "Putative inversion outliers," and the difference between them (the interaction) — were estimated with `lmerTest::contest1D()`, which uses Satterthwaite degrees of freedom (effectively ≈ 3, reflecting the 5-population design) rather than a naive Wald interval based on the full window-level sample size. Models were fit separately for the inbreeding coefficient and `fog_detrend`; the model for the inbreeding coefficient returned a singular fit (the population-level random-slope variance could not be reliably estimated) and its results were treated as unreliable rather than interpreted. The `fog_detrend` model converged without singularity.

Full model summaries and the three fixed-effect contrasts (with Satterthwaite-based 95% CIs and p-values) for both covariates were written to a text summary alongside the figures.

## Final visualization

The final figures plot, for each population and category, the median window-level Δ*D* together with a bootstrap 95% confidence interval (2,000 resamples of the underlying window-level Δ*D* values within that population × category, seed = 42), against the population's inbreeding coefficient or `fog_detrend` value. No fitted regression line or confidence band is overlaid on these plots: with only 5 populations, a continuous regression band implies a level of precision the data do not support, and gives a visually misleading impression of a smooth trend. The bootstrap whiskers instead give an honest, per-point sense of how well each population's summary value is determined by its underlying window count, while the formal slope/interaction test remains available in the mixed-model summary output.
