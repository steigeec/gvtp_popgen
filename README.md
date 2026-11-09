# gvtp_popgen
Structure and selection in Deinandra increscens villosa

Scripts are grouped below by the sections of the paper's Methods. Each script
sets `BASE_DIR` (the project data directory) at the top; inputs and outputs are
read from and written to paths under it, not to this repository. Each script's
header lists the conda environment or R library it needs and how to install it.

## Genetic sampling, sequencing, and pre-processing

Range map (`phylogeny_range/`)

- `phylogeny-range.sh` (runs `phylogeny_range.R`): occurrence records and range map figure

Reads and mapping (`preprocessing/`)

1. `trim_reads.sh`: Trimmomatic adapter and quality trimming of each sample
2. `merge_fastq.sh`: combine a sample's trimmed reads across sequencing runs; `trim_summary.sh`: reads retained per sample
3. `parallelized_mapping.sh` (runs `GVTP_mapping.sh`): bwa mem mapping of paired and orphan reads
4. `parallelized_merging.sh` (runs `merge_bams.sh`): merge BAMs per sample and mark duplicates with Picard
5. `run_coverage.sh` (runs `coverage.sh`): mean depth per sample

Organellar scaffolds (`preprocessing/`)

6. `organelle_nucmer.sh`: NUCmer alignment of all scaffolds to Asteraceae plastomes and mitogenomes
7. `organelle_scaffolds.py`: classify plastid and mitochondrial scaffolds (`chloroplast.txt`, `mito.txt`) from the alignments and read depth

Genotype likelihoods and site filtering (`variant_calling_filtering/`)

8. `launch_variant_call_noMAF.sh` (runs `variant_call_noMAF.sh`): ANGSD genotype likelihoods per genome chunk
9. `preParalog_merge_and_filter_snps.sh`: merge chunks
10. `launch_ngsParalog.sh` (runs `ngsParalog_chunk.sh`), then `site_filter.R`: paralog filter
11. `filterParalogMAF_prepNGSld.sh`: autosomal and MAF-filtered site sets and beagle subsets; runs `ngsLD_launch.sh` (`ngsLD.sh`)
12. `plotLDdecay.sh`: LD decay curves, used to choose the pruning distance
13. `launch_PruneGraph.sh` (runs `PruneGraph_array.sh`), then `make_LD_sites.sh`: LD-pruned site set
14. `launch_heatmapLD.sh` and `launch_heatmap_plot.sh` (runs `heatmap_processed.sh`): LD heatmaps
15. `imputation.sh`: BEAGLE genotype imputation (DR2 > 0.8); `beagle_filter_individuals.sh` subsets beagle files by individual

Invariant sites (`selection/TajD/`)

16. `launch_angsd_w_invar.sh` (runs `angsd_w_invar.sh`): allele frequencies keeping invariant sites
17. `VerifyPolyV11.py`, then `merge_filter.sh`: combine monomorphic and paralog-filtered variable sites into the all-sites set

## Range-wide structure and genetic diversity

Structure (`structure/`)

1. `run_PCAngsd.sh` (runs `PCAngsd_filter.sh`, `PCAngsd_nofilter.sh`), then `plot_PCA.R`: PCA
2. `launch_ohana_year1_strauss.sh` (runs `ohana.sh`), then `OHANA_plot_year1.R` (range-wide, 2020) and `OHANA_plot_Strauss.R` (Tranquillon Mtn): admixture and NMI-based K selection
- `env/stand_map.py`: map of sampling stands within Tranquillon Mtn

Population diversity, Ne and differentiation

3. `pop_SAF.sh`: per-population SAFs over the all-sites set
4. `diversity/heterozygosity/population/pop_het.sh`, then `pop_het_calc.R`: nucleotide diversity
5. `Ne/run_watterson_theta.sh`, then `Ne/summarize_watterson_theta.sh`: Watterson's theta and Ne
6. `structure/calculate_Fst_pairwisepop.sh` (runs `pop_pairwise_array.sh`, `pop_collect_results.sh`), then `structure/plot_Fst_heatmap.R`: pairwise FST
7. `wind/fst_downsample_launch.sh covmatch` (runs `fst_downsample_saf.sh`, `fst_downsample_pair.sh`): FST with coverage-matched downsampling
8. `structure/mantel_IBD.R`: isolation by distance among populations

Individual diversity and differentiation

9. `diversity/heterozygosity/individual/launch_He.sh` (runs `ind_het.sh`), then `het_plot.R`: coverage-corrected individual diversity and tests among localities
10. `structure/calculate_Fst_individuals.sh` (runs `pairwise_array.sh`, `collect_results.sh`): pairwise FST between individuals
11. `structure/mantel_IBD_individuals.R`: isolation by distance within Tranquillon Mtn

Inbreeding and runs of homozygosity (`diversity/`)

12. `inbreeding/doGlf_launch.sh` (runs `createGlf_ngsF.sh`), then `inbreeding/launch_ngsF.sh` (runs `ngsF.sh` and `inbreeding_plot.R`): ngsF inbreeding coefficients
13. `ROHan/prep_samples.sh`, then `ROHan/estimate_rohmu.sh` (runs `rohmu_phase1.sh`, `rohmu_phase2.sh`): rohmu calibration
14. `ROHan/launch_rohan.sh` (runs `rohan_phase1.sh`, `rohan_phase2.sh`; `resubmit_phase2.sh` resubmits unfinished samples)
15. `ROHan/collect_rohan.sh`, then `ROHan/plot_rohan.sh` (runs `roh_plot.R`): ROH burden and tests among years

Effective migration surfaces (`migration/`)

16. `create_grid.sh` (runs `create_grid.R`): spatial grid shared by FEEMS and FRAME
17. `FEEMS.sh`, `FEEMS_year2.sh`, `FEEMS_year3.sh` (run `FEEMS.py`) and `FEEMS_fit.sh` (runs `FEEMS_fit.py`): FEEMS per sampling year, with lambda chosen by cross-validation

## Wind connectivity analysis

Tests whether wind connectivity explains genetic differentiation and diversity
among populations (HRRR wind) and among individuals on Tranquillon Mountain
(WindNinja-downscaled wind), using windscape (Kling & Ackerly 2020, 2021), and
whether gene flow estimated with FRAME is biased downwind. Wind is a pre-turbine
climatology: June-October 2015-2022, 6-hourly (the Strauss Wind Farm began
operating in 2023).

Wind data (`wind/`)

1. `download_station_wind.sh`, `build_station_obs_table.py`: station observations (PTGC1, VGN, TQS)
2. `download_hrrr_wind.py`, then `rotate_hrrr_wind.py`: HRRR 10 m wind, rotated to true north
3. `validate_hrrr_vs_stations.py`: HRRR vs station observations
4. `build_windninja_dem.sh`, `build_coastline.sh`, `build_landfire_lcp.sh`: terrain, coastline and LANDFIRE 2016 landscape files (`build_landfire_lcp.sh` needs `LFPS_EMAIL` set)
5. `download_hrrr_full_subset.py`: HRRR grib2 files used to initialize WindNinja
6. `sbatch run_windninja_landfire_batch.sh region` and `sbatch run_windninja_landfire_batch.sh tranquillon`; `sbatch --array=0-4 run_windninja_tranquillon_sensitivity.sh`
7. `validate_windninja_stations.py`: WindNinja vs station observations

Connectivity and tests (`wind/`)

8. `windscape_build_series.R`: wind time series rasters for windscape
9. `windscape_connectivity.R`: wind roses and pairwise travel times
10. `windscape_tests.R`: isolation and diversity tests (populations and Tranquillon Mtn)
11. `windscape_population_random_walk.R`: random-walk connectivity at the population level
12. `Rscript windscape_tranquillon_sensitivity.R base mesh100 landfire diurnal stability full`: WindNinja configuration sensitivity
13. `windscape_composite_figure.R`: connectivity figure

Directional gene flow (`migration/`; uses the grid from `create_grid.sh` above)

14. `FRAME_fit.sh`, `FRAME_lambda_sensitivity.sh` (run `FRAME_fit.py`): FRAME pooled across years
15. `frame_wind_asymmetry.R`: gene flow asymmetry vs wind flow asymmetry, with a rotated-wind null (run after step 8)

## Environmental variable selection

1. `env/extract_env.R`: environmental values at sample coordinates (`env_var.txt`)
2. `pheno/extract_pheno_env_data.R`: fog (FLCC and Kling), PRISM precipitation and BCM climatic water deficit at each coordinate, as germination-year and 20-year rolling means
3. `selection/select_env_var.R`: candidate variables, resolution and missingness filters, and pruning of collinear variables
4. `selection/gemma/PC_env_creation.R`: PCA of the retained variables; GEMMA phenotype and covariate files
- `env/compare_Werner_Kling_fog.R`: MODIS (Werner) vs Kling fog over their shared years

## Genomewide signatures of local adaptation

GEMMA (`selection/gemma/`)

1. `prep_GEMMA.sh`: PLINK files for the >4x individuals and `scaffold_snp_map.txt` (runs `PC_env_creation.R`)
2. `run_GEMMA.sh`: univariate LMM per environmental variable (`check_qq_lambda.R` gives a quick lambda/QQ check)
3. `selection/Launch_find_outlier_annotate.sh`: per variable, outliers and peaks (`find_outliers.R`), SnpEff annotation (`selection/annotateVCF_Di2.sh`), LD-informed gene assignment (`curate_intergenic_by_ld.py`, `resolve_gene_ids.R`), orthologs (`build_gene_orthology.py`, `blast_chunks.sh`, `lookup_by_gene_id.py`) and gene summary (`compile_gene_summary.py`, `build_snp_gene_map.py`, `compute_gene_snp_stats.py`)
4. `run_40000_3_annotation.sh`: annotation and gene lists for the multi-SNP peaks
5. `build_gene_summary_all.py`: candidate genes across variables

Genome-wide LD scan used for intergenic SNPs (`selection/gemma/`)

- `find_intergenic_candidate_genes.py`, then `launch_ld_scan.sh` (runs `ld_scan_array.sh`, `ld_scan_scaffold.py`)

Permutation tests and gene-set overlap (`selection/gemma/`)

6. `make_permutations.R`, then `run_BIO18_5pct_observed.sbatch` and `run_BIO18_5pct_perm_array.sbatch`, then `summarize_BIO18_5pct.sbatch` (runs `summarize_permutation_null.R`): phenotype-permutation null for lambda and Bonferroni hits
7. `gene_overlap_across_variables.py`: hypergeometric tests of shared candidate genes

Figures

- `selection/gemma/annotated_manhattan.R` (with `select_peak_labels.R`): Manhattan plots
- `selection/compare_env_sites.R`: environmental conditions across localities for variables with outliers

## Gene ontology enrichment & tissue-level expression of candidate genes

GO annotation (`annotation/`)

1. `setup_eggnog_mapper.sh`, then `run_eggnog_mapper.sh` (runs `convert_eggnog_to_geneid_go.py`): eggNOG-mapper GO terms
2. `harvest_blast_go.py`: NCBI Gene GO terms from each gene's best BLASTp hit
3. `fix_gene_mRNA_id.py`, `translate_cds_to_pep.py`, `prep_tair_proteome.sh`, `prep_ensembl_plant_proteome.sh`, `prep_ncbi_datasets_proteome.sh`, `prep_ensembl_plant_gff.sh`: proteomes for orthology
4. `setup_orthofinder.sh`, then `run_orthofinder.sh`: orthogroups
5. `dump_tair_go.R`, then `build_tair_ortholog_go.py`: Arabidopsis GO transfer; `dump_biomart_go.sh`, then `build_ortholog_go.py`: sunflower GO transfer
6. `merge_go_sources.sh`: union of the five GO sources
7. `build_go_ancestor_map.py`: true-path propagation from go-basic.obo; `build_nonplant_go_exclusion_list.py`, `build_taxon_constraint_go_exclusion_list.py`, `build_housekeeping_go_exclusion_list.py`: excluded terms

Enrichment (`annotation/`, `selection/gemma/`)

8. `annotation/build_snp_gene_map.sh` (runs `parse_snp_gene_ann.py`), then `annotation/build_gene_gwas_scores.py`: ACAT gene-level scores
9. `GSEA_GO.R`: fgsea enrichment per variable
10. `cluster_gsea_go_results.R` and `gsea_go_stricter_nominal.R`: nominal terms, reduced with rrvgo; `aggregate_gsea_top_nominal.py`: table across variables
11. `gsea_go_convergence.R`: convergence across variables; `gsea_go_lollipop_3panel.R`: figure

Tissue expression (`RNAseq/`)

12. `run_quant.sh`: salmon quantification of the eight RNA-seq libraries
13. `build_gene_map.py`, then `aggregate_expression.py`: TPM of candidate genes
14. `supp_fig_tissue_heatmap.R`: expression heatmap

## Genomic context and dynamics of the DEINVI_00006 fog-associated haplotype

Windowed PCA and haplotype groups (`inversion/`)

1. `winPCA_maf5.sh`, `winPCA_maf5_DEINVI_00006.sh`, `winPCA_maf5_DEINVI_00006_Strauss.sh`: windowed PCA (genotype likelihoods, MAF >= 0.05); `winPCA_GT.sh`: genotype-based heterozygosity track; `plot_with_retry.sh` (runs `small_chromplot.py`): per-scaffold plots
2. `pc_cluster_consistency.py`: consistency of the PC1/PC2 groups along the scaffold

Inversion tests (`inversion/`)

3. `breakpoint_samples.py`, then `breakpoint_reads.sh`, then `breakpoint_summary.py`: short-read breakpoint search
4. `longread_chr6.sh`, then `longread_inv1.py`: PacBio HiFi and Nanopore reads across the block edges
5. `ld_within_groups.py`: LD decay within homokaryotype groups

Centromeres (`inversion/`, `preprocessing/`)

6. `scaffold_windows.py`, `inversion_analysis.py`: gene density and repeat content in 500-kb windows; `preprocessing/centromere_scaffolds.py`: figure for the 12 largest scaffolds
7. `centromere_scan.py`: tandem satellite arrays; `call_centromeres.py`: centromere calls
8. `outlier_centromere_enrichment.py`: outlier peaks vs centromere proximity
9. `compute_deinvi6_mapq.py`: mapping-quality check along DEINVI_00006

Synteny (`annotation/`, `inversion/`)

10. `annotation/setup_genespace.sh`, `annotation/run_orthofinder_genespace.sh`, then `setup_genespace_inputs.sh` (runs `build_genespace_bed.py`), then `run_genespace.sh` (runs `run_genespace.R` with `genespace_patches.R`)
11. `synteny_junctions.py`: centromeres at junctions between ancestral chromosomes; `gns1_snf1_interval.py`: genes between GNS1 and SNF1
12. `plot_deinandra_scaffold6.R`, `chr6_synteny_panels.py`: synteny figures

Haplotype dynamics within Tranquillon Mtn (`selection/gemma/`)

13. `fog_hypothesis_tests.py` (uses `temporal_null.py`): spatial and temporal tests with permutation and random-SNP nulls
14. `fog_peaks_temporal.py`, `germ_year_polygenic.R`: the same tests for the other fog peaks
15. `fog_peak_figures.py`, `fog_deinvi6_figure.R`, `inversion/DEINVI_0006_stacked.py`: figures

Centromere drive

16. `selection/TajD/submit_TajD_all.sh` (runs `run_TajD_all.sh`, `run_TajD.sh`, `Tajimas_SNPmidpoint.py`): Tajima's D per population, then `inversion/tajd_centromere_matched.py`
17. `inversion/centromere_haplotype_fst.py`, `inversion/centromere_hwe_check.py`: centromeric haplotype frequencies and F
18. `selection/window_Fst/window_Fst.sh` (runs `plot_window_Fst.R`): windowed FST, summarized at centromeres by `centromere_haplotype_fst.py`
19. `inversion/cenh3_collect_cds.py`, then `inversion/run_cenh3_dnds.sh` (runs `cenh3_dnds.py`): CENH3 dN/dS

## Phenotypic responses to environmental moisture

1. `pheno/extract_pheno_env_data.R`: fog and summer precipitation for each specimen (germination year and 20-year rolling mean)
2. `pheno/ray_count_validation.R`: iNaturalist vs dissected ray-floret counts
3. `pheno/phenotype_data.R`: differences among localities, GLMMs of ray-floret proportion and height:width ratio, and trends over time
4. `env/fog_trend.R`: fog trend at each plant's location from the Kling record
