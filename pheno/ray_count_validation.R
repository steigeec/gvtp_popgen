## Paired Wilcoxon signed-rank test of ray-floret counts from dissection vs
## iNaturalist photos, for specimens counted both ways. A non-significant
## result justifies phenotype_data.R treating the two counts as interchangeable.
##
## Input: BASE_DIR/pheno/gtinatandvars.csv
## Usage: Rscript ray_count_validation.R

suppressMessages({
  library(dplyr)
  library(tidyr)
})

BASE_DIR <- "/private/groups/meyerlab/esteiger"
pheno_file <- file.path(BASE_DIR, "pheno/gtinatandvars.csv")
SIG_ALPHA <- 0.05

pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
pheno$recordNumber_CollectorNumber <- trimws(as.character(pheno$recordNumber_CollectorNumber))
pheno$year <- suppressWarnings(as.numeric(pheno$year))
pheno$Count_ray_flourets <- suppressWarnings(as.numeric(pheno$Count_ray_flourets))
pheno$iNat_ray_count     <- suppressWarnings(as.numeric(pheno$iNat_ray_count))

paired <- pheno %>%
  transmute(Count_ray_flourets, iNat_ray_count) %>%
  drop_na(Count_ray_flourets, iNat_ray_count)
cat(sprintf("Ray-count validation: %d specimens with both dissection and iNat counts\n", nrow(paired)))
if (nrow(paired) >= 3) {
  diffs     <- paired$Count_ray_flourets - paired$iNat_ray_count
  shapiro_p <- shapiro.test(diffs)$p.value
  wt <- wilcox.test(paired$Count_ray_flourets, paired$iNat_ray_count, paired = TRUE)
  cat(sprintf("  Shapiro-Wilk on differences: p = %.2g (%s)\n", shapiro_p,
              if (shapiro_p < 0.05) "non-normal, paired Wilcoxon used" else "normal"))
  cat(sprintf("  Paired Wilcoxon signed-rank: V = %.0f, p = %.3f (mean dissection = %.2f, mean iNat = %.2f)\n",
              wt$statistic, wt$p.value, mean(paired$Count_ray_flourets), mean(paired$iNat_ray_count)))
  if (wt$p.value < SIG_ALPHA) {
    warning("Dissection and iNat ray counts differ significantly (p = ", round(wt$p.value, 3),
            ") — fusing them in phenotype_data.R is NOT well-justified by this test.")
  } else {
    cat("  Not significant — dissection and iNat ray counts are treated as interchangeable.\n\n")
  }
}
