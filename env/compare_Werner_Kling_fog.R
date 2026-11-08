# compare_Werner_Kling_fog.R
#
# Scatter of Werner (MODIS FLCC) vs. Kling reconstructed fog, one point per
# coordinate x calendar year for 2000-2014 (the years both products cover),
# from the per-coordinate annual grid in env_wi_range.csv. Writes the plot and
# a Pearson correlation test.
#
# Not germination-year-specific: the per-sample germination-year values
# (fog_germyr, Kling_germyr) overlap for only 3 specimens, since Werner needs
# germination year >= ~2009 and Kling <= 2014. The per-coordinate annual grid
# is the only version with enough points for a regression with CI.
#
# Units differ: Werner is mean days per month with fog (Jun-Sep); Kling is
# percent hours per month with fog. The test asks whether they covary, not
# whether they agree in magnitude.
#
# Usage: Rscript compare_Werner_Kling_fog.R

BASE_DIR <- "/private/groups/meyerlab/esteiger"

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

input_file <- file.path(BASE_DIR, "pheno/env_wi_range.csv")
out_file   <- file.path(BASE_DIR, "gvtp_popgen/env/compare_Werner_Kling_fog.jpg")
test_out   <- file.path(BASE_DIR, "gvtp_popgen/env/compare_Werner_Kling_fog_test.txt")

overlap_years <- 2000:2014

wide <- read.csv(input_file, check.names = FALSE)

werner_cols <- paste0("fog_", overlap_years)
kling_cols  <- paste0("kling_fog_", overlap_years)
id_cols     <- c("Locality", "Lat", "Long")

werner_long <- wide[, c(id_cols, werner_cols)] %>%
  pivot_longer(all_of(werner_cols), names_to = "year", values_to = "werner_fog") %>%
  mutate(year = as.integer(sub("fog_", "", year)))

kling_long <- wide[, c(id_cols, kling_cols)] %>%
  pivot_longer(all_of(kling_cols), names_to = "year", values_to = "kling_fog") %>%
  mutate(year = as.integer(sub("kling_fog_", "", year)),
         kling_fog = kling_fog * 100)  # 0-1 fraction to percent

cmp <- inner_join(werner_long, kling_long, by = c(id_cols, "year")) %>%
  filter(!is.na(werner_fog), !is.na(kling_fog))

cat(sprintf("Points in comparison: %d (coordinates x years, %d-%d)\n",
            nrow(cmp), min(overlap_years), max(overlap_years)))

ct <- cor.test(cmp$werner_fog, cmp$kling_fog)
cat(sprintf("Pearson correlation: r = %.3f, 95%% CI [%.3f, %.3f], t(%d) = %.2f, p = %.3g\n",
            ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$parameter, ct$statistic, ct$p.value))

test_summary <- c(
  sprintf("Pearson correlation test: Werner fog (days/month) vs. Kling fog (%% hours/month)"),
  sprintf("n = %d coordinate x year observations, %d-%d", nrow(cmp), min(overlap_years), max(overlap_years)),
  sprintf("r = %.3f, 95%% CI [%.3f, %.3f]", ct$estimate, ct$conf.int[1], ct$conf.int[2]),
  sprintf("t(%d) = %.2f, p = %.3g", ct$parameter, ct$statistic, ct$p.value)
)
writeLines(test_summary, test_out)
cat(sprintf("Saved: %s\n", test_out))

p <- ggplot(cmp, aes(x = werner_fog, y = kling_fog)) +
  geom_point(color = "black", alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", color = "black", fill = "grey60", linewidth = 1) +
  labs(x = "Days per month\nwith fog (Werner)", y = "Percent hours with fog\n(Kling)") +
  theme_bw(base_size = 24) +
  theme(panel.grid.minor = element_blank())

ggsave(out_file, p, width = 6, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", out_file))
