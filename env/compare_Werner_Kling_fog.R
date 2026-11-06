# compare_Werner_Kling_fog.R
#
# Simple scatter comparing Werner (MODIS FLCC) and Kling reconstructed fog,
# one point per coordinate x overlapping calendar year (2000-2014 -- the
# only years both products cover), from env_wi_range.csv's per-coordinate
# annual grid.
#
# Not germination-year-specific (per user decision, 2026-10-24): fog_germyr
# and Kling_germyr (the per-sample, germination-year values) only overlap
# for 3 specimens total -- Werner's germyr coverage needs germination year
# >=2009ish, Kling's needs <=2014, and almost no specimen's germination
# year falls in that narrow window with both non-missing. The per-
# coordinate annual grid used here (every coordinate's own fog/Kling value
# in every overlapping calendar year, not tied to any one specimen's
# germination year) is the only version with enough points for a real
# regression + CI. See plot_fog_comparison.R (pheno/) for the original,
# locality-colored version this is adapted from.
#
# Styled for display at small size: black points, lm + CI, large base font.
#
# Units (2026-10-24, per user correction): Werner is mean DAYS PER MONTH
# with fog (FLCC = "fog/low cloud cover days", Jun-Sep), not a percentage --
# earlier scripts in this project mislabeled it "Percent hours with fog";
# fixed here. Kling genuinely is percent HOURS per month with fog (its raw
# 0-1 fraction is rescaled to 0-100 below). The two are therefore not the
# same quantity and cannot be expected to agree in absolute units -- the
# correlation test added below asks only whether they track each other,
# not whether they match in magnitude.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

input_file <- "/private/groups/meyerlab/esteiger/pheno/env_wi_range.csv"
out_file   <- "/private/groups/meyerlab/esteiger/gvtp_popgen/env/compare_Werner_Kling_fog.jpg"
test_out   <- "/private/groups/meyerlab/esteiger/gvtp_popgen/env/compare_Werner_Kling_fog_test.txt"

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
         kling_fog = kling_fog * 100)  # 0-1 fraction -> percent hours per month

cmp <- inner_join(werner_long, kling_long, by = c(id_cols, "year")) %>%
  filter(!is.na(werner_fog), !is.na(kling_fog))

cat(sprintf("Points in comparison: %d (coordinates x years, %d-%d)\n",
            nrow(cmp), min(overlap_years), max(overlap_years)))

# Formal test of association (not just a descriptive r): Pearson correlation
# test, two-sided, on the coordinate x year pairs above.
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
