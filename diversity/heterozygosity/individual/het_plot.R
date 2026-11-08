#!/usr/bin/env Rscript
# Individual heterozygosity from per-sample realSFS .ml files: coverage correction,
# per-locality plots, and site-level models against fog, climate water deficit and Ne.
# Run from the directory holding the *.ml files. Writes het.txt, plots, het_statistics_report.txt.
# R library: BASE_DIR/software/R_libs

BASE_DIR <- "/private/groups/meyerlab/esteiger"

.libPaths(file.path(BASE_DIR, "software/R_libs"))

library(ggplot2)
library(ggbeeswarm)
library(dunn.test)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(scales)


gvtp_theme <- theme_bw() + 
  theme(
    text = element_text(size = 70),
    axis.text = element_text(angle = 45),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm"),
    legend.box.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
)


popcolors <- c("CDFW" = "#e69f00", "CDFW_MR" = "#e69f00", "Strauss"="#009e73", "Point_Conception"="#0072b2", "VAFB_Lion_Head"="#f0e442", "VAFB_Arguello"="#d55e00", "Dangermond"="#56b4e9")         # Bang Wang color palette

# Axis labels as "5.0 x 10^-3" plotmath expressions
sci_axis_labels <- function(x) {
  lab <- sprintf("%.1e", x)
  lab <- sub("e([+-])0*(\\d)", "e\\1\\2", lab)
  parse(text = gsub("e", " %*% 10^", lab))
}
			           
files <- list.files(pattern = "*.ml")


het<-data.frame(matrix(ncol=2, nrow=length(files)))
names(het)<-c("Sample", "Het")
het$Sample<-files


for (x in 1:nrow(het)){
a<-scan(files[x])
het$Het[x]<-a[2]/sum(a)
}


write.table(het, "het.txt", row.names=F, col.names=F, quote=F)


# Read in and merge with sample data
sample<-read.table(file.path(BASE_DIR, "field_data.txt"), header=T, sep="\t")
het$Sample <- gsub("_rmdup\\.ml$", "", het$Sample)
het<-merge(het, sample, all.x=T, by.x="Sample", by.y="SampleID")
het$Clusters<-factor(het$Clusters, levels=c("VAFB_Lion_Head", "VAFB_Arguello", "Dangermond", "Point_Conception", "Strauss", "CDFW"))
het$log<-log10(het$Cov)

# Exclude 4531: nearly twice the heterozygosity of the next-highest individual.
# Its coverage is high (~29x), so this is likely contamination or sample mixing.
het<-het[!(het$Sample %in% c("4531")),]


# ---------- Coverage correction ----------

jpeg("He_boxplot_by_location.jpg", height=1000, width=1000)
ggplot(het, aes(x=Clusters, y=Het))+
geom_boxplot(col="black", aes(fill=Clusters), outlier.shape = NA) +
geom_quasirandom(width = 0.2, alpha = 0.5, size = 2) +
scale_x_discrete(labels=c("CDFW"="CDFW", "Dangermond"="Dangermond", "Point_Conception"="Pt. Conception", "Strauss"="Strauss", "VAFB_Arguello"="Pt. Arguello", "VAFB_Lion_Head"="Lion's Head"))+
scale_fill_manual(values=popcolors, guide="none") +
labs(y="Individual He")+
gvtp_theme + 
theme(axis.text = element_text(vjust=0, hjust=1), 
axis.title.x=element_blank())
dev.off()

# Het vs. coverage saturates rather than rising linearly (steep increase at
# low coverage from genotype-likelihood allelic dropout, then a plateau once
# coverage is sufficient for reliable heterozygote calling), so we fit a
# self-starting asymptotic exponential (SSasymp). Fitting against raw Cov
# rather than log10(Cov) tracks the plateau more closely
# (AIC -9799 vs. -9720 for the log10(Cov) version).
mod <- nls(Het ~ SSasymp(Cov, Asym, R0, lrc), data = het)
summary(mod)

het$fitted <- predict(mod)

# Coverage vs. het, with fitted asymptotic curve overlaid
jpeg('het_covg.jpg', units="in", height=5, width=6, res=300)
ggplot() +
  geom_point(data=het, colour="black", size=1.5, shape=21, aes(x=log10(Cov), y=Het, fill=Clusters))+
  geom_line(data=het, aes(x=log10(Cov), y=fitted), color="black", linewidth=1) +
  labs(x="log10(coverage)", y="Nucleotide diversity")+
  scale_fill_manual(values=popcolors, guide="none") +
  theme_bw(base_size=12)
dev.off()

het$residual <- residuals(mod)
inbmed<-het

# Corrected cov vs het
jpeg('het_covg_corrected.jpg', units="in", height=5, width=6, res=300)
ggplot() +
  geom_point(data=inbmed, colour="black", size=1.5, shape=21, aes(x=log10(Cov), y=residual, fill=Clusters))+ 
  labs(x="log10(coverage)", y="Nucleotide diversity")+
  scale_fill_manual(values=popcolors, guide="none") +
  theme_bw(base_size=12)
dev.off()

inbmed$Year<-as.factor(inbmed$Year)
inbmed$Clusters<-factor(inbmed$Clusters, levels=c("VAFB_Lion_Head", "VAFB_Arguello", "Dangermond", "Point_Conception", "Strauss", "CDFW"))

# Corrected He by locality
jpeg("het_boxplot_by_location_correction.jpg", height=1200, width=800)
ggplot(inbmed[inbmed$Year=="2020",], aes(x=Clusters, y=residual))+
geom_boxplot(col="black", aes(fill=Clusters), outlier.shape = NA) +
geom_quasirandom(width = 0.2, alpha = 0.5, size = 2) +
scale_x_discrete(labels=c("CDFW"="CDFW", "Dangermond"="Dangermond", "Point_Conception"="Pt. Conception", "Strauss"="Tranquillon Mtn", "VAFB_Arguello"="Pt. Arguello", "VAFB_Lion_Head"="Lion's Head"))+
scale_fill_manual(values=popcolors, guide="none") +
labs(y="Residual nucleotide diversity")+
gvtp_theme + 
theme(axis.text = element_text(vjust=0, hjust=1), 
axis.title.x=element_blank(),
axis.title.y = element_text(vjust = 0.8))
dev.off()

jpeg("HeResid_boxplot_by_location_forPub.jpg", height=1200, width=1000)
ggplot(inbmed[inbmed$Year=="2020",], aes(x=Clusters, y=residual))+
geom_boxplot(col="black", aes(fill=Clusters),  outlier.size = 4, lwd = 1.5) +
scale_x_discrete(labels=c("CDFW"="CDFW", "Dangermond"="Dangermond", "Point_Conception"="Pt. Conception", "Strauss"="Strauss", "VAFB_Arguello"="Pt. Arguello", "VAFB_Lion_Head"="Lion's Head"))+
scale_fill_manual(values=popcolors, guide="none") +
labs(y="Residual nucleotide diversity")+
scale_y_continuous(labels = sci_axis_labels, breaks = c(0, 0.0025, 0.005)) +
gvtp_theme +
theme(
axis.text.x=element_blank(),
axis.text.y=element_text(angle=45, hjust=1),
axis.title.x=element_blank(),
axis.title.y = element_text(vjust = 0.8, margin = margin(r = -10)))
dev.off()

dunn.test(inbmed$residual, inbmed$Clusters, method = "bonferroni")



# ---------- Environmental variables ----------
env<-read.table(file.path(BASE_DIR, "env_var.txt"), header=T, sep="\t")
inbmed$Sample <- trimws(as.character(inbmed$Sample))
env$SampleID    <- trimws(as.character(env$Sample))
env<-merge(inbmed, env, by.x="Sample", by.y="SampleID")
env<-env[env$Year.x=="2020",]
edna <- read.table(file.path(BASE_DIR, "eDNA_env_data.csv"), header=T, sep=",")
edna$SampleID <- trimws(as.character(edna$SampleID))
env <- merge(env, edna[, c("SampleID", "terraclimate_climate_water_deficit_mm_30y",
                                       "terraclimate_precip_warmest_quarter_30y")],
             by.x="Sample", by.y="SampleID", all.x=TRUE)
names(env)[names(env)=="terraclimate_climate_water_deficit_mm_30y"]   <- "deficit"
names(env)[names(env)=="terraclimate_precip_warmest_quarter_30y"]     <- "precip_warm"



# Order clusters by mean fog and color boxes by mean fog
fog_order <- aggregate(fog ~ Clusters, data=env, FUN=mean)
fog_order <- fog_order[order(fog_order$fog), ]
env$Clusters_fog <- factor(env$Clusters, levels=fog_order$Clusters)
env$mean_fog <- fog_order$fog[match(env$Clusters, fog_order$Clusters)]

jpeg("het_by_cluster_fog_ordered.jpg", height=1200, width=1100)
ggplot(env, aes(x=Clusters_fog, y=residual)) +
geom_boxplot(aes(fill=mean_fog), col="black", outlier.shape=NA) +
geom_quasirandom(width=0.2, alpha=0.5, size=2) +
scale_x_discrete(labels=c("CDFW"="Gaviota", "Dangermond"="Dangermond", "Point_Conception"="Pt. Conception", "Strauss"="Tranquillon Mtn", "VAFB_Arguello"="Pt. Arguello", "VAFB_Lion_Head"="Lion's Head")) +
scale_fill_viridis_c(direction=-1, name="Mean fog\n(hrs/day)") +
labs(y="Residual nucleotide diversity") +
gvtp_theme +
theme(axis.text=element_text(vjust=0, hjust=1),
      axis.title.x=element_blank(),
      legend.text=element_text(size=40),
      legend.title=element_text(size=40))
dev.off()


# Site means: residual He and climate water deficit, one row per locality (n=6)
site_means <- aggregate(cbind(residual, deficit) ~ Clusters, data=env, FUN=mean)

cat("=== Deficit ~ residual He (site level, n=6) ===\n")
mod_def_site <- lm(residual ~ deficit, data=site_means)
summary(mod_def_site)


# ---- Plot: site-mean residual He vs climate water deficit (n=6) -------------
jpeg("het_vs_deficitt.jpg", height=1000, width=1000)
ggplot(data=site_means) +
geom_point(aes(x=deficit, y=residual, fill=Clusters), color="black", size=8, shape=21) +
labs(x="Climate water deficit (30-yr mean, mm)", y="Residual nucleotide diversity") +
scale_fill_manual(values=popcolors, guide="none") +
gvtp_theme + theme(axis.text.x = element_text(vjust=1, hjust=1),
axis.title.y = element_text(hjust=-0.5))
dev.off()


# Site means: residual He and BCM (30-yr mean), one row per locality (n=6)
site_means_bcm <- aggregate(cbind(residual, BCM_30y) ~ Clusters, data=env, FUN=mean)

cat("=== BCM_30y ~ residual He (site level, n=6) ===\n")
mod_bcm_site <- lm(residual ~ BCM_30y, data=site_means_bcm)
summary(mod_bcm_site)




# ---- Plot: site-mean residual He vs BCM (n=6) --------------------------------
jpeg("het_vs_BCM.jpg", height=1000, width=1000)
ggplot(data=site_means_bcm) +
geom_point(aes(x=BCM_30y, y=residual, fill=Clusters), color="black", size=8, shape=21) +
labs(x="BCM (30-yr mean)", y="Residual nucleotide diversity") +
scale_fill_manual(values=popcolors, guide="none") +
gvtp_theme + theme(axis.text.x = element_text(vjust=1, hjust=1),
axis.title.y = element_text(hjust=-0.5))
dev.off()


# ==============================================================================
# SITE-LEVEL ANALYSIS: DISENTANGLING FOG FROM POPULATION HISTORY (LATITUDE)
# ==============================================================================
# env has per-individual rows; aggregate to 6-row site table.
# Lat may appear as "Lat" or "Lat.x" depending on merge column overlap.
library(vegan)

lat_col   <- if ("Lat" %in% names(env)) "Lat" else "Lat.x"
site_full <- aggregate(cbind(residual, fog, deficit, precip_warm) ~ Clusters, data = env, FUN = mean)
site_lat  <- aggregate(env[[lat_col]] ~ env$Clusters, FUN = mean)
names(site_lat) <- c("Clusters", "Lat")
site_full <- merge(site_full, site_lat, by = "Clusters")

# Add effective population size (Ne) estimated externally for each site
Ne_vals <- c(
  CDFW             =  42,
  Dangermond       = 165,
  Point_Conception = 827,
  Strauss          = 712,
  VAFB_Arguello    = 165,
  VAFB_Lion_Head   = 455
)
site_full$Ne <- Ne_vals[as.character(site_full$Clusters)]

site_labels <- c(
  "CDFW"             = "CDFW",
  "Dangermond"       = "Dangermond",
  "Point_Conception" = "Pt. Conception",
  "Strauss"          = "Strauss",
  "VAFB_Arguello"    = "Pt. Arguello",
  "VAFB_Lion_Head"   = "Lion's Head"
)

# ---- 1. Latitude-fog correlation and plot ------------------------------------
lat_fog_cor <- cor.test(site_full$Lat, site_full$fog, method = "pearson")
cat(sprintf("\n=== Latitude-fog Pearson correlation (n=6 sites) ===\n r = %.3f, t = %.3f, p = %.4f\n",
            lat_fog_cor$estimate, lat_fog_cor$statistic, lat_fog_cor$p.value))

jpeg("lat_vs_fog.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(site_full, aes(x = Lat, y = fog, fill = Clusters)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("r = %.3f\np = %.4f",
                           lat_fog_cor$estimate, lat_fog_cor$p.value),
           size = 7) +
  labs(x = "Latitude (degrees N)", y = "Mean fog (days/month)") +
  theme_bw(base_size = 24)
dev.off()
cat("Latitude-fog plot written to lat_vs_fog.jpg\n")


# ---- 2. Ne-fog correlation and Ne ~ He model ---------------------------------
Ne_fog_cor <- cor.test(site_full$Ne, site_full$fog, method = "pearson")
cat(sprintf("\n=== Ne-fog Pearson correlation (n=6 sites) ===\n r = %.3f, t = %.3f, p = %.4f\n",
            Ne_fog_cor$estimate, Ne_fog_cor$statistic, Ne_fog_cor$p.value))

jpeg("Ne_vs_fog.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(site_full, aes(x = Ne, y = fog, fill = Clusters)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("r = %.3f\np = %.4f",
                           Ne_fog_cor$estimate, Ne_fog_cor$p.value),
           size = 7) +
  labs(x = "Effective population size (Ne)", y = "Mean fog (days/month)") +
  theme_bw(base_size = 24)
dev.off()
cat("Ne-fog plot written to Ne_vs_fog.jpg\n")


# ---- 3. Site-level regressions: fog alone, Ne alone, fog + Ne ----------------
mod_site_fog  <- lm(residual ~ fog,      data = site_full)
mod_site_Ne   <- lm(residual ~ Ne,       data = site_full)
mod_site_both <- lm(residual ~ fog + Ne, data = site_full)

jpeg("het_vs_fog.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(site_full, aes(x = fog, y = residual, fill = Clusters)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "Mean fog (days/month)", y = "Residual nucleotide diversity") +
  theme_bw(base_size = 24)
dev.off()
cat("het vs fog plot written to het_vs_fog.jpg\n")

# ---- Scatterplot: individual residual He vs site-mean fog (all individuals) --
jpeg("het_vs_fog_individuals.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(env, aes(x = fog, y = residual, fill = Clusters)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 2, shape = 21, alpha = 0.7) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "Mean fog (days/month)", y = "Residual nucleotide diversity") +
  theme_bw(base_size = 24)
dev.off()
cat("Individual-level het vs fog plot written to het_vs_fog_individuals.jpg\n")


# ==============================================================================
# INDIVIDUAL-LEVEL PLOT: fog vs residual He, per-cluster smooth lines
# GLM: beta family (Het bounded [0,1]), logit link, fog + Clusters
# Smooth lines shown only for clusters with significant fog slope (p < 0.05)
# ==============================================================================

# Beta regression requires Het strictly in (0,1); apply small-sample squish
n_env       <- nrow(env)
env$Het_adj <- (env$Het * (n_env - 1) + 0.5) / n_env

# ==============================================================================
# COMPETING BETA GLMs: Het_adj ~ [env_var] + Clusters
# All models fit on the same complete cases so AIC values are comparable
# deficit = terraclimate_climate_water_deficit_mm_30y (renamed above)
# ==============================================================================
glm_candidates_het <- c(
  "fog", "BIO15", "CHILI", "slope", "deficit", "PercVolSoilH20_10cm", "BIO18"
)

env_glm <- env[complete.cases(env[, c("Het_adj", "Clusters", glm_candidates_het)]), ]

cat(sprintf("\nIndividual beta GLM dataset: %d samples (complete cases, all candidates)\n",
            nrow(env_glm)))

# Fit each candidate beta GLM
beta_fits <- lapply(setNames(glm_candidates_het, glm_candidates_het), function(v) {
  glmmTMB(as.formula(paste("Het_adj ~", v, "+ Clusters")),
          family = beta_family(link = "logit"),
          data   = env_glm)
})

# AIC comparison table
aic_tab_het <- data.frame(
  predictor = glm_candidates_het,
  AIC       = round(sapply(beta_fits, AIC), 2),
  row.names = NULL
)
aic_tab_het           <- aic_tab_het[order(aic_tab_het$AIC, na.last = TRUE), ]
aic_tab_het$delta_AIC <- round(aic_tab_het$AIC - min(aic_tab_het$AIC, na.rm = TRUE), 2)

cat("\n=== Competing Beta GLMs: Het_adj ~ [env] + Clusters — AIC table ===\n")
print(aic_tab_het, row.names = FALSE)

cat("\n=== Top 3 models — summaries ===\n")
for (i in seq_len(3)) {
  v <- aic_tab_het$predictor[i]
  cat(sprintf("\n--- Rank %d: Het_adj ~ %s + Clusters  (AIC = %.2f, dAIC = %.2f) ---\n",
              i, v, aic_tab_het$AIC[i], aic_tab_het$delta_AIC[i]))
  print(summary(beta_fits[[v]]))
}

# Retain fog model for downstream per-cluster analysis and report
mod_het_fog_beta <- beta_fits[["fog"]]

mod_het_fog_beta_noclust <- glmmTMB(Het_adj ~ fog,
                                     family = beta_family(link = "logit"),
                                     data   = env_glm)
cat("\n=== Beta GLM: Het ~ fog only (individual-level, 2020) ===\n")
print(summary(mod_het_fog_beta_noclust))

# ---- Per-cluster fog slopes to determine which clusters get smoothed lines ---
clusters_testable <- names(which(table(env$Clusters) > 2))
fog_slope_pvals_het <- sapply(clusters_testable, function(cl) {
  d <- env[env$Clusters == cl, ]
  if (var(d$fog, na.rm = TRUE) == 0) return(NA_real_)
  fit <- lm(residual ~ fog, data = d)
  s   <- coef(summary(fit))
  if (!"fog" %in% rownames(s)) return(NA_real_)
  s["fog", "Pr(>|t|)"]
})
cat("\nPer-cluster fog slope p-values (lm on residual He):\n")
print(round(fog_slope_pvals_het, 4))

sig_clusters_het <- names(fog_slope_pvals_het)[
  !is.na(fog_slope_pvals_het) & fog_slope_pvals_het < 0.05
]
cat("Significant clusters (p < 0.05):",
    if (length(sig_clusters_het)) paste(sig_clusters_het, collapse = ", ") else "none", "\n")

# ---- Build plot --------------------------------------------------------------
clusters_with_lines <- setdiff(levels(env$Clusters), "CDFW")

p_het_fog_ind <- ggplot(env, aes(x = fog, y = residual)) +
  geom_point(aes(fill = Clusters), color = "black", size = 2, shape = 21, alpha = 0.7) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "Mean fog (days/month)", y = "Residual nucleotide diversity") +
  theme_bw(base_size = 24)

for (cl in clusters_with_lines) {
  p_het_fog_ind <- p_het_fog_ind +
    geom_smooth(data        = env[env$Clusters == cl, ],
                aes(x       = fog, y = residual),
                method      = "lm", se = TRUE,
                color       = popcolors[cl],
                fill        = popcolors[cl],
                linetype    = "dashed",
                alpha       = 0.15,
                linewidth   = 0.9,
                inherit.aes = FALSE)
}

jpeg("het_vs_fog_individuals_bycluster.jpg", units = "in", height = 5, width = 5, res = 300)
print(p_het_fog_ind)
dev.off()
cat("Per-cluster fog plot written to het_vs_fog_individuals_bycluster.jpg\n")


cat("\n=== Site-level model: residual He ~ fog (n=6) ===\n")
print(summary(mod_site_fog))
cat("\n=== Site-level model: residual He ~ Ne (n=6) ===\n")
print(summary(mod_site_Ne))
cat("\n=== Site-level model: residual He ~ fog + Ne (n=6) ===\n")
print(summary(mod_site_both))
cat("\n=== AIC comparison: fog models (lower = better fit) ===\n")
print(AIC(mod_site_fog, mod_site_Ne, mod_site_both))


# ---- 4. Variance partitioning of residual He: fog vs. Ne --------------------
# Fractions: [a] unique fog, [b] shared, [c] unique Ne, [d] unexplained
vp <- varpart(site_full$residual, ~ fog, ~ Ne, data = site_full)
cat("\n=== Variance partitioning: fog vs. Ne -> residual He ===\n")
print(vp)

jpeg("varpart_fog_Ne.jpg", units = "in", height = 5, width = 5, res = 300)
plot(vp, digits = 2, Xnames = c("Fog", "Ne"),
     bg = c("steelblue", "tomato"),
     main = "Variance partitioning of residual He\n(n = 6 sites)")
dev.off()
cat("Variance partitioning plot written to varpart_fog_Ne.jpg\n")


# ==============================================================================
# SITE-LEVEL ANALYSIS: DISENTANGLING CLIMATE WATER DEFICIT FROM Ne
# ==============================================================================

# ---- 5. Ne-deficit correlation and plot --------------------------------------
Ne_def_cor <- cor.test(site_full$Ne, site_full$deficit, method = "pearson")
cat(sprintf("\n=== Ne-deficit Pearson correlation (n=6 sites) ===\n r = %.3f, t = %.3f, p = %.4f\n",
            Ne_def_cor$estimate, Ne_def_cor$statistic, Ne_def_cor$p.value))

jpeg("Ne_vs_deficit.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(site_full, aes(x = Ne, y = deficit, fill = Clusters)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("r = %.3f\np = %.4f",
                           Ne_def_cor$estimate, Ne_def_cor$p.value),
           size = 7) +
  labs(x = "Effective population size (Ne)", y = "Climate water deficit (30-yr mean, mm)") +
  theme_bw(base_size = 24)
dev.off()
cat("Ne-deficit plot written to Ne_vs_deficit.jpg\n")


# ---- Added-variable plot: He vs fog, controlling for Ne ---------------------
e_He_Ne  <- residuals(lm(residual ~ Ne, data = site_full))
e_fog_Ne <- residuals(lm(fog      ~ Ne, data = site_full))

av_data <- data.frame(
  e_fog    = e_fog_Ne,
  e_He     = e_He_Ne,
  Clusters = site_full$Clusters
)

av_mod <- lm(e_He ~ e_fog, data = av_data)

jpeg("added_variable_fog_controlNe.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(av_data, aes(x = e_fog, y = e_He, fill = Clusters)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "Fog | Ne",
       y = "Residual He | Ne") +
  theme_bw(base_size = 24)
dev.off()
cat("Added-variable plot written to added_variable_fog_controlNe.jpg\n")


# ---- 6. Site-level regressions: deficit alone, Ne alone, deficit + Ne --------
mod_site_def      <- lm(residual ~ deficit,      data = site_full)
mod_site_def_both <- lm(residual ~ deficit + Ne, data = site_full)

cat("\n=== Site-level model: residual He ~ deficit (n=6) ===\n")
print(summary(mod_site_def))
cat("\n=== Site-level model: residual He ~ deficit + Ne (n=6) ===\n")
print(summary(mod_site_def_both))
cat("\n=== AIC comparison: deficit models (lower = better fit) ===\n")
print(AIC(mod_site_def, mod_site_Ne, mod_site_def_both))


# ---- 7a. Precip warmest quarter: regression and plot -------------------------
Ne_precip_warm_cor <- cor.test(site_full$Ne, site_full$precip_warm, method = "pearson")
cat(sprintf("\n=== Ne-precip_warm Pearson correlation (n=6 sites) ===\n r = %.3f, t = %.3f, p = %.4f\n",
            Ne_precip_warm_cor$estimate, Ne_precip_warm_cor$statistic, Ne_precip_warm_cor$p.value))

mod_site_precip_warm <- lm(residual ~ precip_warm, data = site_full)
cat("\n=== Site-level model: residual He ~ precip_warmest_quarter_30y (n=6) ===\n")
print(summary(mod_site_precip_warm))

jpeg("het_vs_precip_warm.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(site_full, aes(x = precip_warm, y = residual, fill = Clusters)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("r = %.3f\np = %.4f",
                           sqrt(summary(mod_site_precip_warm)$r.squared) *
                             sign(coef(mod_site_precip_warm)["precip_warm"]),
                           summary(mod_site_precip_warm)$coefficients["precip_warm", 4]),
           size = 7) +
  labs(x = "Precip, warmest quarter (30-yr mean, mm)", y = "Residual nucleotide diversity") +
  theme_bw(base_size = 24)
dev.off()
cat("Precip warmest quarter plot written to het_vs_precip_warm.jpg\n")

# ---- 7. Variance partitioning: deficit vs. Ne --------------------------------
vp_def <- varpart(site_full$residual, ~ deficit, ~ Ne, data = site_full)
cat("\n=== Variance partitioning: deficit vs. Ne -> residual He ===\n")
print(vp_def)

jpeg("varpart_deficit_Ne.jpg", units = "in", height = 5, width = 5, res = 300)
plot(vp_def, digits = 2, Xnames = c("Deficit", "Ne"),
     bg = c("steelblue", "tomato"),
     main = "Variance partitioning of residual He\n(n = 6 sites)")
dev.off()
cat("Variance partitioning plot written to varpart_deficit_Ne.jpg\n")


# ==============================================================================
# STATISTICS REPORT
# ==============================================================================
report_file <- "het_statistics_report.txt"
sep         <- paste(rep("=", 65), collapse = "")

write_section <- function(title, content) {
  cat(sep, "\n", title, "\n", sep, "\n", sep = "", file = report_file, append = TRUE)
  cat(paste(content, collapse = "\n"), "\n\n", file = report_file, append = TRUE)
}

cat("He Analysis - Statistical Report\n",
    "Generated: ", format(Sys.time()), "\n\n",
    sep = "", file = report_file)

write_section(
  "1. Asymptotic Model (SSasymp): nucleotide diversity ~ Coverage [all samples]",
  capture.output(summary(mod))
)

write_section(
  "2. Dunn Test: Residual He by Population Cluster (Bonferroni correction)",
  capture.output(dunn.test(inbmed$residual, inbmed$Clusters, method = "bonferroni"))
)

write_section(
  "3. Latitude-fog Pearson correlation (n=6 sites)",
  sprintf("r = %.3f, t = %.3f, df = %d, p = %.4f",
          lat_fog_cor$estimate, lat_fog_cor$statistic,
          lat_fog_cor$parameter, lat_fog_cor$p.value)
)

write_section(
  "4. Ne-fog Pearson correlation (n=6 sites)",
  sprintf("r = %.3f, t = %.3f, df = %d, p = %.4f",
          Ne_fog_cor$estimate, Ne_fog_cor$statistic,
          Ne_fog_cor$parameter, Ne_fog_cor$p.value)
)

write_section(
  "5. Site-level models: residual He ~ fog (n=6)",
  c(capture.output(summary(mod_site_fog)),
    "", "AIC comparison (fog / Ne / fog+Ne):",
    capture.output(AIC(mod_site_fog, mod_site_Ne, mod_site_both)))
)

write_section(
  "6. Site-level models: residual He ~ Ne (n=6)",
  capture.output(summary(mod_site_Ne))
)

write_section(
  "7. Site-level models: residual He ~ fog + Ne (n=6)",
  capture.output(summary(mod_site_both))
)

write_section(
  "8. Variance partitioning: fog vs. Ne -> residual He",
  capture.output(print(vp))
)

write_section(
  "9. Site-level model: residual He ~ precip_warmest_quarter_30y (n=6)",
  c(sprintf("Ne-precip_warm correlation: r = %.3f, p = %.4f",
            Ne_precip_warm_cor$estimate, Ne_precip_warm_cor$p.value),
    "", capture.output(summary(mod_site_precip_warm)))
)

write_section(
  "10. Ne-deficit Pearson correlation (n=6 sites)",
  sprintf("r = %.3f, t = %.3f, df = %d, p = %.4f",
          Ne_def_cor$estimate, Ne_def_cor$statistic,
          Ne_def_cor$parameter, Ne_def_cor$p.value)
)

write_section(
  "11. Site-level models: residual He ~ deficit (n=6)",
  c(capture.output(summary(mod_site_def)),
    "", "AIC comparison (deficit / Ne / deficit+Ne):",
    capture.output(AIC(mod_site_def, mod_site_Ne, mod_site_def_both)))
)

write_section(
  "12. Site-level models: residual He ~ deficit + Ne (n=6)",
  capture.output(summary(mod_site_def_both))
)

write_section(
  "13. Variance partitioning: deficit vs. Ne -> residual He",
  capture.output(print(vp_def))
)

write_section(
  "14. Competing Beta GLMs: Het_adj ~ [env] + Clusters (individual-level, 2020; beta_family logit link)",
  c(sprintf("n = %d samples (complete cases across all candidates)", nrow(env_glm)),
    "",
    "Model family: beta_family(link = 'logit') via glmmTMB.",
    "Het is a continuous proportion in (0,1); beta regression respects the bounded",
    "support and the logit link ensures predictions remain within (0,1).",
    "Het_adj = (Het*(n-1)+0.5)/n applied before fitting.",
    "",
    "AIC comparison table:",
    capture.output(print(aic_tab_het, row.names = FALSE)),
    "",
    "Top 3 model summaries:")
)
for (i in seq_len(3)) {
  v <- aic_tab_het$predictor[i]
  write_section(
    sprintf("14.%d Rank %d: Het_adj ~ %s + Clusters  (AIC = %.2f, dAIC = %.2f)",
            i, i, v, aic_tab_het$AIC[i], aic_tab_het$delta_AIC[i]),
    capture.output(print(summary(beta_fits[[v]])))
  )
}

write_section(
  "15. Beta GLM: Het ~ fog + Clusters vs fog only (individual-level, 2020)",
  c(capture.output(print(summary(mod_het_fog_beta))),
    "",
    "Per-cluster fog slope p-values (lm on residual He):",
    capture.output(print(round(fog_slope_pvals_het, 4))),
    sprintf("Significant clusters (p < 0.05): %s",
            if (length(sig_clusters_het)) paste(sig_clusters_het, collapse = ", ") else "none"),
    "",
    capture.output(print(summary(mod_het_fog_beta_noclust))),
    "",
    sprintf("AIC comparison (fog only vs fog + Clusters): %.2f vs %.2f",
            AIC(mod_het_fog_beta_noclust), AIC(mod_het_fog_beta)))
)

cat("Report written to", report_file, "\n")
