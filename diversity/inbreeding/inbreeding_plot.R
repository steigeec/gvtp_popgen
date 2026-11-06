#!/usr/bin/env Rscript


.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")


library(ggplot2)
library(data.table)
library(sf)
library(ggbeeswarm)
library(patchwork)
library(lmerTest)
library(glmmTMB)
library(ggeffects)
library(dplyr)
library(AICcmodavg)
library(vegan)




# ESTABLISH PLOTTING PARAMS


gvtp_theme <- theme_bw() + 
  theme(
    text = element_text(size = 70), # Sets size for all text elements
    axis.text = element_text(angle = 45, hjust=1), # Only specify unique modifications
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


popcolors <- c("CDFW" = "#e69f00", "CDFW_MR" = "#e69f00", "Strauss"="#009e73", "Point_Conception"="#0072b2", "VAFB_Lion_Head"="#f0e442", "VAFB_Arguello"="#d55e00", "Dangermond"="#56b4e9", "Between"="black")         # Bang Wang color palette


location_labels <- c(
"CDFW"="Gaviota",             
"Strauss"="Tranquillon Mtn",
"Dangermond"="Dangermond",  
 "Point_Conception"="Pt. Conception",
"VAFB_Lion_Head"="Lion's Head",
"VAFB_Arguello"="Pt. Arguello"
)




# Rebuilt 2026-10-21 (per user decision) from field_data.txt directly, after
# samples were reassigned to a revised set of Strauss/Tranquillon Mtn
# stands. inbmed (post coverage-filter, all years) now only ever has 9
# distinct Strauss OHANA_order values: 72,73,75,78,80,85,87,88,89
# (confirmed empirically, not assumed) — the other 10 stands from the old
# 19-stand set (74,79,81,82,83,86,90,91,392,486) no longer appear at all.
# Per user decision, each surviving stand KEEPS its original color
# (just dropping the entries for stands that no longer exist), rather
# than reassigning fresh colors.
# Re-adjusted 2026-10-22 per user decision: 72/73/75/78 keep their
# original colors; 80 moves from yellow to grey (joining 85); 87 and 89
# move from grey to the darker blue (joining 88).
clustercolors <- c(
  "72" = "#EE8866",
  "73" = "#99DDFF",
  "75" = "#EEDD88", "78" = "#EEDD88",
  "80" = "#DDDDDD", "85" = "#DDDDDD",
  "87" = "#77AADD", "88" = "#77AADD", "89" = "#77AADD"
)






LABELS<-c("Strauss" = "Strauss", "Between" = "Between", "Point_Conception" = "Point Conception", "CDFW" = "CDFW", "VAFB_Lion_Head" = "VAFB Lion's Head", "Dangermond" = "Dangermond", "VAFB_Arguello" = "VAFB Arguello")


inbreeding <- read.table("ngsF_results.txt", header=T)
inbreeding$SampleID<-gsub("_rmdup.bam", "", inbreeding$BAM)
pop<-read.table("/private/groups/meyerlab/esteiger/field_data.txt", header=T, sep="\t")
pop$SampleID<-as.character(pop$SampleID)
inbreeding<-merge(inbreeding, pop, by="SampleID")
inbreeding$log<-log10(inbreeding$Cov)


# Outlier check for sample 4531 within its cluster
{
  s4531 <- inbreeding[inbreeding$SampleID == "4531", ]
  if (nrow(s4531) > 0) {
    cluster_4531 <- s4531$Clusters[1]
    cluster_data <- inbreeding[inbreeding$Clusters == cluster_4531 & inbreeding$SampleID != "4531", "NGSF"]
    q1  <- quantile(cluster_data, 0.25, na.rm=TRUE)
    q3  <- quantile(cluster_data, 0.75, na.rm=TRUE)
    iqr <- q3 - q1
    upper_fence <- q3 + 1.5 * iqr
    lower_fence <- q1 - 1.5 * iqr
    val_4531 <- s4531$NGSF[1]
    is_outlier <- val_4531 > upper_fence | val_4531 < lower_fence
    message(sprintf(
      "Sample 4531 | Cluster: %s | NGSF: %.4f | Cluster Q1=%.4f Q3=%.4f IQR=%.4f | Fences: [%.4f, %.4f] | Outlier: %s",
      cluster_4531, val_4531, q1, q3, iqr, lower_fence, upper_fence, is_outlier
    ))
  } else {
    message("Sample 4531 not found in dataset after merge.")
  }
}


# We decide to remove individual 4531, which has neary twice the heterozygosity of the next-highest individual in the dataset
# 	this individual has high coverage (nearly 29x), so that's not the issue
#	it is potentially impacted by contamination of sample mixing. 
inbreeding<-inbreeding[!(inbreeding$Sample %in% c("4531")),]


# add shapefile for 
# Define the cropping extent
crop_extent <- st_bbox(c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9), crs = st_crs(proj))
# crop the coast off
CA<-st_read("/private/groups/meyerlab/esteiger/env/study_area/GVTP_clip.shp")
CA<-st_make_valid(CA)
crop_extent_sf <- st_as_sfc(st_bbox(crop_extent))
CA <- st_crop(CA, crop_extent)
CA<- st_simplify(CA)


jpeg("log_cov_by_inb.jpg", height=1000, width=2000)
ggplot(data=inbreeding)+
geom_point(aes(x=log10(Cov), y=NGSF, color=Clusters), size=2) +
geom_vline(xintercept=.85)+
labs(x="log10(coverage)", y="Inbreeding coef.")+
scale_color_manual(values=popcolors) +
gvtp_theme
dev.off()

# Just filter out samples below reliable coverage
inbmed <- inbreeding[inbreeding$log > 0.85, ]
inbmed$Year<-as.factor(inbmed$Year)
inbmed$Clusters<-factor(inbmed$Clusters, levels=c("VAFB_Lion_Head", "VAFB_Arguello", "Dangermond", "Point_Conception", "Strauss", "CDFW"))

# Coverage correction model — residuals used downstream for site-level analyses
mod_cov         <- lm(NGSF ~ log10(Cov), data = inbmed)
cat("Coverage correction model:\n"); print(summary(mod_cov))
inbmed$residual <- residuals(mod_cov)


jpeg("log_cov_by_inb_trustworthy.jpg", height=1000, width=2000)
ggplot(data=inbmed)+
geom_point(aes(x=log10(Cov), y=NGSF, color=Clusters), size=2) +
geom_vline(xintercept=.20)+
geom_vline(xintercept=.85)+
labs(x="log10(coverage)", y="Inbreeding coef.")+
scale_color_manual(values=popcolors) +
gvtp_theme
dev.off()


# INBREEDING AS BOXPLOY BY LOCATION WITH CORRECTION
jpeg("inbreeding_boxplot_by_location.jpg", height=1200, width=800)
ggplot(inbmed[inbmed$Year=="2020",], aes(x=Clusters, y=residual))+
geom_boxplot(col="black", aes(fill=Clusters), outlier.shape = NA) +
geom_quasirandom(width = 0.2, alpha = 0.5, size = 2) +
scale_x_discrete(labels=c("CDFW"="CDFW", "Dangermond"="Dangermond", "Point_Conception"="Pt. Conception", "Strauss"="Tranquillon Mtn", "VAFB_Arguello"="Pt. Arguello", "VAFB_Lion_Head"="Lion's Head"))+
scale_fill_manual(values=popcolors, guide="none") +
labs(y="Residual inbreeding coef.")+
gvtp_theme + 
theme(axis.text = element_text(vjust=0, hjust=1), 
axis.title.x=element_blank())
dev.off()


# FOR PUBLICATION VERSION, WITH CORRECTION
jpeg("inbreeding_boxplot_by_location_forPub.jpg", height=1200, width=1000)
ggplot(inbmed[inbmed$Year=="2020",], aes(x=Clusters, y=residual))+
geom_boxplot(col="black", aes(fill=Clusters),  outlier.size = 4, lwd = 1.5) +
scale_fill_manual(values=popcolors, guide="none") +
labs(y="Residual inbreeding coef.")+
gvtp_theme + 
theme(
axis.text.x=element_blank(),
#axis.text = element_text(vjust=0, hjust=1), 
axis.title.x=element_blank())
dev.off()






# boxplots of inbreeding for each year
jpeg("inbreeding_violin_by_year.jpg", height=1100, width=1100)
ggplot(data=inbmed[inbmed$Clusters=="Strauss",],) +
  geom_boxplot(aes(x = Year, y =residual), fill = "gray80", color = "black", outlier.size=5, linewidth=2) +
  labs(x = "Year", y = "Residual inbreeding coef.") +
  gvtp_theme+
  theme(axis.title.y = element_text(hjust = 0.3))
dev.off()


inbmed$Climate <- ifelse(inbmed$Year == "2021", "Drought", "Non_drought")
inbmed$Climate <- factor(inbmed$Climate)


kruskal.test(NGSF ~ Year, data = inbmed)
jpeg("inbreeding_map.jpg", height=1000, width=1700)
ggplot()+
 geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", size = 0.7) +  
  geom_jitter(data=inbmed, aes(x=Long, y=Lat, fill=residual), colour="black", size=3, shape=21, width=.008, height=.008)+
  scale_fill_gradient(low = "white", high = "red",
                       guide = guide_colorbar(title="Residual inbreeding coef.")) +    
coord_sf(crs= 4326) + 
theme_bw(base_size=12) +gvtp_theme
dev.off()






jpeg("inbreeding_map_Strauss.jpg", height=1000, width=1700)
ggplot()+
  geom_jitter(data=inbmed[inbmed$Clusters=="Strauss",], aes(x=Long, y=Lat, fill=residual), colour="black", size=4, shape=21, width=.0008, height=.0008)+
  scale_fill_gradient(low = "white", high = "red",
                       guide = guide_colorbar(title="Residual inbreeding coef.")) +  
coord_sf(crs= 4326) +   
theme_bw(base_size=12) 
dev.off()




# OHANA_order barplot
inbmed$OHANA_order <- factor(inbmed$OHANA_order, levels=c("72", "73", "75", "78", "80", "85", "87", "89", "88"))
jpeg("inbreeding_bars_Strauss.jpg", height=1000, width=1200)
ggplot()+
  geom_boxplot(data=inbmed[inbmed$Clusters=="Strauss",], aes(x=OHANA_order, y=residual, fill=OHANA_order), colour="black", linewidth=2, outlier.size=5)+
  scale_fill_manual(values=clustercolors, guide="none") +
  scale_x_discrete(labels=as.character(seq_along(levels(inbmed$OHANA_order)))) +
  scale_y_continuous(breaks=c(0.00, 0.05, 0.10)) +
labs(x="Stand", y="Residual inbreeding coef.")+
gvtp_theme+
theme(axis.text.x=element_text(angle=0, hjust=0.5),
axis.text.y=element_text(hjust=0))
dev.off()


library(dunn.test)
dunn.test(inbmed$NGSF, inbmed$Clusters, method = "bonferroni")

# Calculate median and IQR per locality
aggregate(NGSF ~ Clusters, data = inbmed, 
          FUN = function(x) c(median = median(x, na.rm = TRUE),
                              IQR = IQR(x, na.rm = TRUE)))
# Does repeated drought drive inbreeding?


env<-read.table("/private/groups/meyerlab/esteiger/env_var.txt", header=T, sep="\t")
inbmed$SampleID <- trimws(as.character(inbmed$SampleID))
env$SampleID    <- trimws(as.character(env$SampleID))
env<-merge(inbmed, env, by="SampleID")
env<-env[env$Year.x=="2020",]

# Also merge eDNA environmental data (many additional env variables per sample)
edna_env <- read.csv("/private/groups/meyerlab/esteiger/eDNA_env_data.csv", header=TRUE)
edna_env$SampleID <- trimws(as.character(edna_env$SampleID))
env <- merge(env, edna_env, by="SampleID", suffixes=c("", "_edna"))


# Summarize to cluster means for env plots
site_env <- env %>%
  group_by(Clusters) %>%
  summarize(
    residual                                    = mean(residual, na.rm = TRUE),
    fog                                         = mean(fog, na.rm = TRUE),
    BCM_30y                                     = mean(BCM_30y, na.rm = TRUE),
    terraclimate_severe_drought_month_count_30y = mean(terraclimate_severe_drought_month_count_30y, na.rm = TRUE),
    terraclimate_palmer_drought_severity_index_30y = mean(terraclimate_palmer_drought_severity_index_30y, na.rm = TRUE),
    terraclimate_soil_moisture_mm_30y           = mean(terraclimate_soil_moisture_mm_30y, na.rm = TRUE)
  )

Ne_vals <- c(
  CDFW             =  42,
  Dangermond       = 165,
  Point_Conception = 827,
  Strauss          = 712,
  VAFB_Arguello    = 165,
  VAFB_Lion_Head   = 455
)
site_env$Ne <- Ne_vals[as.character(site_env$Clusters)]

mod <- lm(NGSF ~ fog + Cov.x + Clusters, data = env)
summary(mod)

# ---- Site-level scatterplots: env variable vs. inbreeding, controlling for Ne ----
# (fog and BCM_30y get their own dedicated added-variable plots further below;
#  this covers the remaining site_env candidates using the same added-variable approach)
site_env_controlNe <- function(var, xlab, outfile) {
  e_inb <- residuals(lm(residual        ~ Ne, data = site_env))
  e_var <- residuals(lm(site_env[[var]] ~ Ne, data = site_env))
  d <- data.frame(e_var = e_var, e_inb = e_inb, Clusters = site_env$Clusters)

  jpeg(outfile, units = "in", height = 5, width = 5, res = 300)
  print(ggplot(d, aes(x = e_var, y = e_inb, fill = Clusters)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_smooth(method = "lm", se = TRUE, color = "grey50", linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
    geom_point(color = "black", size = 4, shape = 21) +
    labs(x = paste0(xlab, " | Ne"), y = "Residual inbreeding coef. | Ne") +
    scale_fill_manual(values = popcolors, guide = "none") +
    theme_bw(base_size = 24))
  dev.off()
  cat(sprintf("Ne-controlled scatter written to %s\n", outfile))
}

site_env_controlNe("terraclimate_severe_drought_month_count_30y",    "Months drought in 30yrs",        "added_variable_inb_drought_controlNe.jpg")
site_env_controlNe("terraclimate_palmer_drought_severity_index_30y", "Palmer drought severity, 30yrs", "added_variable_inb_droughtseverity_controlNe.jpg")
site_env_controlNe("terraclimate_soil_moisture_mm_30y",              "Soil moisture, 30yrs",           "added_variable_inb_soilmoist_controlNe.jpg")


# ==============================================================================
# SITE-LEVEL ANALYSIS: FOG AND Ne AS PREDICTORS OF INBREEDING
# ==============================================================================
# env has per-individual rows filtered to Year 2020 and good coordinates.
# residual column (coverage-corrected NGSF) comes through from inbmed via merge.

lat_col_inb   <- if ("Lat" %in% names(env)) "Lat" else "Lat.x"
site_full_inb <- aggregate(cbind(residual, fog) ~ Clusters, data = env, FUN = mean)
site_bcm_inb <- aggregate(cbind(residual, BCM_30y) ~ Clusters, data = env, FUN = mean)
site_lat_inb  <- aggregate(env[[lat_col_inb]] ~ env$Clusters, FUN = mean)
names(site_lat_inb) <- c("Clusters", "Lat")
site_full_inb <- merge(site_full_inb, site_lat_inb, by = "Clusters")

# Ne_vals already defined above (used for site_env$Ne)
site_full_inb$Ne <- Ne_vals[as.character(site_full_inb$Clusters)]
site_bcm_inb$Ne  <- Ne_vals[as.character(site_bcm_inb$Clusters)]

site_labels_inb <- c(
  "CDFW"             = "CDFW",
  "Dangermond"       = "Dangermond",
  "Point_Conception" = "Pt. Conception",
  "Strauss"          = "Strauss",
  "VAFB_Arguello"    = "Pt. Arguello",
  "VAFB_Lion_Head"   = "Lion's Head"
)

# ---- Ne vs fog scatterplot ---------------------------------------------------
jpeg("Ne_vs_fog_inb.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(site_full_inb, aes(x = fog, y = Ne, fill = Clusters)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "Fog (hrs/day)", y = "Effective population size (Ne)") +
  theme_bw(base_size = 24)
dev.off()
cat("Ne vs fog plot written to Ne_vs_fog_inb.jpg\n")

# ---- Ne-fog correlation ------------------------------------------------------
Ne_fog_cor_inb <- cor.test(site_full_inb$Ne, site_full_inb$fog, method = "pearson")
cat(sprintf("\n=== Ne-fog Pearson correlation (n=6 sites) ===\n r = %.3f, t = %.3f, p = %.4f\n",
            Ne_fog_cor_inb$estimate, Ne_fog_cor_inb$statistic, Ne_fog_cor_inb$p.value))

# ---- Site-level regressions --------------------------------------------------
mod_inb_fog  <- lm(residual ~ fog,      data = site_full_inb)
mod_inb_Ne   <- lm(residual ~ Ne,       data = site_full_inb)
mod_inb_both <- lm(residual ~ fog + Ne, data = site_full_inb)

cat("\n=== Site-level: inbreeding residual ~ fog (n=6) ===\n")
print(summary(mod_inb_fog))
cat("\n=== Site-level: inbreeding residual ~ Ne (n=6) ===\n")
print(summary(mod_inb_Ne))
cat("\n=== Site-level: inbreeding residual ~ fog + Ne (n=6) ===\n")
print(summary(mod_inb_both))
cat("\n=== AIC comparison (lower = better fit) ===\n")
print(AIC(mod_inb_fog, mod_inb_Ne, mod_inb_both))

# ---- Variance partitioning: fog vs. Ne -> inbreeding ------------------------
vp_inb <- varpart(site_full_inb$residual, ~ fog, ~ Ne, data = site_full_inb)
cat("\n=== Variance partitioning: fog vs. Ne -> inbreeding residual ===\n")
print(vp_inb)

jpeg("varpart_inb_fog_Ne.jpg", units = "in", height = 5, width = 5, res = 300)
plot(vp_inb, digits = 2, Xnames = c("Fog", "Ne"),
     bg = c("steelblue", "tomato"),
     main = "Variance partitioning of inbreeding\n(n = 6 sites)")
dev.off()
cat("Variance partitioning plot written to varpart_inb_fog_Ne.jpg\n")

# ---- Added-variable plot: fog -> inbreeding, controlling for Ne --------------
e_inb_Ne     <- residuals(lm(residual ~ Ne, data = site_full_inb))
e_fog_Ne_inb <- residuals(lm(fog      ~ Ne, data = site_full_inb))

av_data_inb <- data.frame(
  e_fog    = e_fog_Ne_inb,
  e_inb    = e_inb_Ne,
  Clusters = site_full_inb$Clusters
)

jpeg("added_variable_inb_fog_controlNe.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(av_data_inb, aes(x = e_fog, y = e_inb, fill = Clusters)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "Fog | Ne",
       y = "Inbreeding | Ne") +
  theme_bw(base_size = 24)
dev.off()
cat("Added-variable plot written to added_variable_inb_fog_controlNe.jpg\n")


# ---- Variance partitioning: BCM (CWD) vs. Ne -> inbreeding ------------------
vp_inb_bcm <- varpart(site_bcm_inb$residual, ~ BCM_30y, ~ Ne, data = site_bcm_inb)
cat("\n=== Variance partitioning: BCM CWD vs. Ne -> inbreeding residual ===\n")
print(vp_inb_bcm)

jpeg("varpart_inb_BCM_Ne.jpg", units = "in", height = 5, width = 5, res = 300)
plot(vp_inb_bcm, digits = 2, Xnames = c("CWD", "Ne"),
     bg = c("steelblue", "tomato"),
     main = "Variance partitioning of inbreeding\n(n = 6 sites)")
dev.off()
cat("Variance partitioning plot written to varpart_inb_BCM_Ne.jpg\n")

# ---- Added-variable plot: BCM (CWD) -> inbreeding, controlling for Ne -------
e_inb_Ne     <- residuals(lm(residual ~ Ne, data = site_bcm_inb))
e_bcm_Ne_inb <- residuals(lm(BCM_30y  ~ Ne, data = site_bcm_inb))

av_data_inb <- data.frame(
  e_bcm    = e_bcm_Ne_inb,
  e_inb    = e_inb_Ne,
  Clusters = site_bcm_inb$Clusters
)

jpeg("added_variable_inb_BCM_controlNe.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(av_data_inb, aes(x = e_bcm, y = e_inb, fill = Clusters)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "CWD | Ne",
       y = "Inbreeding | Ne") +
  theme_bw(base_size = 24)
dev.off()
cat("Added-variable plot written to added_variable_inb_BCM_controlNe.jpg\n")




# ==============================================================================
# SITE-LEVEL ANALYSIS: PRECIP WARMEST QUARTER (BIO18) AND Ne AS PREDICTORS
# ==============================================================================

# Aggregate BIO18 (precipitation of warmest quarter) to site level
site_full_inb$BIO18 <- aggregate(
  BIO18 ~ Clusters, data = env, FUN = mean
)$BIO18[match(site_full_inb$Clusters,
  aggregate(BIO18 ~ Clusters, data = env, FUN = mean)$Clusters)]

# ---- Ne-BIO18 correlation -----------------------------------------------
Ne_BIO18_cor_inb <- cor.test(site_full_inb$Ne, site_full_inb$BIO18, method = "pearson")
cat(sprintf("\n=== Ne-BIO18 (precip, warmest qtr.) Pearson correlation (n=6 sites) ===\n r = %.3f, t = %.3f, p = %.4f\n",
            Ne_BIO18_cor_inb$estimate, Ne_BIO18_cor_inb$statistic, Ne_BIO18_cor_inb$p.value))

# ---- Site-level regressions --------------------------------------------------
mod_inb_BIO18      <- lm(residual ~ BIO18,      data = site_full_inb)
mod_inb_BIO18_both <- lm(residual ~ BIO18 + Ne, data = site_full_inb)

cat("\n=== Site-level: inbreeding residual ~ BIO18 (n=6) ===\n")
print(summary(mod_inb_BIO18))
cat("\n=== Site-level: inbreeding residual ~ BIO18 + Ne (n=6) ===\n")
print(summary(mod_inb_BIO18_both))
cat("\n=== AIC comparison: BIO18 / Ne / BIO18+Ne ===\n")
print(AIC(mod_inb_BIO18, mod_inb_Ne, mod_inb_BIO18_both))

# ---- Variance partitioning: BIO18 vs. Ne --------------------------------
vp_inb_BIO18 <- varpart(site_full_inb$residual, ~ BIO18, ~ Ne, data = site_full_inb)
cat("\n=== Variance partitioning: BIO18 vs. Ne -> inbreeding ===\n")
print(vp_inb_BIO18)

jpeg("varpart_inb_BIO18_Ne.jpg", units = "in", height = 5, width = 5, res = 300)
plot(vp_inb_BIO18, digits = 2, Xnames = c("Precip\n(warmest qtr.)", "Ne"),
     bg = c("steelblue", "tomato"),
     main = "Variance partitioning of inbreeding\n(n = 6 sites)")
dev.off()
cat("Variance partitioning plot written to varpart_inb_BIO18_Ne.jpg\n")

# ---- Added-variable plot: BIO18 -> inbreeding, controlling for Ne -------
e_inb_Ne_BIO18 <- residuals(lm(residual ~ Ne,    data = site_full_inb))
e_BIO18_Ne     <- residuals(lm(BIO18    ~ Ne,    data = site_full_inb))

av_data_inb_BIO18 <- data.frame(
  e_BIO18  = e_BIO18_Ne,
  e_inb    = e_inb_Ne_BIO18,
  Clusters = site_full_inb$Clusters
)

jpeg("added_variable_inb_BIO18_controlNe.jpg", units = "in", height = 5, width = 5, res = 300)
ggplot(av_data_inb_BIO18, aes(x = e_BIO18, y = e_inb, fill = Clusters)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
  geom_point(color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = popcolors, guide = "none") +
  labs(x = "Precip (warmest qtr.) | Ne",
       y = "Residual inbreeding coef. | Ne") +
  theme_bw(base_size = 24)
dev.off()
cat("Added-variable plot written to added_variable_inb_BIO18_controlNe.jpg\n")

# ==============================================================================
# SITE-LEVEL MODEL COMPARISON: fog + BIO18 (precip, warmest qtr.) + Ne
# ------------------------------------------------------------------------------
# Ne is kept as an ordinary fixed-effect covariate, not a random effect: a
# random effect needs a grouping factor with multiple observations per level
# to estimate a variance component, but at the site level (n=6) every
# cluster/site has exactly one row and one Ne value -- there's no within-group
# replication, and Ne is a continuous covariate here, not a grouping factor at
# all. Fitting it as a random-effect grouping term would be uninterpretable
# (and, in an individual-level model, would just be aliased with a per-cluster
# random intercept, since Ne varies only at the cluster level).
#
# With only 6 site-level observations we also can't fit every possible
# interaction: the full fog*BIO18*Ne three-way model has 8 parameters, more
# than there are data points, and is rank-deficient, so it's excluded below.
# We fit the rest -- all three single terms, all three additive pairs, the
# full additive triple, and all three pairwise interactions -- and rank them
# by AICc rather than plain AIC, since AIC is known to be biased toward
# complex models once parameter count approaches sample size.
#
# AICc's own small-sample correction term is undefined once a model's
# coefficient count gets too close to n=6: every 4-coefficient model below
# (the additive triple, and all three pairwise interactions) comes back as
# AICc = Inf. That isn't a bug -- it's AICc correctly refusing to rank models
# we don't have the degrees of freedom to support. The 3-coefficient additive
# pairs do return a finite AICc, but it's already 24-33 points worse than the
# best single-predictor model -- nowhere close to competitive. In practice,
# only the three single-predictor models (fog, BIO18, Ne) are worth comparing
# at this sample size; everything else is included in the table for
# completeness but should not be used to declare a "best" multi-predictor model.
# ==============================================================================

site_models <- list(
  fog            = lm(residual ~ fog,             data = site_full_inb),
  BIO18          = lm(residual ~ BIO18,           data = site_full_inb),
  Ne             = lm(residual ~ Ne,              data = site_full_inb),
  `fog+BIO18`    = lm(residual ~ fog + BIO18,     data = site_full_inb),
  `fog+Ne`       = lm(residual ~ fog + Ne,        data = site_full_inb),
  `BIO18+Ne`     = lm(residual ~ BIO18 + Ne,      data = site_full_inb),
  `fog+BIO18+Ne` = lm(residual ~ fog + BIO18 + Ne, data = site_full_inb),
  `fog*BIO18`    = lm(residual ~ fog * BIO18,     data = site_full_inb),
  `fog*Ne`       = lm(residual ~ fog * Ne,        data = site_full_inb),
  `BIO18*Ne`     = lm(residual ~ BIO18 * Ne,      data = site_full_inb)
)

site_model_tab <- data.frame(
  model = names(site_models),
  df    = sapply(site_models, function(m) length(coef(m))),
  AICc  = sapply(site_models, AICc),
  stringsAsFactors = FALSE
)
site_model_tab            <- site_model_tab[order(site_model_tab$AICc), ]
site_model_tab$delta_AICc <- round(site_model_tab$AICc - min(site_model_tab$AICc), 2)
site_model_tab$AICc       <- round(site_model_tab$AICc, 2)

cat("\n=== Site-level model comparison: fog + BIO18 + Ne, additive & interactive (n=6, AICc) ===\n")
print(site_model_tab, row.names = FALSE)

best_site_model_name <- site_model_tab$model[1]
cat(sprintf("\nBest model by AICc: residual ~ %s\n", best_site_model_name))
print(summary(site_models[[best_site_model_name]]))


# ==============================================================================
# SITE-LEVEL ANALYSIS: fog_germ_yr + BIO12_germ_yr + Ne -> inbreeding
# Cluster means (2020), Ne values from Ne_vals defined above
# ==============================================================================

site_fogprecip <- env %>%
  group_by(Clusters) %>%
  summarize(
    residual        = mean(residual,        na.rm = TRUE),
    fog_germ_yr     = mean(fog_germ_yr,     na.rm = TRUE),
    BIO12_germ_yr = mean(BIO12_germ_yr, na.rm = TRUE)
  )
site_fogprecip$Ne <- Ne_vals[as.character(site_fogprecip$Clusters)]

cat("\n=== Site means: fog_germ_yr, BIO12_germ_yr, Ne (2020) ===\n")
print(as.data.frame(site_fogprecip))

# ---- Variance partitioning: fog_germ_yr + BIO12_germ_yr + Ne --------------
vp_fogprecip <- varpart(site_fogprecip$residual,
                        ~ fog_germ_yr, ~ BIO12_germ_yr, ~ Ne,
                        data = site_fogprecip)
cat("\n=== Variance partitioning: fog_germ_yr + BIO12_germ_yr + Ne -> inbreeding (n=6) ===\n")
print(vp_fogprecip)

jpeg("varpart_inb_foggerm_BIO12germ_Ne.jpg", units = "in", height = 5, width = 5, res = 300)
plot(vp_fogprecip, digits = 2,
     Xnames = c("Fog\n(water yr)", "Precip\n(water yr)", "Ne"),
     bg     = c("steelblue", "tomato", "gold"),
     main   = "Variance partitioning of inbreeding\n(n = 6 sites, 2020)")
dev.off()
cat("Variance partitioning plot written to varpart_inb_foggerm_BIO12germ_Ne.jpg\n")

# ---- Site-level scatterplots, controlling for Ne -----------------------------
site_fogprecip_controlNe <- function(var, xlab, outfile) {
  e_inb <- residuals(lm(residual              ~ Ne, data = site_fogprecip))
  e_var <- residuals(lm(site_fogprecip[[var]] ~ Ne, data = site_fogprecip))
  d <- data.frame(e_var = e_var, e_inb = e_inb, Clusters = site_fogprecip$Clusters)

  jpeg(outfile, units = "in", height = 5, width = 5, res = 300)
  print(ggplot(d, aes(x = e_var, y = e_inb, fill = Clusters)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_smooth(method = "lm", se = TRUE, color = "grey50", linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
    geom_point(color = "black", size = 4, shape = 21) +
    labs(x = paste0(xlab, " | Ne"), y = "Residual inbreeding coef. | Ne") +
    scale_fill_manual(values = popcolors, guide = "none") +
    theme_bw(base_size = 24))
  dev.off()
  cat(sprintf("Ne-controlled scatter written to %s\n", outfile))
}

site_fogprecip_controlNe("fog_germ_yr",   "Fog (hrs/day)",              "added_variable_inb_foggerm_controlNe.jpg")
site_fogprecip_controlNe("BIO12_germ_yr", "Annual precipitation (mm)", "added_variable_inb_BIO12germ_controlNe.jpg")


# ==============================================================================
# BCM_30y AND BCM_germ_yr VS INBREEDING (individual and population level)
# ==============================================================================

env10 <- env[!is.na(env$Cov.x) & env$Cov.x > 20, ]
cat(sprintf("\n>20x env dataset: %d individuals (Year 2020)\n", nrow(env10)))

# NOTE: env variable vs. inbreeding relationships are only shown at site-level scale
# (see the site_env / site_full_inb / site_bcm_inb / site_fogprecip added-variable plots
# above). Strauss is the exception: it has enough individuals sampled across years to
# support individual-level plots, so those are kept below.

# Strauss-only scatter: BCM_30y vs inbreeding
jpeg("inb_vs_BCM_30y_strauss.jpg", units = "in", height = 5, width = 5, res = 300)
print(ggplot(env10[env10$Clusters == "Strauss", ], aes(x = BCM_30y, y = residual)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", fill = popcolors["Strauss"], size = 2, shape = 21, alpha = 0.7) +
  labs(x = "CWD (30y avg)", y = "Inbreeding coef. (residual)") +
  theme_bw(base_size = 24))
dev.off()
cat("Strauss scatter written to inb_vs_BCM_30y_strauss.jpg\n")

# Strauss-only scatter: BCM_germ_yr vs inbreeding
jpeg("inb_vs_BCM_germ_yr_strauss.jpg", units = "in", height = 5, width = 5, res = 300)
print(ggplot(env10[env10$Clusters == "Strauss", ], aes(x = BCM_germ_yr, y = residual)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", fill = popcolors["Strauss"], size = 2, shape = 21, alpha = 0.7) +
  labs(x = "CWD (germ. yr.)", y = "Inbreeding coef. (residual)") +
  theme_bw(base_size = 24))
dev.off()
cat("Strauss scatter written to inb_vs_BCM_germ_yr_strauss.jpg\n")

# Strauss-only scatter: fog (30y avg) vs inbreeding
jpeg("inb_vs_fog_30y_strauss.jpg", units = "in", height = 5, width = 5, res = 300)
print(ggplot(env10[env10$Clusters == "Strauss", ], aes(x = fog, y = residual)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", fill = popcolors["Strauss"], size = 2, shape = 21, alpha = 0.7) +
  labs(x = "Fog (30y avg, hrs/day)", y = "Inbreeding coef. (residual)") +
  theme_bw(base_size = 24))
dev.off()
cat("Strauss scatter written to inb_vs_fog_30y_strauss.jpg\n")

# Strauss-only scatter: fog_germ_yr vs inbreeding
jpeg("inb_vs_fog_germ_yr_strauss.jpg", units = "in", height = 5, width = 5, res = 300)
print(ggplot(env10[env10$Clusters == "Strauss", ], aes(x = fog_germ_yr, y = residual)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", fill = popcolors["Strauss"], size = 2, shape = 21, alpha = 0.7) +
  labs(x = "Fog (germ. yr., hrs/day)", y = "Inbreeding coef. (residual)") +
  theme_bw(base_size = 24))
dev.off()
cat("Strauss scatter written to inb_vs_fog_germ_yr_strauss.jpg\n")

# Strauss-only scatter: BIO12 (annual precip) vs inbreeding
jpeg("inb_vs_BIO12_strauss.jpg", units = "in", height = 5, width = 5, res = 300)
print(ggplot(env10[env10$Clusters == "Strauss", ], aes(x = BIO12, y = residual)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", fill = popcolors["Strauss"], size = 2, shape = 21, alpha = 0.7) +
  labs(x = "Annual precip (30y avg, mm)", y = "Inbreeding coef. (residual)") +
  theme_bw(base_size = 24))
dev.off()
cat("Strauss scatter written to inb_vs_BIO12_strauss.jpg\n")

# Strauss-only scatter: BIO18 (precip, warmest qtr.) vs inbreeding
jpeg("inb_vs_BIO18_strauss.jpg", units = "in", height = 5, width = 5, res = 300)
print(ggplot(env10[env10$Clusters == "Strauss", ], aes(x = BIO18, y = residual)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(color = "black", fill = popcolors["Strauss"], size = 2, shape = 21, alpha = 0.7) +
  labs(x = "Precip, warmest qtr. (30y avg, mm)", y = "Inbreeding coef. (residual)") +
  theme_bw(base_size = 24))
dev.off()
cat("Strauss scatter written to inb_vs_BIO18_strauss.jpg\n")


# ==============================================================================
# STRAUSS-ONLY MODEL COMPARISON: BIO18 + fog + BIO15 predicting residual F
# ------------------------------------------------------------------------------
# Individual-level (n=57, >20x coverage, 2020 only), restricted to Strauss.
# 2020-only, like the site-level analyses above: we checked whether Year needs
# to be controlled for here (random-intercept models across all 3 years), and
# found Year's variance component was estimated at exactly 0 -- singular fit,
# even in a bare intercept-only model with no predictors at all. There's no
# detectable between-year variation in Strauss's residual F to control for, so
# there's no need to bring in the other years or a random effect; plain
# 2020-only fixed-effect lm() is appropriate.
#
# With n=57 there's plenty of residual df even for the fully saturated 3-way
# interaction model (8 parameters), so we fit the full candidate set: all
# three single terms, all three additive pairs, the additive triple, all three
# pairwise-interaction models, the all-pairwise-interactions model (no 3-way),
# and the fully saturated model. Ranked by AICc.
# ==============================================================================

strauss10 <- env10[env10$Clusters == "Strauss", ]
cat(sprintf("\nStrauss individual-level model dataset: n = %d\n", nrow(strauss10)))

strauss_models <- list(
  BIO18                    = lm(residual ~ BIO18,                   data = strauss10),
  fog                      = lm(residual ~ fog,                     data = strauss10),
  BIO15                    = lm(residual ~ BIO15,                   data = strauss10),
  `BIO18+fog`              = lm(residual ~ BIO18 + fog,             data = strauss10),
  `BIO18+BIO15`            = lm(residual ~ BIO18 + BIO15,           data = strauss10),
  `fog+BIO15`               = lm(residual ~ fog + BIO15,             data = strauss10),
  `BIO18+fog+BIO15`        = lm(residual ~ BIO18 + fog + BIO15,     data = strauss10),
  `BIO18*fog`               = lm(residual ~ BIO18 * fog,             data = strauss10),
  `BIO18*BIO15`             = lm(residual ~ BIO18 * BIO15,           data = strauss10),
  `fog*BIO15`               = lm(residual ~ fog * BIO15,             data = strauss10),
  `all pairwise (no 3-way)` = lm(residual ~ (BIO18 + fog + BIO15)^2, data = strauss10),
  `BIO18*fog*BIO15`         = lm(residual ~ BIO18 * fog * BIO15,     data = strauss10)
)

strauss_model_tab <- data.frame(
  model = names(strauss_models),
  df    = sapply(strauss_models, function(m) length(coef(m))),
  AICc  = sapply(strauss_models, AICc),
  stringsAsFactors = FALSE
)
strauss_model_tab            <- strauss_model_tab[order(strauss_model_tab$AICc), ]
strauss_model_tab$delta_AICc <- round(strauss_model_tab$AICc - min(strauss_model_tab$AICc), 2)
strauss_model_tab$AICc       <- round(strauss_model_tab$AICc, 2)

cat("\n=== Strauss-only model comparison: BIO18 + fog + BIO15, additive & interactive (n=57, AICc) ===\n")
print(strauss_model_tab, row.names = FALSE)

best_strauss_model_name <- strauss_model_tab$model[1]
cat(sprintf("\nBest model by AICc: residual ~ %s\n", best_strauss_model_name))
print(summary(strauss_models[[best_strauss_model_name]]))


# ==============================================================================
# INDIVIDUAL-LEVEL COMPETING GLMs: residual ~ [env_var] + Clusters
# Candidate env predictors swapped in one at a time; all models fit on the
# intersection of complete cases so AIC values are directly comparable.
# Uses all sampling years (2020-2022); Clusters absorbs spatial structure.
# ==============================================================================

env_var_all <- read.table("/private/groups/meyerlab/esteiger/env_var.txt", header = T, sep = "\t")
env_var_all$SampleID <- trimws(as.character(env_var_all$SampleID))
edna_env_glm <- read.csv("/private/groups/meyerlab/esteiger/eDNA_env_data.csv", header = TRUE)
edna_env_glm$SampleID <- trimws(as.character(edna_env_glm$SampleID))

inb_all_env <- merge(inbmed, env_var_all, by = "SampleID")
inb_all_env <- merge(inb_all_env, edna_env_glm, by = "SampleID", suffixes = c("", "_edna"))

glm_candidates <- c(
  "fog_germ_yr",
  "BIO15",
  "CHILI",
  "slope",
  "BCM_30y",
  "PercVolSoilH20_10cm",
  "BIO18"
)

# Filter to complete cases across all candidates for valid AIC comparison
inb_glm <- inb_all_env[complete.cases(inb_all_env[, c("residual", "Clusters", glm_candidates)]), ]

cat(sprintf("\nIndividual GLM dataset: %d samples (complete cases, all candidates), %d clusters, %d years\n",
            nrow(inb_glm),
            length(unique(inb_glm$Clusters)),
            length(unique(inb_glm$Year.x))))

# Fit all candidate models
glm_fits <- lapply(setNames(glm_candidates, glm_candidates), function(v) {
  lm(as.formula(paste("residual ~", v, "+ Clusters")), data = inb_glm)
})

# AIC comparison table
aic_tab <- data.frame(
  predictor = glm_candidates,
  AIC       = round(sapply(glm_fits, AIC), 2),
  row.names = NULL
)
aic_tab           <- aic_tab[order(aic_tab$AIC), ]
aic_tab$delta_AIC <- round(aic_tab$AIC - min(aic_tab$AIC), 2)

cat("\n=== Competing GLMs: residual ~ [env] + Clusters — AIC table ===\n")
print(aic_tab, row.names = FALSE)

cat("\n=== Top 3 models — summaries ===\n")
for (i in seq_len(3)) {
  v <- aic_tab$predictor[i]
  cat(sprintf("\n--- Rank %d: residual ~ %s + Clusters  (AIC = %.2f, dAIC = %.2f) ---\n",
              i, v, aic_tab$AIC[i], aic_tab$delta_AIC[i]))
  print(summary(glm_fits[[v]]))
}

# Retain fog_germ_yr model as mod_ind_glm for downstream added-variable and beta GLM sections
mod_ind_glm <- glm_fits[["fog_germ_yr"]]
inb_all_env <- inb_glm
cat("\nAIC (fog_germ_yr model):", AIC(mod_ind_glm), "\n")

# NOTE: individual-level env-vs-inbreeding plot removed here (all clusters combined);
# fog vs. inbreeding, controlling for Ne, is shown at site-level scale further above
# (added_variable_inb_foggerm_controlNe.jpg).

# ==============================================================================
# INDIVIDUAL-LEVEL BETA GLM + PER-CLUSTER FOG PLOT
# Equivalent to het_vs_fog_individuals_bycluster.jpg in het_plot.R
#
# Family: beta_family(link = "logit")
# Rationale: NGSF (inbreeding coefficient) is bounded strictly in [0, 1],
# making the Gaussian assumption of lm inappropriate (residuals are not
# symmetric and variance is not constant near the boundaries). The beta
# distribution is the natural choice for a continuous proportion response.
# The logit link maps the linear predictor from (-Inf, Inf) onto (0, 1),
# matching the support of the beta distribution and providing an
# interpretable log-odds scale for the predictors.
# ==============================================================================

n_inb <- nrow(inb_all_env)
inb_all_env$NGSF_adj <- (inb_all_env$NGSF * (n_inb - 1) + 0.5) / n_inb

mod_inb_fog_beta <- glmmTMB(NGSF_adj ~ fog_germ_yr + Clusters,
                              family = beta_family(link = "logit"),
                              data   = inb_all_env)
cat("\n=== Beta GLM: NGSF ~ fog_germ_yr + Clusters (individual-level, all years) ===\n")
print(summary(mod_inb_fog_beta))

# ---- Per-cluster fog slopes (lm on residual, to decide which clusters get lines)
clusters_testable_inb <- names(which(table(inb_all_env$Clusters) > 2))
fog_slope_pvals_inb <- sapply(clusters_testable_inb, function(cl) {
  d <- inb_all_env[inb_all_env$Clusters == cl, ]
  if (var(d$fog_germ_yr, na.rm = TRUE) == 0) return(NA_real_)
  fit <- lm(residual ~ fog_germ_yr, data = d)
  s   <- coef(summary(fit))
  if (!"fog_germ_yr" %in% rownames(s)) return(NA_real_)
  s["fog_germ_yr", "Pr(>|t|)"]
})
cat("\nPer-cluster fog_germ_yr slope p-values (lm on residual inbreeding):\n")
print(round(fog_slope_pvals_inb, 4))

sig_clusters_inb <- names(fog_slope_pvals_inb)[
  !is.na(fog_slope_pvals_inb) & fog_slope_pvals_inb < 0.05
]
cat("Significant clusters (p < 0.05):",
    if (length(sig_clusters_inb)) paste(sig_clusters_inb, collapse = ", ") else "none", "\n")

# NOTE: individual-level per-cluster fog plot removed here (all clusters combined,
# non-Strauss individuals included); the per-cluster slope p-values above are still
# computed and reported.


# ==============================================================================
# TWO-PANEL PLOT: Soil moisture + Inbreeding by Year, Strauss only
# ==============================================================================

# Build Strauss-only multi-year dataset with soil moisture
strauss_env_raw <- read.table("/private/groups/meyerlab/esteiger/env_var.txt", header=T, sep="\t")
strauss_env_raw$SampleID <- trimws(as.character(strauss_env_raw$SampleID))

strauss_all <- inbmed[inbmed$Clusters == "Strauss", ]
strauss_all$SampleID <- trimws(as.character(strauss_all$SampleID))
strauss_all <- merge(strauss_all, strauss_env_raw, by = "SampleID")


strauss_all$Year    <- as.factor(strauss_all$Year.x)
strauss_all$Climate <- factor(ifelse(strauss_all$Year == "2021", "Drought", "Non_drought"))

# --- Statistical tests: inbreeding ~ Year (Strauss) --------------------------
kw_strauss_year    <- kruskal.test(NGSF ~ Year, data = strauss_all)
dunn_strauss_year  <- dunn.test(strauss_all$NGSF, strauss_all$Year, method = "bonferroni")

# Drought vs non-drought (Wilcoxon rank-sum, Strauss only)
wilcox_drought_inb <- wilcox.test(NGSF ~ Climate, data = strauss_all)

cat("\n=== Kruskal-Wallis: Strauss inbreeding ~ Year ===\n")
print(kw_strauss_year)
cat("\n=== Wilcoxon: Strauss inbreeding, Drought vs Non-drought ===\n")
print(wilcox_drought_inb)

# --- Build plots -------------------------------------------------------------
p_soil <- ggplot(strauss_all, aes(x = Year, y = PercVolSoilH20_10cm)) +
  geom_boxplot(fill = "gray80", color = "black", outlier.size = 5, linewidth = 2) +
  labs(x = "Year", y = "Soil moisture (% vol, 10 cm)") +
  gvtp_theme +
  theme(axis.title.y = element_text(hjust = 0.3))

p_inb_year <- ggplot(strauss_all, aes(x = Year, y = residual)) +
  geom_boxplot(fill = "gray80", color = "black", outlier.size = 5, linewidth = 2) +
  labs(x = "Year", y = "Residual inbreeding coef.") +
  gvtp_theme +
  theme(axis.title.y = element_text(hjust = 0.3))

jpeg("soilmoist_inbreeding_by_year_Strauss.jpg", height = 1100, width = 2200)
print(p_soil + p_inb_year + plot_annotation(tag_levels = "A"))
dev.off()
cat("Two-panel plot written to soilmoist_inbreeding_by_year_Strauss.jpg\n")

# =============================================================================
# OUTPUT REPORT
# =============================================================================

report_lines <- c()
rcat <- function(...) { report_lines <<- c(report_lines, paste0(...)) }

rcat("# Inbreeding Analysis Report")
rcat("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

rcat("---\n")
rcat("## Sample sizes (post-filtering)\n")
rcat("```")
rcat(capture.output(print(table(inbmed$Clusters, inbmed$Year))))
rcat("```\n")

rcat("## Median and IQR of inbreeding coefficient by locality\n")
rcat("```")
rcat(capture.output(
  print(aggregate(NGSF ~ Clusters, data = inbmed,
    FUN = function(x) c(median = median(x, na.rm=TRUE), IQR = IQR(x, na.rm=TRUE))))
))
rcat("```\n")

rcat("## Kruskal-Wallis test: inbreeding ~ Year (all populations)\n")
rcat("```")
rcat(capture.output(print(kruskal.test(NGSF ~ Year, data = inbmed))))
rcat("```\n")

rcat("## Kruskal-Wallis test: Strauss inbreeding ~ Year\n")
rcat("```")
rcat(capture.output(print(kw_strauss_year)))
rcat("```\n")

rcat("## Dunn test: Strauss inbreeding ~ Year (Bonferroni correction)\n")
rcat("```")
rcat(capture.output(dunn_strauss_year))
rcat("```\n")

rcat("## Wilcoxon test: Strauss inbreeding, Drought (2021) vs Non-drought\n")
rcat("```")
rcat(capture.output(print(wilcox_drought_inb)))
rcat(sprintf("n drought = %d, n non-drought = %d",
             sum(strauss_all$Climate == "Drought"),
             sum(strauss_all$Climate == "Non_drought")))
rcat("```\n")

rcat("## Dunn test: inbreeding ~ Clusters (Bonferroni correction)\n")
rcat("```")
rcat(capture.output(dunn.test(inbmed$NGSF, inbmed$Clusters, method = "bonferroni")))
rcat("```\n")

rcat("## Linear model: NGSF ~ fog + Cov + Clusters\n")
rcat("```")
rcat(capture.output(print(summary(mod))))
rcat("```\n")

rcat("## env dataset dimensions (after all merges and filters)\n")
rcat("```")
rcat(capture.output(cat("Rows:", nrow(env), "  Columns:", ncol(env), "\n")))
rcat(capture.output(print(table(env$Clusters))))
rcat("```\n")

rcat("## Coverage correction model: NGSF ~ log10(Cov)\n")
rcat("```")
rcat(capture.output(print(summary(mod_cov))))
rcat("```\n")

rcat("## Ne-fog Pearson correlation (n=6 sites)\n")
rcat("```")
rcat(sprintf("r = %.3f, t = %.3f, df = %d, p = %.4f",
             Ne_fog_cor_inb$estimate, Ne_fog_cor_inb$statistic,
             Ne_fog_cor_inb$parameter, Ne_fog_cor_inb$p.value))
rcat("```\n")

rcat("## Site-level models: inbreeding residual ~ fog / Ne / both (n=6)\n")
rcat("```")
rcat(capture.output(print(summary(mod_inb_fog))))
rcat("")
rcat(capture.output(print(summary(mod_inb_Ne))))
rcat("")
rcat(capture.output(print(summary(mod_inb_both))))
rcat("")
rcat("AIC comparison:")
rcat(capture.output(print(AIC(mod_inb_fog, mod_inb_Ne, mod_inb_both))))
rcat("```\n")

rcat("## Variance partitioning: fog vs. Ne -> inbreeding residual\n")
rcat("```")
rcat(capture.output(print(vp_inb)))
rcat("```\n")

rcat("## Variance partitioning: BCM CWD vs. Ne -> inbreeding residual\n")
rcat("```")
rcat(capture.output(print(vp_inb_bcm)))
rcat("```\n")

rcat("## Ne-BIO18 (precip, warmest qtr.) Pearson correlation (n=6 sites)\n")
rcat("```")
rcat(sprintf("r = %.3f, t = %.3f, df = %d, p = %.4f",
             Ne_BIO18_cor_inb$estimate, Ne_BIO18_cor_inb$statistic,
             Ne_BIO18_cor_inb$parameter, Ne_BIO18_cor_inb$p.value))
rcat("```\n")

rcat("## Site-level models: inbreeding residual ~ BIO18 / Ne / both (n=6)\n")
rcat("```")
rcat(capture.output(print(summary(mod_inb_BIO18))))
rcat("")
rcat(capture.output(print(summary(mod_inb_BIO18_both))))
rcat("")
rcat("AIC comparison (BIO18 / Ne / BIO18+Ne):")
rcat(capture.output(print(AIC(mod_inb_BIO18, mod_inb_Ne, mod_inb_BIO18_both))))
rcat("```\n")

rcat("## Variance partitioning: BIO18 (precip, warmest qtr.) vs. Ne -> inbreeding residual\n")
rcat("```")
rcat(capture.output(print(vp_inb_BIO18)))
rcat("```\n")

rcat("## Site-level model comparison: fog + BIO18 + Ne, additive & interactive (n=6, AICc)\n")
rcat("")
rcat("Ne is a fixed-effect covariate here, not a random effect -- see comment in script")
rcat("for why a random effect isn't identifiable at n=6 sites with one Ne value each.")
rcat("The full 3-way interaction (8 parameters, 6 observations) is excluded as")
rcat("rank-deficient. Every 4-coefficient model below (the additive triple and all")
rcat("three pairwise interactions) returns AICc = Inf -- AICc's small-sample")
rcat("correction is undefined at that parameter count for n=6. The 3-coefficient")
rcat("additive pairs are finite but already 24+ AICc points worse than the best")
rcat("single predictor, so in practice only the three single-predictor models")
rcat("(fog, BIO18, Ne) are meaningfully comparable at this sample size.")
rcat("")
rcat("```")
rcat(capture.output(print(site_model_tab, row.names = FALSE)))
rcat("")
rcat(sprintf("Best model by AICc: residual ~ %s", best_site_model_name))
rcat(capture.output(print(summary(site_models[[best_site_model_name]]))))
rcat("```\n")

rcat("## Site-level variance partitioning: fog_germ_yr + BIO12_germ_yr + Ne -> inbreeding (n=6, 2020)\n")
rcat("```")
rcat(capture.output(print(as.data.frame(site_fogprecip))))
rcat("")
rcat(capture.output(print(vp_fogprecip)))
rcat("```\n")

rcat("## Strauss-only model comparison: BIO18 + fog + BIO15, additive & interactive (n=57, AICc)\n")
rcat("")
rcat(sprintf("Individual-level (Strauss only, >20x coverage, 2020 only), n = %d.",
             nrow(strauss10)))
rcat("2020-only, matching the site-level analyses above: a separate check (all 3 years,")
rcat("Year as a random intercept) found Year's variance was estimated at exactly 0 --")
rcat("singular fit even with no predictors at all -- so there's no detectable between-year")
rcat("variation to control for, and no need for the other years or a random effect here.")
rcat("Full candidate set fit (single terms, additive pairs/triple, pairwise interactions,")
rcat("all-pairwise, and the fully saturated model) since n=57 supports it.")
rcat("")
rcat("```")
rcat(capture.output(print(strauss_model_tab, row.names = FALSE)))
rcat("")
rcat(sprintf("Best model by AICc: residual ~ %s", best_strauss_model_name))
rcat(capture.output(print(summary(strauss_models[[best_strauss_model_name]]))))
rcat("```\n")

rcat("## Individual-level competing GLMs: residual ~ [env] + Clusters (all years)\n")
rcat(sprintf("n = %d samples (complete cases across all candidates), %d clusters, %d years\n",
             nrow(inb_all_env),
             length(unique(inb_all_env$Clusters)),
             length(unique(inb_all_env$Year.x))))
rcat("### AIC comparison table\n")
rcat("```")
rcat(capture.output(print(aic_tab, row.names = FALSE)))
rcat("```\n")
rcat("### Top 3 model summaries\n")
for (i in seq_len(3)) {
  v <- aic_tab$predictor[i]
  rcat(sprintf("#### Rank %d: residual ~ %s + Clusters  (AIC = %.2f, dAIC = %.2f)\n",
               i, v, aic_tab$AIC[i], aic_tab$delta_AIC[i]))
  rcat("```")
  rcat(capture.output(print(summary(glm_fits[[v]]))))
  rcat("```\n")
}

rcat("## Individual-level Beta GLM: NGSF ~ fog_germ_yr + Clusters (all years)\n")
rcat("")
rcat("### Model family and link function")
rcat("")
rcat("Family: **beta_family(link = \"logit\")** via glmmTMB.")
rcat("")
rcat("NGSF (the inbreeding coefficient) is a continuous proportion strictly bounded")
rcat("in (0, 1). The Gaussian assumption of ordinary lm is inappropriate here because")
rcat("residuals are not symmetric and variance is heteroscedastic near the boundaries.")
rcat("The beta distribution is the canonical choice for a continuous proportion response:")
rcat("it is flexible enough to capture both symmetric and skewed shapes over (0, 1).")
rcat("The logit link maps the unbounded linear predictor onto the (0, 1) interval,")
rcat("giving coefficients interpretable on a log-odds scale and ensuring predictions")
rcat("never exceed the biological bounds of the inbreeding coefficient.")
rcat("")
rcat("```")
rcat(capture.output(print(summary(mod_inb_fog_beta))))
rcat("")
rcat("Per-cluster fog_germ_yr slope p-values (lm on residual inbreeding):")
rcat(capture.output(print(round(fog_slope_pvals_inb, 4))))
rcat(sprintf("Significant clusters (p < 0.05): %s",
             if (length(sig_clusters_inb)) paste(sig_clusters_inb, collapse = ", ") else "none"))
rcat("```\n")

writeLines(report_lines, "inbreeding_analysis_report.md")
message("Report written to inbreeding_analysis_report.md")
