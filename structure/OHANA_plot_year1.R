#!/usr/bin/env Rscript

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(ggplot2)
library(tidyr)
library(scatterpie)
library(raster)
#library(ggnewscale)
library(sp)
library(dplyr)
#library(ggtree)
library(rgdal)
library(ggspatial)
library(maptools)
library(rgeos)
library(sf)
library(patchwork)
library(terra)


#Plotting params
gvtp_theme <- theme_bw() + 
  theme(
    text = element_text(size = 50), # Sets size for all text elements
    axis.text = element_text(angle = 45), # Only specify unique modifications
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm"),
    legend.box.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank())


proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


# add shapefile for 
# Define the cropping extent
crop_extent <- st_bbox(c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9), crs = st_crs(proj))
# crop the coast off
CA<-st_read("../PCAngsd/GVTP_clip.shp")
CA<-st_make_valid(CA)
crop_extent_sf <- st_as_sfc(st_bbox(crop_extent))
CA <- st_crop(CA, crop_extent)
CA<- st_simplify(CA)


# Read in the DEM
DEM <- rast("../PCAngsd/tranquillon_mountain_30m.dem")
crs(DEM) <- "EPSG:26710"
DEM_transformed <- project(DEM, proj)

# Calculate Evanno's K


# First, I need to compute the delta K� 


# Suppose your data has columns: iter, duration, log_likelihood, delta_lle, K
year1 <- read.table("year1/log_likelihoods_across_runs.txt", header=TRUE, sep="\t")


# Calculate mean and SD per K
summary <- year1 %>%
  group_by(K) %>%
  summarise(my_K = first(K),
    mean_lnL = mean(loglik),
    sd_lnL = sd(loglik))


# Compute ?K (Evanno method)
stats <- summary %>%
  mutate(lnL_lead = lead(mean_lnL),
    lnL_lag  = lag(mean_lnL),
    delta_K  = abs(lnL_lead - 2 * mean_lnL + lnL_lag) / sd_lnL)


# Plot ?K
jpeg("year1/deltaK_year1.jpg")
ggplot(stats, aes(x=K, y=delta_K)) +
  geom_line() +
  geom_point() +
  labs(x="K", y="Delta K (Evanno method)", title="Evanno Method for Choosing K") +
  theme_minimal()
dev.off()


best_k <- stats$K[which.max(stats$delta_K)]
cat("Best K according to Evanno's method:", best_k, "\n")
#     Best one 9 July is 4
#     Best one 31 Oct is 3 


label_vec <- c(
  CDFW = "Gaviota",
  Strauss = "Tranquillon Mountain",
  Point_Conception = "Conception",
  VAFB_Lion_Head = "Lion's Head",
  VAFB_Arguello = "Arguello",
  Dangermond = "D."
)


# Automatic ancestry-component coloring: for a given K, figure out which
# known population each component (C1..CK) is most tied to, then color it
# via popcolors. Strauss (by far the largest population) tends to split
# into multiple locally-dominant components as K increases; only the
# single largest of those keeps the Strauss popcolor (priority 1). Any
# additional Strauss-dominant components are colored by geography: the
# most-western one (lowest mean longitude) is treated as West Tranquillon
# Mtn -> black, and any others -> grey.
popcolors <- c(
  CDFW = "#e69f00",
  CDFW_MR = "#e69f00",
  Strauss = "#009e73",
  Point_Conception = "#0072b2",
  VAFB_Lion_Head = "#f0e442",
  VAFB_Arguello = "#d55e00",
  Dangermond = "#56b4e9"
)
priority_pops <- c("Strauss", "Point_Conception", "CDFW", "VAFB_Lion_Head", "VAFB_Arguello")

assign_ancestry_colors <- function(dat, k) {
  comp_cols <- paste0("C", seq_len(k))
  assignment <- apply(as.matrix(dat[, comp_cols]), 1, which.max)

  dominant_pop <- vapply(seq_len(k), function(c) {
    idx <- which(assignment == c)
    if (length(idx) == 0) return(NA_character_)
    tab <- table(as.character(dat$Clusters[idx]))
    names(tab)[which.max(tab)]
  }, character(1))

  colors <- setNames(rep(NA_character_, k), comp_cols)
  claimed <- rep(FALSE, k)

  # Each population, in priority order, claims its single best-supported
  # unclaimed component.
  for (pop in priority_pops) {
    candidates <- which(dominant_pop == pop & !claimed)
    if (length(candidates) == 0) next
    sizes <- vapply(candidates, function(c) sum(assignment == c), numeric(1))
    top <- candidates[which.max(sizes)]
    colors[top] <- popcolors[[pop]]
    claimed[top] <- TRUE
  }

  # Any further components still dominant within Strauss reflect a
  # within-Strauss split. The most-western one is West Tranquillon Mtn ->
  # black; any others -> grey.
  other_strauss <- which(dominant_pop == "Strauss" & !claimed)
  if (length(other_strauss) > 0) {
    mean_long <- vapply(other_strauss, function(c) mean(dat$Long[assignment == c], na.rm=TRUE), numeric(1))
    west_comp <- other_strauss[which.min(mean_long)]
    colors[west_comp] <- "black"
    claimed[west_comp] <- TRUE
    remaining <- setdiff(other_strauss, west_comp)
    if (length(remaining) > 0) {
      colors[remaining] <- "#999999"
      claimed[remaining] <- TRUE
    }
  }

  # Fallback for any component with no clear dominant population among the
  # above (e.g. a Dangermond-dominant component, or a genuinely ambiguous one).
  leftover <- which(!claimed)
  if (length(leftover) > 0) {
    fallback_palette <- c("#56b4e9", "white", "#cc79a7", "#a6761d")
    if (length(leftover) > length(fallback_palette)) {
      warning("assign_ancestry_colors: more unclaimed components than fallback colors; recycling")
    }
    colors[leftover] <- rep(fallback_palette, length.out = length(leftover))
  }

  colors
}


# Always read in the best-likelihood run (per best_runs.txt) for each K,
# rather than hard-coding which iteration was best.
best_runs <- read.table("year1/best_runs.txt", header=FALSE, sep="\t",
                         col.names=c("K", "iteration", "loglik"))

for (k in 2:8) {
  best_iter <- best_runs$iteration[best_runs$K == k]
  if (length(best_iter) != 1) {
    stop(sprintf("Expected exactly one best run for K=%d in best_runs.txt, found %d", k, length(best_iter)))
  }
  q_file <- sprintf("year1/q_angsd_paralog_LD_year1_%d_%d.matrix", k, best_iter)
  q <- read.table(q_file, skip=1, sep="\t", row.names=NULL)
  names(q) <- paste0("C", seq_len(k))
  assign(paste0("q_", k), q)
}


order<-read.table("/private/groups/meyerlab/esteiger/mapping/merged/bamlist_year1.txt", header=F)
order$index <- seq(1, nrow(order))
order$V1 <- sub(".*/([^/]+)$", "\\1", order$V1)
to_remove <- c("_rmdup.bam")
order$V1 <- gsub(paste(to_remove, collapse = "|"), "", order$V1)


#What's the population info for these samples?
pop<-read.table("/private/groups/meyerlab/esteiger/field_data.txt", header=T, sep="\t")
pop<-unique(pop)
pop<-merge(pop[,c("SampleID", "Year", "Lat", "Long", "Cov", "Clusters", "OHANA_order")], order, by.x="SampleID", by.y="V1", all.y=F)
my_pop<-pop[order(pop$index),]
my_pop$Clusters <- factor(my_pop$Clusters, levels=c("VAFB_Lion_Head", "VAFB_Arguello", "Dangermond", "Point_Conception", "Strauss", "CDFW"))


#merge pop data together with q matrices
data_2<-cbind(my_pop, q_2)
data_2<-data_2[!is.na(data_2$Clusters),]
# Sort by environ variable of interest and then re-index
data_2<-data_2[order(data_2$OHANA_order),]
data_2$index <- seq(1, nrow(data_2))
#hoizontalto vertical
data_long_2<-gather(data_2, Cluster, Percentage, C1:C2, factor_key=T)
data_3<-cbind(my_pop, q_3)
data_3<-data_3[!is.na(data_3$Clusters),]
data_3<-data_3[order(data_3$OHANA_order),]
data_3$index <- seq(1, nrow(data_3))
data_long_3<-gather(data_3, Cluster, Percentage, C1:C3, factor_key=T)
data_4<-cbind(my_pop, q_4)
data_4<-data_4[!is.na(data_4$Clusters),]
data_4<-data_4[order(data_4$OHANA_order),]
data_4$index <- seq(1, nrow(data_4))
data_long_4<-gather(data_4, Cluster, Percentage, C1:C4, factor_key=T)
data_5<-cbind(my_pop, q_5)
data_5<-data_5[!is.na(data_5$Clusters),]
data_5<-data_5[order(data_5$OHANA_order),]
data_5$index <- seq(1, nrow(data_5))
data_long_5<-gather(data_5, Cluster, Percentage, C1:C5, factor_key=T)
data_6<-cbind(my_pop, q_6)
data_6<-data_6[!is.na(data_6$Clusters),]
data_6<-data_6[order(data_6$OHANA_order),]
data_6$index <- seq(1, nrow(data_6))
data_long_6<-gather(data_6, Cluster, Percentage, C1:C6, factor_key=T)
data_7<-cbind(my_pop, q_7)
data_7<-data_7[!is.na(data_7$Clusters),]
data_7<-data_7[order(data_7$OHANA_order),]
data_7$index <- seq(1, nrow(data_7))
data_long_7<-gather(data_7, Cluster, Percentage, C1:C7, factor_key=T)
data_8<-cbind(my_pop, q_8)
data_8<-data_8[!is.na(data_8$Clusters),]
data_8<-data_8[order(data_8$OHANA_order),]
data_8$index <- seq(1, nrow(data_8))
data_long_8<-gather(data_8, Cluster, Percentage, C1:C8, factor_key=T)


# Evaluate how well the inferred clusters at each K recover the known
# field populations (Clusters), as a complement to Evanno's delta K.
# ARI and NMI are computed from scratch (base R) since neither mclust
# nor aricode is available in the project's R library.

adjusted_rand_index <- function(known, assignment) {
  tab <- table(known, assignment)
  n <- sum(tab)
  a <- rowSums(tab)
  b <- colSums(tab)
  sum_comb <- function(v) sum(v * (v - 1) / 2)
  index <- sum(tab * (tab - 1) / 2)
  expected_index <- sum_comb(a) * sum_comb(b) / (n * (n - 1) / 2)
  max_index <- (sum_comb(a) + sum_comb(b)) / 2
  if (max_index == expected_index) return(1)
  (index - expected_index) / (max_index - expected_index)
}

normalized_mutual_info <- function(known, assignment) {
  tab <- table(known, assignment)
  n <- sum(tab)
  p_xy <- tab / n
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)
  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      if (p_xy[i, j] > 0) {
        mi <- mi + p_xy[i, j] * log(p_xy[i, j] / (p_x[i] * p_y[j]))
      }
    }
  }
  h_x <- -sum(p_x[p_x > 0] * log(p_x[p_x > 0]))
  h_y <- -sum(p_y[p_y > 0] * log(p_y[p_y > 0]))
  if (h_x == 0 || h_y == 0) return(0)
  unname(mi / sqrt(h_x * h_y))
}

concordance <- data.frame(K = integer(), ARI = numeric(), NMI = numeric())
for (k in 2:8) {
  dat <- get(paste0("data_", k))
  comp_cols <- paste0("C", seq_len(k))
  assignment <- apply(as.matrix(dat[, comp_cols]), 1, which.max)
  known <- as.character(dat$Clusters)
  concordance <- rbind(concordance, data.frame(
    K = k,
    ARI = adjusted_rand_index(known, assignment),
    NMI = normalized_mutual_info(known, assignment)
  ))
}

P_deltaK <- ggplot(stats, aes(x=K, y=delta_K)) +
  geom_line() +
  geom_point() +
  labs(x="K", y="Delta K (Evanno method)", title="Evanno Method for Choosing K") +
  theme_minimal()

P_concordance <- concordance %>%
  gather(Metric, Value, ARI:NMI, factor_key=T) %>%
  ggplot(aes(x=K, y=Value, color=Metric)) +
  geom_line() +
  geom_point() +
  labs(x="K", y="Concordance with known populations",
       title="Cluster assignment vs. known populations") +
  theme_minimal()

jpeg("year1/ARI_NMI_year1.jpg")
print(P_concordance)
dev.off()

jpeg("year1/K_selection_year1.jpg", width=1600, height=800)
print(P_deltaK | P_concordance)
dev.off()


# By Clusters only:


P2<-ggplot(data=data_long_2, aes(x=as.character(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black") + 
facet_grid(cols=vars(Clusters), scales="free_x", space="free_x",  labeller = labeller(Clusters = label_vec)) +
  scale_fill_manual("", values=assign_ancestry_colors(data_2, 2))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
theme(strip.text=element_text(size=60, angle=45),
        axis.text = element_blank(), axis.title=element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"))  


P3<-ggplot(data=data_long_3, aes(x=as.character(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black") + 
facet_grid(cols=vars(Clusters), scales="free_x", space="free_x",  labeller = labeller(Clusters = label_vec)) +
  scale_fill_manual("", values=assign_ancestry_colors(data_3, 3))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
theme(strip.text=element_blank(),
        axis.text = element_blank(), axis.title=element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"))  


P4<-ggplot(data=data_long_4, aes(x=as.character(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black") + 
facet_grid(cols=vars(Clusters), scales="free_x", space="free_x",  labeller = labeller(Clusters = label_vec)) +
  scale_fill_manual("", values=assign_ancestry_colors(data_4, 4))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
theme(strip.text=element_blank(),
        axis.text = element_blank(), axis.title=element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"))  


P5<-ggplot(data=data_long_5, aes(x=as.character(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black") + 
facet_grid(cols=vars(Clusters), scales="free_x", space="free_x",  labeller = labeller(Clusters = label_vec)) +
  scale_fill_manual("", values=assign_ancestry_colors(data_5, 5))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
theme(strip.text=element_blank(),
        axis.text = element_blank(), axis.title=element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"))  


P6<-ggplot(data=data_long_6, aes(x=as.character(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black") + 
facet_grid(cols=vars(Clusters), scales="free_x", space="free_x",  labeller = labeller(Clusters = label_vec)) +
  scale_fill_manual("", values=assign_ancestry_colors(data_6, 6))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
theme(strip.text=element_blank(),
        axis.text = element_blank(), axis.title=element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"))  


P7<-ggplot(data=data_long_7, aes(x=as.character(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black") + 
facet_grid(cols=vars(Clusters), scales="free_x", space="free_x",  labeller = labeller(Clusters = label_vec)) +
  scale_fill_manual("", values=assign_ancestry_colors(data_7, 7))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
theme(strip.text=element_blank(),
        axis.text = element_blank(), axis.title=element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"))  


P8<-ggplot(data=data_long_8, aes(x=as.character(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black") +
facet_grid(cols=vars(Clusters), scales="free_x", space="free_x",  labeller = labeller(Clusters = label_vec)) +
  scale_fill_manual("", values=assign_ancestry_colors(data_8, 8))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
theme(strip.text=element_blank(),
        axis.text = element_blank(), axis.title=element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"))  






jpeg("year1/ohana_7.jpg", width=6400, height=1000)
print(P7 + theme(strip.text=element_text(size=75)))
dev.off()


# show the complete stack
jpeg("year1/ohana_stack.jpg", width=6000, height=4500)
print(P2 / P3 / P4 / P5 / P6 / P7 + plot_layout(heights = c(1, 1, 1, 1, 1,1)))
dev.off()








# Plot as pie charts on the landscape--INCOMPLETE
#average ancestry components by site
pie2 <- data_2 %>%
	dplyr::group_by(Clusters) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), Long=mean(Long), Lat=mean(Lat))


pie3 <- data_3 %>%
	dplyr::group_by(Clusters) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), Long=mean(Long), Lat=mean(Lat))


pie4 <- data_4 %>%
	dplyr::group_by(Clusters) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), Long=mean(Long), Lat=mean(Lat))


pie5 <- data_5 %>%
	dplyr::group_by(Clusters) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), Long=mean(Long), Lat=mean(Lat))


pie6 <- data_6 %>%
	dplyr::group_by(Clusters) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), F=mean(as.numeric(C6)), Long=mean(Long), Lat=mean(Lat))


pie7 <- data_7 %>%
	dplyr::group_by(Clusters) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), F=mean(as.numeric(C6)), G=mean(as.numeric(C7)), Long=mean(Long), Lat=mean(Lat))


pie8 <- data_8 %>%
	dplyr::group_by(Clusters) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), F=mean(as.numeric(C6)), G=mean(as.numeric(C7)), H=mean(as.numeric(C8)), Long=mean(Long), Lat=mean(Lat))




###################################


BASEPLOT<-ggplot2::ggplot() +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
xlim(-120.7, -120.1) +
ylim(34.4, 34.9) +
coord_sf(crs= 4326) + 
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 


jpeg("Map-pie2.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long, y=Lat, group=Clusters), data=pie2, cols=LETTERS[1:2], color="black", pie_scale=2)+
scale_fill_manual("", values=c(A="blue", B="red", C="yellow", D="green", E="orange", F="purple"))
dev.off()


jpeg("Map-pie3.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long, y=Lat, group=Clusters), data=pie3, cols=LETTERS[1:3], color="black", pie_scale=2)+
scale_fill_manual("", values=c(A="blue", B="red", C="yellow", D="green", E="orange", F="purple"))
dev.off()


jpeg("Map-pie4.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long, y=Lat, group=Clusters), data=pie4, cols=LETTERS[1:4], color="black", pie_scale=2)+
scale_fill_manual("", values=c(A="blue", B="red", C="yellow", D="green", E="orange", F="purple"))
dev.off()


jpeg("Map-pie5.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long, y=Lat, group=Clusters), data=pie5, cols=LETTERS[1:5], color="black", pie_scale=2)+
scale_fill_manual("", values=c(A="blue", B="red", C="yellow", D="green", E="orange", F="purple"))
dev.off()


jpeg("Map-pie6.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long, y=Lat, group=Clusters), data=pie6, cols=LETTERS[1:6], color="black", pie_scale=3)+
scale_fill_manual("", values=c(F="#d55e00", B="#999999", C="#009e73", D="#000000", E="#f0e442", A="#e69f00"))
dev.off()


jpeg("year1/Map-pie7.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long, y=Lat, group=Clusters), data=pie7, cols=LETTERS[1:7], color="black", pie_scale=2)+
scale_fill_manual("", values=c(C="#0072b2", E="#f0e442", F="#e69f00",  D="#000000", A="#009e73", B="#999999", G="#d55e00"))
dev.off()


jpeg("Map-pie8.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long, y=Lat, group=Clusters), data=pie8, cols=LETTERS[1:8], color="black", pie_scale=2)+
scale_fill_manual("", values=c(C="#d55e00", F="#0072b2", G="#f0e442", B="#56b4e9", A="black", D="white", E="#009e73", H="#e69f00"))
dev.off()


#popcolors <- c("CDFW" = "#e69f00", "CDFW_MR" = "#e69f00", "Strauss"="#009e73", "Point_Conception"="#0072b2", "VAFB_Lion_Head"="#f0e442", "VAFB_Arguello"="#d55e00", "Dangermond"="#56b4e9")         # Bang Wang color palette
# TANGERINE - "#e69f00" 
# GREEN - "#009e73"
# BLUE - "#0072b2"
# YELLOW - "#f0e442"
# ORANGE - "#d55e00"
# LIGHTBLUE - "#56b4e9"






#####################






data_7 <- data_7 %>%
  mutate(
    Long_jitter = jitter(Long, amount = 0.025),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = 0.025)
  )




pie7 <- data_7 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), F=mean(as.numeric(C6)), G=mean(as.numeric(C7)))
pie7<-merge(pie7, data_7[,c(1:7,16,17)], by.x="SampleID", by.y="SampleID")
pie7<-unique(pie7)


jpeg("year1/Map-pie7-jittered.jpg", width=1500, height=1500)
BASEPLOT +
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie7, cols=LETTERS[1:7], color="black", pie_scale=.5)+
scale_fill_manual("", values=c(C="#0072b2", E="#f0e442", F="#e69f00",  D="#000000", A="#009e73", B="#999999", G="#d55e00")) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()



