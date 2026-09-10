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


# First, I need to compute the delta K… 
#LOCALITY<-"strauss"
#LOCALITY<-"VAFB_Lion_Head"
# LOCALITY<-"CDFW"
#LOCALITY<-"other"


# Suppose your data has columns: iter, duration, log_likelihood, delta_lle, K
year1 <- read.table(paste0(LOCALITY, "/log_likelihoods_across_runs.txt"), header=TRUE, sep="\t")


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
jpeg(paste0(LOCALITY,"/deltaK_year1.jpg"))
ggplot(stats, aes(x=K, y=delta_K)) +
  geom_line() +
  geom_point() +
  labs(x="K", y="Delta K (Evanno method)", title="Evanno Method for Choosing K") +
  theme_minimal()
dev.off()


best_k <- stats$K[which.max(stats$delta_K)]
cat("Best K according to Evanno's method:", best_k, "\n")
#     strauss : 5
#     Lion's Head: 3
#     CDFW: 5
#     other:3


# Find the best iteration (max LogLik) for each K
best_iters <- year1 %>%
  group_by(K) %>%
  slice_max(loglik, n = 1, with_ties = FALSE) %>%  # Don't keep ties
  dplyr::select(K, best_iter = iteration, best_LogLik = loglik)


# Create a function to read Q matrix for best iteration
read_best_Q <- function(k_val, best_iters_df, dir) {
  best_iter <- best_iters_df$best_iter[best_iters_df$K == k_val]
  filename <- paste0(LOCALITY, "/q_angsd_paralog_LD_", LOCALITY, "_", k_val, "_", best_iter, ".matrix") 
  if (!file.exists(filename)) {
    warning(paste("File not found:", filename))
    return(NULL)
  }  
  q_mat <- read.table(filename, skip=1, sep="\t", row.names=NULL)
  col_names <- paste0("C", 1:k_val)
  names(q_mat) <- col_names  
  cat("Read K=", k_val, " from iteration ", best_iter, "\n")
  return(q_mat)
}


# Read best Q matrices for K=2 to K=6
q_2 <- read_best_Q(2, best_iters, LOCALITY)
q_3 <- read_best_Q(3, best_iters, LOCALITY)
q_4 <- read_best_Q(4, best_iters, LOCALITY)
q_5 <- read_best_Q(5, best_iters, LOCALITY)
q_6 <- read_best_Q(6, best_iters, LOCALITY)


my_order<-read.table(paste0("/private/groups/meyerlab/esteiger/mapping/merged/bamlist_", LOCALITY, ".txt"), header=F)    # There are 736 bams
my_order$index <- seq(1, nrow(my_order))
my_order$V1 <- sub(".*/([^/]+)$", "\\1", my_order$V1)
to_remove <- c("_rmdup.bam")
my_order$V1 <- gsub(paste(to_remove, collapse = "|"), "", my_order$V1)


#What's the population info for these samples?
pop<-read.table("/private/groups/meyerlab/esteiger/field_data.txt", header=T, sep="\t")
pop<-unique(pop)
pop<-merge(pop[,c("SampleID", "Year", "Lat", "Long", "Cov", "OHANA_order")], my_order, by.x="SampleID", by.y="V1", all.y=F)
my_pop<-pop[order(pop$index),]


#merge pop data together with q matrices
data_2<-cbind(my_pop, q_2)
# Sort by environ variable of interest and then re-index
data_2<-data_2[order(data_2$OHANA_order, data_2$Long), ]
data_2$index <- seq(1, nrow(data_2))
#hoizontalto vertical
data_long_2<-gather(data_2, Cluster, Percentage, C1:C2, factor_key=T)


data_3<-cbind(my_pop, q_3)
data_3<-data_3[order(data_3$OHANA_order, data_3$Long), ]
data_3$index <- seq(1, nrow(data_3))
data_long_3<-gather(data_3, Cluster, Percentage, C1:C3, factor_key=T)


data_4<-cbind(my_pop, q_4)
data_4<-data_4[order(data_4$OHANA_order, data_4$Long), ]
data_4$index <- seq(1, nrow(data_4))
data_long_4<-gather(data_4, Cluster, Percentage, C1:C4, factor_key=T)


data_5<-cbind(my_pop, q_5)
data_5<-data_5[order(data_5$OHANA_order, data_5$Long), ]
data_5$index <- seq(1, nrow(data_5))
data_long_5<-gather(data_5, Cluster, Percentage, C1:C5, factor_key=T)


data_6<-cbind(my_pop, q_6)
data_6<-data_6[order(data_6$OHANA_order, data_6$Long), ]
data_6$index <- seq(1, nrow(data_6))
data_long_6<-gather(data_6, Cluster, Percentage, C1:C6, factor_key=T)


data_7<-cbind(my_pop, q_7)   # differing number of rows: 736, 749
data_7<-data_7[order(data_7$OHANA_order, data_7$Long), ]
data_7$index <- seq(1, nrow(data_7))
data_long_7<-gather(data_7, Cluster, Percentage, C1:C7, factor_key=T)


data_8<-cbind(my_pop, q_8)
data_8<-data_8[order(data_8$OHANA_order, data_8$Long), ]
data_8$index <- seq(1, nrow(data_8))
data_long_8<-gather(data_8, Cluster, Percentage, C1:C8, factor_key=T)




# Remove individuals whose coordinates I do not trust. 
EXCLUDE<-c("4548", "4549", "4550", "4551", "4552", "4553", "4554", "4555", "4556", "4557", "4559", "4561", "4562")
data_2<-data_2[!(data_2$SampleID %in% EXCLUDE),]
data_3<-data_3[!(data_3$SampleID %in% EXCLUDE),]
data_4<-data_4[!(data_4$SampleID %in% EXCLUDE),]
data_5<-data_5[!(data_5$SampleID %in% EXCLUDE),]
data_6<-data_6[!(data_6$SampleID %in% EXCLUDE),]
data_7<-data_7[!(data_7$SampleID %in% EXCLUDE),]
data_8<-data_8[!(data_8$SampleID %in% EXCLUDE),]


# ONLY SORTED BY CLUMPS OF PLANTS


P2<-ggplot(data=data_long_2, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black", width=1) + 
#facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
scale_fill_manual("", values=c("C1"="#EEDD88", "C2"="#77AADD"))+
scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")  


P3<-ggplot(data=data_long_3, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black", width=1) + 
#facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C2"="#EEDD88", "C3"="#77AADD", "C1"="#AAAA00"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
#scale_x_continuous(expand = c(0,0))+
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


P4<-ggplot(data=data_long_4, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black", width=1) + 
#facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C4"="#77AADD", "C3"="#EEDD88", "C2"="#EE8866", "C1"="#AAAA00"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


P5<-ggplot(data=data_long_5, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black", width=1) + 
#facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C3"="#EE8866", "C2"="#EEDD88", "C4"="#77AADD",  "C1"="#99DDFF", "C5"="#AAAA00"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none") 


# yellow=#EEDD88
# blue=#77AADD
# green=#AAAA00
# orange=#EE8866
# sky=#99DDFF
# grey=#DDDDDD
# turquoise="#44BB99"


P6<-ggplot(data=data_long_6, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black", width=1) + 
#facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C5"="#AAAA00",  "C6"="#77AADD", "C2"="#DDDDDD", "C4"="#99DDFF", "C3"="#EE8866", "C1"="#EEDD88" ))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


P7<-ggplot(data=data_long_7, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black", width=1) + 
#facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C4"="#EEDD88", "C5"="#77AADD", "C3"="#99DDFF", "C6"="#AAAA00", "C2"="#44BB99", "C1"="#EE8866", "C7"="#DDDDDD"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


P8<-ggplot(data=data_long_8, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black", width=1) + 
#facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C1"="white", "C6"="#BBCC33", "C3"="#AAAA00", "C4"="#99DDFF", "C5"="#44BB99", "C2"="#DDDDDD", "C7"="#EEDD88", "C8"="#77AADD"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


# show the complete stack
jpeg(paste0(LOCALITY, "/", LOCALITY, "_stack.jpg"), width=4000, height=5000)
print(P2 / P3 / P4 / P5 / P6  + plot_layout(heights = c(1, 1, 1, 1, 1)))
dev.off()


jpeg("year1/ohana_7.jpg", width=6400, height=1000)
print(P7 + theme(strip.text=element_text(size=80, angle=40)))
dev.off()




#  SPLIT BY GERMINATION YEAR


process_data <- function(q_data, n_clusters) {
  cluster_cols <- paste0("C", 1:n_clusters)
  
  result <- cbind(my_pop, q_data) %>%
    arrange(OHANA_order) %>%
    group_by(Year) %>%
    mutate(index = row_number()) %>%
    ungroup()
  
  result_long <- gather(result, Cluster, Percentage, all_of(cluster_cols), factor_key = TRUE)
  
  return(list(wide = result, long = result_long))
}


# Use the function
output_2 <- process_data(q_2, 2)
data_2 <- output_2$wide
data_long_2 <- output_2$long


output_3 <- process_data(q_3, 3)
data_3 <- output_3$wide
data_long_3 <- output_3$long


output_4 <- process_data(q_4, 4)
data_4 <- output_4$wide
data_long_4 <- output_4$long


output_5 <- process_data(q_5, 5)
data_5 <- output_5$wide
data_long_5 <- output_5$long








P2<-ggplot(data=data_long_2, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black", width=1) + 
facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
scale_fill_manual("", values=c("C1"="#EEDD88", "C2"="#77AADD"))+
scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")  


P3<-ggplot(data=data_long_3, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black", width=1) + 
facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C2"="#EEDD88", "C3"="#77AADD", "C1"="#AAAA00"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
gvtp_theme +
#scale_x_continuous(expand = c(0,0))+
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


P4<-ggplot(data=data_long_4, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black", width=1) + 
facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C4"="#77AADD", "C3"="#EEDD88", "C2"="#EE8866", "C1"="#AAAA00"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


P5<-ggplot(data=data_long_5, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
geom_bar(stat="identity", color="black", width=1) + 
facet_grid(rows=vars(as.factor(Year)), scales="free_x", space="free_x")+
  scale_fill_manual("", values=c("C3"="#EE8866", "C2"="#EEDD88", "C4"="#77AADD",  "C1"="#99DDFF", "C5"="#AAAA00"))+
  scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2))+
#scale_x_continuous(expand = c(0,0))+
gvtp_theme +
theme(strip.text=element_text(size=40),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none") 






# show the complete stack
jpeg("strauss_all/ohana_stack_year.jpg", width=4000, height=4000)
print(P2 / P3 / P4 / P5  + plot_layout(heights = c(1, 1, 1, 1)))
dev.off()










# NOW PIE CHART MAPS


# Calculate extent from scatterpie locations with some padding
x_range <- range(pie2$Long, na.rm = TRUE)
y_range <- range(pie2$Lat, na.rm = TRUE)
# Add padding (adjust the multiplier as needed)
x_padding <- diff(x_range) * 0.2  # 10% padding
y_padding <- diff(y_range) * 0.2
# Calculate plot extent
x_extent <- diff(x_range) + 2 * x_padding
y_extent <- diff(y_range) + 2 * y_padding
# Calculate aspect ratio
aspect_ratio <- x_extent / y_extent
# Set base dimension and calculate the other
base_height <- 1500
plot_width <- base_height * aspect_ratio
plot_height <- base_height
base_scale <- 0.20 
jitter<-x_extent/1000


#   create some jittered coordinates
data_2 <- data_2 %>%
  mutate(
    Long_jitter = jitter(Long, amount = jitter),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = jitter)
  )
data_3 <- data_3 %>%
  mutate(
    Long_jitter = jitter(Long, amount = jitter),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = jitter)
  )
data_4 <- data_4 %>%
  mutate(
    Long_jitter = jitter(Long, amount = jitter),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = jitter)
  )
data_5 <- data_5 %>%
  mutate(
    Long_jitter = jitter(Long, amount = jitter),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = jitter)
  )
data_6 <- data_6 %>%
  mutate(
    Long_jitter = jitter(Long, amount = jitter),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = jitter)
  )
data_7 <- data_7 %>%
  mutate(
    Long_jitter = jitter(Long, amount = jitter),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = jitter)
  )
data_8 <- data_8 %>%
  mutate(
    Long_jitter = jitter(Long, amount = jitter),  # adjust amount as needed
    Lat_jitter  = jitter(Lat, amount = jitter)
  )


# Plot as pie charts on the landscape--INCOMPLETE
#average ancestry components by site
pie2 <- data_2 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)))
pie2<-merge(pie2, data_2[,c(1:7,10,11)], by.x="SampleID", by.y="SampleID")
pie2<-unique(pie2)


pie3 <- data_3 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)))
pie3<-merge(pie3, data_3[,c(1:7,11,12)], by.x="SampleID", by.y="SampleID")
pie3<-unique(pie3)


pie4 <- data_4 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)))
pie4<-merge(pie4, data_4[,c(1:7,12,13)], by.x="SampleID", by.y="SampleID")
pie4<-unique(pie4)


pie5 <- data_5 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)))
pie5<-merge(pie5, data_5[,c(1:7,13,14)], by.x="SampleID", by.y="SampleID")
pie5<-unique(pie5)


pie6 <- data_6 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), F=mean(as.numeric(C6)))
pie6<-merge(pie6, data_6[,c(1:7,14,15)], by.x="SampleID", by.y="SampleID")
pie6<-unique(pie6)


pie7 <- data_7 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), F=mean(as.numeric(C6)), G=mean(as.numeric(C7)))
pie7<-merge(pie7, data_7[,c(1:7,15,16)], by.x="SampleID", by.y="SampleID")
pie7<-unique(pie7)


pie8 <- data_8 %>%
	dplyr::group_by(SampleID) %>%
dplyr::summarise(A=mean(as.numeric(C1)), B=mean(as.numeric(C2)), C=mean(as.numeric(C3)), D=mean(as.numeric(C4)), E=mean(as.numeric(C5)), F=mean(as.numeric(C6)), G=mean(as.numeric(C7)), H=mean(as.numeric(C8)))
pie8<-merge(pie8, data_8[,c(1:7,16,17)], by.x="SampleID", by.y="SampleID")
pie8<-unique(pie8)




# crop to strauss extent
crop_extent <- st_bbox(c(xmin =-120.55, xmax = -120.50, ymin = 34.56, ymax = 34.605), crs = proj)
crop_extent_terra <- ext(crop_extent[1], crop_extent[3], crop_extent[2], crop_extent[4])
DEM_cropped <- crop(DEM_transformed, crop_extent_terra)
DEM_cropped <- DEM_transformed
contour_lines <- as.contour(DEM_cropped, levels = seq(0, 800, 25))  # Adjust levels as needed
contour_sf <- st_as_sf(contour_lines)


jpeg(paste0(LOCALITY, "/", LOCALITY, "_Map-pie2.jpg"), width=plot_width, height=plot_height)
ggplot2::ggplot() +
 geom_sf(data=contour_sf, color="gray70", size=1, inherit.aes=FALSE, linewidth=2) +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie2, cols=LETTERS[1:2], color="black", pie_scale=base_scale * (2 / 3) / x_extent)+
scale_fill_manual("", values=c("A"="#EEDD88", "B"="#77AADD"), guide="none")+
 coord_sf(crs = 4326, 
           xlim = c(x_range[1] - x_padding, x_range[2] + x_padding), 
           ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()


jpeg(paste0(LOCALITY, "/", LOCALITY, "_Map-pie3.jpg"), width=plot_width, height=plot_height)
ggplot2::ggplot() +
 geom_sf(data=contour_sf, color="gray70", size=0.3, inherit.aes=FALSE, linewidth=2) +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie3, cols=LETTERS[1:3], color="black", pie_scale=base_scale * (3/ 3) / x_extent)+
  scale_fill_manual("", values=c("A"="#EEDD88", "B"="#77AADD", "C"="#AAAA00"), guide="none")+
 coord_sf(crs = 4326, 
           xlim = c(x_range[1] - x_padding, x_range[2] + x_padding), 
           ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()


jpeg(paste0(LOCALITY, "/", LOCALITY, "_Map-pie4.jpg"), width=plot_width, height=plot_height)
ggplot2::ggplot() +
 geom_sf(data=contour_sf, color="gray70", size=0.3, inherit.aes=FALSE, linewidth=2) +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie4, cols=LETTERS[1:4], color="black", pie_scale=base_scale * (4/ 3) / x_extent)+
  scale_fill_manual("", values=c("A"="#77AADD", "B"="#EEDD88", "C"="#AAAA00", "D"="#EE8866"), guide="none")+
 coord_sf(crs = 4326, 
           xlim = c(x_range[1] - x_padding, x_range[2] + x_padding), 
           ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()


jpeg(paste0(LOCALITY, "/", LOCALITY, "_Map-pie5.jpg"), width=plot_width, height=plot_height)
ggplot2::ggplot() +
 geom_sf(data=contour_sf, color="gray70", size=0.3, inherit.aes=FALSE, linewidth=2) +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie5, cols=LETTERS[1:5], color="black", pie_scale=base_scale * (5 / 3) / x_extent)+
  scale_fill_manual("", values=c("A"="#EE8866", "B"="#EEDD88", "C"="#99DDFF", "D"="#77AADD", "E"="#AAAA00"), guide="none")+
 coord_sf(crs = 4326, 
           xlim = c(x_range[1] - x_padding, x_range[2] + x_padding), 
           ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()


jpeg(paste0(LOCALITY, "/", LOCALITY, "_Map-pie6.jpg"), width=plot_width, height=plot_height)
ggplot2::ggplot() +
 geom_sf(data=contour_sf, color="gray70", size=0.3, inherit.aes=FALSE, linewidth=2) +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie6, cols=LETTERS[1:6], color="black", pie_scale=base_scale * (6 / 3) / x_extent)+
  scale_fill_manual("", values=c("E"="#AAAA00", "D"="#77AADD", "C"="#EE8866", "A"="#EEDD88", "B"="#99DDFF", "F"="#DDDDDD"), guide="none")+
 coord_sf(crs = 4326, 
           xlim = c(x_range[1] - x_padding, x_range[2] + x_padding), 
           ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()


jpeg(paste0(LOCALITY, "/", LOCALITY, "_Map-pie7.jpg"), width=plot_width, height=plot_height)
ggplot2::ggplot() +
 geom_sf(data=contour_sf, color="gray70", size=0.3, inherit.aes=FALSE, linewidth=2) +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie7, cols=LETTERS[1:7], color="black", pie_scale=base_scale * (7 / 3) / x_extent)+
  scale_fill_manual("", values=c(D="#EEDD88", E="#77AADD", C="#99DDFF", F="#AAAA00", B="#44BB99", A="#EE8866", G="#DDDDDD"), guide="none")+
 coord_sf(crs = 4326, 
           xlim = c(x_range[1] - x_padding, x_range[2] + x_padding), 
           ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()


jpeg(paste0(LOCALITY, "/", LOCALITY, "_Map-pie8.jpg"), width=plot_width, height=plot_height)
ggplot2::ggplot() +
geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) + 
geom_scatterpie(aes(x=Long_jitter, y=Lat_jitter, group=SampleID), data=pie8, cols=LETTERS[1:8], color="black", pie_scale=base_scale * (8 / 3) / x_extent)+
scale_fill_manual("", values=c(A="white", B="grey", C="yellow", D="orange", E="#44BB99", F="#DDDDDD", G="#EEDD88", H="#AAAA00"))+
 coord_sf(crs = 4326, 
           xlim = c(x_range[1] - x_padding, x_range[2] + x_padding), 
           ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
theme(panel.background = element_blank()) +
theme(axis.title = element_blank(),
    	axis.ticks = element_blank(),
    	axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
    	panel.grid.minor=element_blank(),
    	plot.margin=grid::unit(c(0,0,0,0), "mm")) 
dev.off()
