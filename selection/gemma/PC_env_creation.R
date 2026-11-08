#!/usr/bin/env Rscript
# PCA of the PCAngsd covariance matrix (>4x samples). Writes GEMMA inputs env.txt and
# covariates_intercept.txt (intercept, PC1-5, coverage), plus scree, collinearity and PCA plots.
# R library: BASE_DIR/software/R_libs/

BASE_DIR <- "/private/groups/meyerlab/esteiger"

.libPaths(file.path(BASE_DIR, "software/R_libs/"))
library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)
library(psych)

gvtp_theme <- theme_bw() + 
  theme(text = element_text(size = 50),
    axis.text = element_text(angle = 45, hjust=1),
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

custom_labels <- c("2020" = "Year 1", "2021" = "Year 2",  "2022" = "Year 3")

clustercolors<-c("76"="#99DDFF", #sky
"77"="#EEDD88", #yellow
"79"="#77AADD", #blue
"73"="#AAAA00", #green
"74"="#AAAA00", #green
"75"="#AAAA00", #green
"71"="#EE8866",  #orange
"72"= "#DDDDDD", #grey
"78"= "#99DDFF",
"80"=  "#77AADD", #blue
"81"=  "#77AADD", #blue
"82"= "#77AADD", #blue
"83"=  "#77AADD", #blue,
"84"= "#77AADD", #blue 
"85"= "#77AADD", #blue
"86"= "#77AADD", #blue
"87"=  "#77AADD", #blue
"88"= "#77AADD")

location_labels <- c(
"CDFW"="Gaviota",             
"Strauss"="Tranquillon Mtn",
"Dangermond"="Dangermond",  
 "Point_Conception"="Pt, Conception",
"VAFB_Lion_Head"="Lion's Head",
"VAFB_Arguello"="Pt. Arguello"
)

#What's the population info for these samples?
pop<-read.table(file.path(BASE_DIR, "field_data.txt"), header=T, sep="\t")
pop<-pop[order(pop$index),]
pop$Location <- as.factor(pop$Location)
pop$Year <-as.factor(pop$Year)
pop$Cov<-as.numeric(pop$Cov)

# Sample order of the bamlist (matches the PCAngsd covariance matrix)
sample_order<-read.table(file.path(BASE_DIR, "mapping/merged/bamlist_4x.txt"), header=F, sep="\t")
sample_order$index<-seq(1, nrow(sample_order))
to_remove <- c("_rmdup.bam")
sample_order$V1 <- gsub(paste(to_remove, collapse = "|"), "", sample_order$V1)
my_pop<-merge(pop, sample_order, by.y="V1", by.x="SampleID", all=F)

# Read in env variables 
env<-read.table(file.path(BASE_DIR, "env_var.txt"), header=T, sep="\t")
colnames(env)[colnames(env)=="PercVolSoilH20_10cm"] <- "PercH20"
my_pop<-merge(my_pop, env[, !names(env) %in% c("Year","Cov")], by="SampleID")
# Read in new env variables from Cal-eDNA
env<-read.table(file.path(BASE_DIR, "eDNA_env_data.csv"), header=T, sep=",")
colnames(env)[colnames(env) == "terraclimate_palmer_drought_severity_index_30y"] <- "severity"
colnames(env)[colnames(env) == "terraclimate_severe_drought_month_count_30y"] <- "drought_mn"
colnames(env)[colnames(env) == "terraclimate_soil_moisture_mm_30y"] <- "moisture"
colnames(env)[colnames(env) == "terraclimate_climate_water_deficit_mm_30y"] <- "deficit"
env <- env[, c("SampleID", "severity", "drought_mn", "moisture", "deficit")]
my_pop<-merge(my_pop, env, all.x=T, by="SampleID")
my_pop[my_pop == -999999] <- NA

# eigen decomposition of covariance matrix
C <- as.matrix(read.table(file.path(BASE_DIR, "structure/PCAngsd/postParalog_4x_LD.cov")))
e <- eigen(C)

# Scree plot, to decide how many PCs to retain
scree <- data.frame(
  PC = seq_along(e$values),
  Eigenvalue = e$values
)

jpeg(file.path(BASE_DIR, "selection/gemma/scree.jpg"), height=1000, width=1000)
ggplot(scree, aes(x = PC, y = Eigenvalue)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(breaks = scree$PC) +
  labs(x = "Principal Component", y = "Eigenvalue", title = "Scree Plot") + theme_bw()
dev.off()

PC1<- round((e$values/sum(e$values))[1], 2)*100
PC2<- round((e$values/sum(e$values))[2], 2)*100
PC3<- round((e$values/sum(e$values))[3], 2)*100
PC4<- round((e$values/sum(e$values))[4], 2)*100
PC5<- round((e$values/sum(e$values))[5], 2)*100

myPCs <- e$vectors[,1:5]
myPCs <- as.data.frame(myPCs)
names(myPCs)<- c("PC1", "PC2", "PC3", "PC4", "PC5")

# Bind PCs to population data before filtering
my_pop <- my_pop[order(my_pop$index.y),]
my_pop <- cbind(my_pop, myPCs)
stopifnot(all(diff(my_pop$index.y) > 0))

# Remove 13 individuals with untrustworthy coordinates, then export env variables and PCs
EXCLUDE<-c("4548", "4549", "4550", "4551", "4552", "4553", "4554", "4555", "4556", "4557", "4559", "4561", "4562")
my_pop_filtered <- my_pop[!(my_pop$SampleID %in% EXCLUDE),]
write.table(my_pop_filtered[,c("PercH20",  "CHILI",  "BIO12", "BIO15", "BIO18", "fog", "fog_germ_yr", "BIO12_germ_yr", "BIO15_germ_yr", "BIO18_germ_yr", "BCM_30y", "BCM_germ_yr", "deficit", "severity", "drought_mn", "moisture")], file.path(BASE_DIR, "selection/gemma/env.txt"), row.names=F, col.names=F)
# GEMMA does not add an intercept when -c is given, so include a column of 1s
covariates <- cbind(Intercept = 1, my_pop_filtered[,c("PC1", "PC2", "PC3", "PC4", "PC5", "Cov")])
write.table(covariates, file.path(BASE_DIR, "selection/gemma/covariates_intercept.txt"), row.names=F, col.names=F)

# Collinearity of environmental variables
jpeg(file.path(BASE_DIR, "selection/gemma/collin.jpg"), height=1000, width=1000)
pairs.panels(my_pop[,c("PercH20",  "CHILI",  "BIO12", "BIO15", "BIO18", "fog", "fog_germ_yr", "BIO12_germ_yr", "BIO15_germ_yr", "BIO18_germ_yr", "BCM_30y", "BCM_germ_yr", "deficit", "severity", "drought_mn", "moisture")], scale=T)
dev.off()

# Colored by locality
P1<-ggplot(data=my_pop, aes(x=PC1, y=PC2, fill=Clusters)) +
  geom_point(colour="black", size=9, shape=21)+
  scale_fill_manual(values=popcolors, guide="none")+
  labs(x=paste("PC1 (", PC1,"% variance)"), y=paste("PC2 (",PC2,"% variance)"), fill="Locality", shape="Age")+
 gvtp_theme + coord_cartesian(expand = TRUE)
P3<-ggplot(my_pop, aes(PC3, PC4, fill=Clusters)) +
  geom_point(colour="black", size=9, shape=21)+
  scale_fill_manual(values=popcolors, guide="none")+
  labs(x=paste("PC3 (", PC3,"% variance)"), y=paste("PC4 (",PC4,"% variance)"), fill="Region", shape="Age")+
 gvtp_theme + coord_cartesian(expand = TRUE)

jpeg(file.path(BASE_DIR, "selection/gemma/PCA_12.jpg"), height=1000, width=1000)
print(P1)
dev.off()
jpeg(file.path(BASE_DIR, "selection/gemma/PCA_34.jpg"), height=1000, width=1000)
print(P3)
dev.off()
