setwd("~/Library/CloudStorage/GoogleDrive-sscheim@ucsc.edu/.shortcut-targets-by-id/1th4k9MOXb4NFZdnAlKncS7yQaC1BQnAT/GaviotaTarplant/Analysis")
gt<-read.table("Data/phenotype/gvtp_phenotype_6nov2025.tsv", sep='\t', header=T)
library(tidyverse)
library(ggpubr)
library(rstatix)

subset<-gt[,c("Count_ray_flourets", "iNat_ray_count", "recordNumber_CollectorNumber")]
subset <- subset %>% drop_na(Count_ray_flourets)
subset <- subset %>% drop_na(iNat_ray_count)
subset.long <- subset %>%
  gather(key = "group", value = "weight", Count_ray_flourets, iNat_ray_count)
head(subset.long, 3)
subset <- subset %>% mutate(differences = Count_ray_flourets - iNat_ray_count)
head(subset, 3)
subset %>% identify_outliers(differences)
subset %>% shapiro_test(differences) 
#not normal
stat.test <- subset.long %>% 
  wilcox_test(weight ~ group, paired=T) %>%
  add_significance()
stat.test

subset <- subset(gt, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality, BIO12, BIO15, 
                                BIO18, fog, CHILI, recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd, 
                                diff, burned, number, cause))
modern <- subset(subset, year>2000)
modern2 <- modern %>% drop_na(Count_ray_flourets)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- as.factor(modern2$Locality)
modern2$Locality <- droplevels(modern2$Locality)

krusk<- kruskal.test(Count_ray_flourets ~ Locality, data = modern2) 
print(krusk)
#ns
modern2 <- modern2 %>% drop_na(Count_disk_flourets)
krusk<- kruskal.test(Count_disk_flourets ~ Locality, data = modern2) 
print(krusk)
#sig
library(FSA)
dunnTest(modern2$Count_disk_flourets, modern2$Locality, method='bonferroni')

modern2$total_florets <- with(modern2, Count_ray_flourets + Count_disk_flourets)
modern2$ray_prop <- with(modern2, Count_ray_flourets/total_florets)
krusk<- kruskal.test(ray_prop ~ Locality, data = modern2) 
print(krusk)


lm <- lm(Count_ray_flourets ~ year, data = gt)
#ns
lmdata <- gt %>% drop_na(Count_disk_flourets)
lm <- lm(Count_disk_flourets ~ year, data = lmdata)
#sig

summary(lm)

