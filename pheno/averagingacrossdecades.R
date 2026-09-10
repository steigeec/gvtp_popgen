setwd("~/Library/CloudStorage/GoogleDrive-sscheim@ucsc.edu/.shortcut-targets-by-id/1th4k9MOXb4NFZdnAlKncS7yQaC1BQnAT/GaviotaTarplant/Analysis/Data/phenotype/climateNA")
library(psych)
library(tidyverse)
library(dplyr)

#read all files in as one dataframe

tbl <-
  list.files(pattern = "\\.csv$") %>% 
  map_df(~read_csv(.))

#filter to only variables of interest - MAT MAP MSP and CMI
tbl <- tbl[ , c("ID1", "Latitude", "Longitude", "MAT", "MAP",
                        "MSP", "CMI")]

#calculate mean -> group by ID, mean for each column
means <- tbl %>% 
  group_by(ID1) %>% 
  summarize(across(everything(), list(mean)))

#stuff to see collinearities - decided on variables
subset<-means[,c("MAT_1", "CMI_1", "MSP_1", "MAP_1")] # drop MAP or MSP?
par(mar = c(1, 1, 1, 1))
pairs.panels(subset, scale=T) 
dev.off()

setwd("~/Library/CloudStorage/GoogleDrive-sscheim@ucsc.edu/.shortcut-targets-by-id/1th4k9MOXb4NFZdnAlKncS7yQaC1BQnAT/GaviotaTarplant/Analysis")

#join gt and means
#rename ID1
colnames(means)[colnames(means) == 'ID1'] <- 'recordNumber_CollectorNumber'

merged <- merge(gt, means, by="recordNumber_CollectorNumber", all=T)

#MODELS

filter <- merged %>%
  filter(recordNumber_CollectorNumber <=8484 | recordNumber_CollectorNumber>=8501 )
filter <- filter %>% drop_na(Count_disk_flourets)


subset <- subset(filter, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality,
                                    recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd, 
                                    diff, burned, number, cause))

subset$Locality <- as.factor(subset$Locality)
subset <- filter[subset$Locality != "Dangermond",]
vars3 <- c("MAT_1", "MAP_1", "MSP_1" ,"CMI_1") #
models3 <- list()

#rays
# All additive combinations
for(i in 1:length(vars3)) {
  combos <- combn(vars3, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " + "))
    models3[[length(models3) + 1]] <- glm(as.formula(formula_str), data = subset, family=poisson)
  }
}

# Add interaction models -- 22nd
for(i in 2:length(vars3)) {
  i<-4
  combos <- combn(vars3, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " * "))
    models3[[length(models3) + 1]] <- glm(as.formula(formula_str), data = subset, family=poisson)
  }
}
aictab(models3)
summary(models3[[4]])

 #disks
vars4 <- c("MAT_1", "MAP_1", "MSP_1" ,"CMI_1", "Locality") #
models4 <- list()



# All additive combinations
for(i in 1:length(vars4)) {
  combos <- combn(vars4, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " + "))
    models4[[length(models4) + 1]] <- glm(as.formula(formula_str), data = subset, family=poisson)
  }
}

# Add interaction models -- 22nd
for(i in 2:length(vars4)) {
  i<-4
  combos <- combn(vars4, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " * "))
    models4[[length(models4) + 1]] <- glm(as.formula(formula_str), data = subset, family=poisson)
  }
}
aictab(models4)
summary(models4[[13]])
