setwd("~/Library/CloudStorage/GoogleDrive-sscheim@ucsc.edu/.shortcut-targets-by-id/1th4k9MOXb4NFZdnAlKncS7yQaC1BQnAT/GaviotaTarplant/Analysis")


library(tidyverse)
library(AICcmodavg)
library(tidyr)
library(dplyr)
library(lme4)
library(stringr)


gt<-read.csv("Data/phenotype/gtinatandvars.csv", header=T)
# Clean up this dataframe
gt <- gt[rowSums(gt == "" | is.na(gt)) < ncol(gt), ]
gt <- gt[rowSums(is.na(gt)) != ncol(gt), ]
#drop weird ones
gt <- gt[gt$recordNumber_CollectorNumber != "707",]
gt <-gt[gt$recordNumber_CollectorNumber != "233000618",]
gt <- gt %>%
  filter(recordNumber_CollectorNumber <=8484 | recordNumber_CollectorNumber >=8501 )
gt <-gt[gt$recordNumber_CollectorNumber != "8481A",]
gt <-gt[gt$recordNumber_CollectorNumber != "8481B",]
gt <- gt %>%
  filter(recordNumber_CollectorNumber <=8501 | recordNumber_CollectorNumber >=8505 )

gt <- gt[, -which(names(gt) == "BIO12")]
gt <- gt[, -which(names(gt) == "BIO15")]
gt <- gt[, -which(names(gt) == "BIO18")]
gt <- gt[, -which(names(gt) == "fog")]
gt <- gt[, -which(names(gt) == "CHILI")]

newvars <- read.csv("Data/phenotype/emma_gee_results_final.csv", header=TRUE)

newvars <- newvars %>% mutate(SampleID = str_remove_all(SampleID, "https://www.inaturalist.org/observations/"))
colnames(newvars)[colnames(newvars) == 'SampleID'] <- 'recordNumber_CollectorNumber'
nv <- subset(newvars, select = c(recordNumber_CollectorNumber, fog_germyr, fog_25yr, CHILI, BIO12_germyr, 
                                 BIO12_25yr, BIO12_65yr, BIO15_germyr, BIO15_25yr, BIO15_65yr, BIO18_germyr, 
                                 BIO18_25yr, BIO18_65yr, Kling_germyr, Kling_65yr))
numbs <- as.vector(gt$recordNumber_CollectorNumber)
nv2 <- subset(nv, recordNumber_CollectorNumber %in% numbs)   
df_list <- list(gt, nv2)
newgt <- df_list %>% reduce(full_join, by='recordNumber_CollectorNumber')

# write.csv(newgt, "~/Library/CloudStorage/GoogleDrive-sscheim@ucsc.edu/.shortcut-targets-by-id/1th4k9MOXb4NFZdnAlKncS7yQaC1BQnAT/GaviotaTarplant/Analysis/Data/phenotype/updatedvarsandpheno.csv", row.names=FALSE)


# rays germ year ######

subset <- subset(newgt, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality, fog_germyr, fog_25yr, CHILI, BIO12_germyr, 
                                   BIO12_25yr, BIO12_98yr, BIO15_germyr, BIO15_25yr, BIO15_98yr, BIO18_germyr, 
                                   BIO18_25yr, BIO18_98yr, Kling_germyr, Kling_65yr, recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd, 
                                diff, burned, number, cause))
modern <- subset(subset, year>2000)
modern <- modern %>% drop_na(Count_ray_flourets)
modern <- modern %>% drop_na(BIO12_germyr)

modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)
modern2$BIO12_scaled <- scale(modern2$BIO12_germyr)
modern2$CHILI_scaled <- scale(modern2$CHILI)
modern2$BIO18_scaled <- scale(modern2$BIO18_germyr)
modern2$BIO15_scaled <- scale(modern2$BIO15_germyr)
modern2$fog_scaled <- scale(modern2$fog_germyr)
vars <- c("BIO12_scaled", "BIO18_scaled", "fog_scaled", "BIO15_scaled", "CHILI_scaled")
models <- list()

# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), data = modern2, family=poisson)
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), control= glmerControl(optimizer="bobyqa",
                                                                                         optCtrl = list(maxfun = 100000)), data = modern2, family=poisson
    )
    
  }
}

# ### REPORT SINGULARITY WARNING
## coefficent biologically reasonable
# see if leaving out random effect changes coefficient (estimate)

aictab(models)
s#mod2 wins
summary(models[[2]])
summary(models[[4]])
summary(models[[11]])
mod <- glmer(Count_ray_flourets ~ BIO12 + (1 | Locality), data=modern2, family=poisson)
# make sure all fits work out the same--you can ignore singularity warning
allFit <- allFit(mod)

dispersion_glmer(mod)
#underdispersed, but poisson log link fine for that

sim_resid <- simulateResiduals(mod2)

plot(sim_resid) 

#resids not good

testZeroInflation(sim_resid) 


# Visual check of predictions
# Does the fitted model predict realistic values?
#fix this
modern2$pred <- predict(mod2)
modern2$pred <- exp(predict(mod2))  # back-transform from log scale
ggplot(modern2, aes(x = Count_ray_flourets) +
         geom_histogram(aes(y = ..density..), alpha = 0.5) +
         geom_density(aes(x = pred_gaussian), color = "blue") +
         geom_density(aes(x = pred_gamma), color = "red") +
         labs(title = "Observed vs predicted distributions"))


# Shapiro-Wilk test on residuals (for small datasets):
rshapiro.test(residuals(m_gaussian))
# p > 0.05 means residuals are normal


#Levene's test for homogeneity of variance:
leveneTest(height ~ factor1, data = dat)
# p > 0.05 means variances are equal across groups


#disks germyr ########

#exclude potentially problematic specimens?

modern <- subset(subset, year>2000)
modern <- modern %>% drop_na(Count_disk_flourets)
modern <- modern %>% drop_na(BIO12_germyr)

modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)
modern2$BIO12_scaled <- scale(modern2$BIO12_germyr)
modern2$CHILI_scaled <- scale(modern2$CHILI)
modern2$BIO18_scaled <- scale(modern2$BIO18_germyr)
modern2$BIO15_scaled <- scale(modern2$BIO15_germyr)
modern2$fog_scaled <- scale(modern2$fog_germyr)

vars2 <- c("BIO18_scaled", "fog_scaled", "BIO12_scaled", "BIO15_scaled", "CHILI_scaled") 
models2 <- list()

# All additive combinations
for(i in 1:length(vars2)) {
  combos <- combn(vars2, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models2[[length(models2) + 1]] <- glmer(as.formula(formula_str), data = modern2, family=poisson)
  }
}

# Add interaction models
for(i in 2:length(vars2)) {
  combos <- combn(vars2, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models2[[length(models2) + 1]] <- glmer(as.formula(formula_str), data = modern2, control= glmerControl(optimizer="bobyqa",
                                                                                                           optCtrl = list(maxfun = 100000)), family=poisson)
  }
}

test <- glmer(Count_disk_flourets ~ BIO18_scaled * BIO12_scaled  * BIO15_scaled + (1 | Locality), data=modern2, family = poisson)
summary(test)
formula(test)
aictab(models2)

fix.formula <- Count_disk_flourets ~ BIO18_scaled * BIO12_scaled  *BIO15_scaled 

X <- model.matrix (fix.formula, modern2)

#mod 32 wins 
summary(models2[[45]])
summary(models2[[37]])
summary(models2[[36]])
#
diskmod <- glmer(Count_disk_flourets ~ BIO18 * fog + (1 | Locality), data=modern2, family = poisson)
# make sure all fits work out the same--you can ignore singularity warning
allFit <- allFit(mod10)

library(blmeco)
dispersion_glmer(mod10)
#between .75 and 1.4
#not overdispersed

sim_resid <- simulateResiduals(mod10)

plot(sim_resid) 

# resids fine

################## proportion #######
library(glmmTMB)

modern <- subset(subset, year>2000)
modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)
modern2$total_florets <- with(modern2, Count_ray_flourets + Count_disk_flourets)
modern2$ray_prop <- with(modern2, Count_ray_flourets/total_florets)
modern2$BIO12_scaled <- scale(modern2$BIO12_germyr)
modern2$CHILI_scaled <- scale(modern2$CHILI)
modern2$BIO18_scaled <- scale(modern2$BIO18_germyr)
modern2$BIO15_scaled <- scale(modern2$BIO15_germyr)
modern2$fog_scaled <- scale(modern2$fog_germyr)
rayprop <- modern2 %>% drop_na(ray_prop)
rayprop <- modern2 %>% drop_na(BIO12_scaled)
rayprop <- rayprop[rayprop$Locality != "Dangermond",]
rayprop$Locality <- as.factor(rayprop$Locality)
vars <- c("BIO12_scaled", "CHILI_scaled", "BIO18_scaled", "BIO15_scaled", "fog_scaled")
models <- list()

#additive
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_prop ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmmTMB(as.formula(formula_str), family=beta_family(link="logit"), data = modern2)
  }
}

# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_prop ~", paste(combo, collapse = " * "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmmTMB(as.formula(formula_str), family=beta_family(link="logit"), data = modern2)
    
  }
}

aictab(models)
summary(models[[3]])
summary(models[[5]])
summary(models[[4]])

plot(residuals(beta) ~ fitted(beta))