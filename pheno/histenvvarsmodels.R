#modern and historic specimens
library(AICcmodavg)
library(tidyr)
library(lme4)


subset <- subset(newgt, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality, fog_germyr, fog_25yr, CHILI, BIO12_germyr, 
                                   BIO12_25yr, BIO12_65yr, BIO15_germyr, BIO15_25yr, BIO15_65yr, BIO18_germyr, 
                                   BIO18_25yr, BIO18_65yr, Kling_germyr, Kling_65yr, recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd, 
                                   diff, burned, number, cause))
                                   #diff, alarmdate, burned, acreage, Height_in, Width_in, Height_width, Elevation._masl))

#### rays decadal ####

rays <- subset %>% drop_na(Count_ray_flourets)
rays <- rays[rays$Locality != "Dangermond",]
rays$Locality <- as.factor(rays$Locality)
rays$BIO12_scaled <- scale(rays$BIO12_65yr)
rays$CHILI_scaled <- scale(rays$CHILI)
rays$BIO18_scaled <- scale(rays$BIO18_65yr)
rays$BIO15_scaled <- scale(rays$BIO15_65yr)
rays$fog_scaled <- scale(rays$Kling_65yr)
vars <- c("BIO12_scaled", "BIO18_scaled", "fog_scaled", "BIO15_scaled", "CHILI_scaled")
models <- list()
# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), control = glmerControl(optimizer="bobyqa",
                                                                                         optCtrl = list(maxfun = 100000)), data = rays, family=poisson)
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), control= glmerControl(optimizer="bobyqa",
                                                                                                  optCtrl = list(maxfun = 100000)), data = rays, family=poisson
    )
    
  }
}

# ### REPORT SINGULARITY WARNING
## coefficent biologically reasonable
# see if leaving out random effect changes coefficient (estimate)

aictab(models)
summary(models[[3]])
summary(models[[14]])
summary(models[[7]])



#### disks ####

disks <- subset %>% drop_na(Count_disk_flourets)
disks <- disks[disks$Locality != "Dangermond",]
disks$Locality <- as.factor(disks$Locality)
disks$BIO12_scaled <- scale(disks$BIO12_65yr)
disks$CHILI_scaled <- scale(disks$CHILI)
disks$BIO18_scaled <- scale(disks$BIO18_65yr)
disks$BIO15_scaled <- scale(disks$BIO15_65yr)
disks$fog_scaled <- scale(disks$Kling_65yr)
vars <- c("BIO12_scaled", "BIO18_scaled", "fog_scaled", "BIO15_scaled", "CHILI_scaled")
models <- list()
# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), control=glmerControl(optimizer="bobyqa",
                                          optCtrl = list(maxfun = 100000)), data = disks, family=poisson)
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), control=glmerControl(optimizer="bobyqa",
                                                                                        optCtrl = list(maxfun = 100000)), data = disks, family=poisson
    )
    
  }
}

# ### REPORT SINGULARITY WARNING



aictab(models)
summary(models[[13]])
summary(models[[19]])
summary(models[[22]])

### ray prop ###



library(glmmTMB)
subset$total_florets <- with(subset, Count_ray_flourets + Count_disk_flourets)
subset$ray_prop <- with(subset, Count_ray_flourets/total_florets)
rayprop <- subset %>% drop_na(ray_prop)
rayprop <- rayprop[rayprop$Locality != "Dangermond",]
rayprop$Locality <- as.factor(rayprop$Locality)
rayprop$BIO12_scaled <- scale(rayprop$BIO12_65yr)
rayprop$CHILI_scaled <- scale(rayprop$CHILI)
rayprop$BIO18_scaled <- scale(rayprop$BIO18_65yr)
rayprop$BIO15_scaled <- scale(rayprop$BIO15_65yr)
rayprop$fog_scaled <- scale(rayprop$Kling_65yr)
vars <- c("BIO12_scaled", "BIO18_scaled", "fog_scaled", "BIO15_scaled", "CHILI_scaled")
models <- list()

#additive
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_prop ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmmTMB(as.formula(formula_str), family=beta_family(link="logit"), data = rayprop)
  }
}

# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_prop ~", paste(combo, collapse = " * "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmmTMB(as.formula(formula_str), family=beta_family(link="logit"), data = rayprop)
    
  }
}
aictab(models)
summary(models[[11]])
summary(models[[8]])
summary(models[[15]])

