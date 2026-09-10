setwd("~/Library/CloudStorage/GoogleDrive-sscheim@ucsc.edu/.shortcut-targets-by-id/1th4k9MOXb4NFZdnAlKncS7yQaC1BQnAT/GaviotaTarplant/Analysis")

library(AICcmodavg)
library(tidyr)
library(dplyr)
library(lme4)


##UPDATE WITH NEW VARS ####
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

# rays ######



subset <- subset(newgt, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality, BIO12, BIO15, 
                                    BIO18, fog, CHILI, recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd, 
                                    diff, burned, number, cause))
modern <- subset(subset, year>2000)
modern <- modern %>% drop_na(Count_ray_flourets)
modern <- modern %>% drop_na(BIO12)

modern$Locality <- as.factor(modern$Locality)
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)
modern2$BIO12_scaled <- scale(modern2$BIO12)
modern2$CHILI_scaled <- scale(modern2$CHILI)
modern2$BIO18_scaled <- scale(modern2$BIO18)
modern2$BIO15_scaled <- scale(modern2$BIO15)
modern2$fog_scaled <- scale(modern2$fog)
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
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), control= glmerControl(autoscale=TRUE, optimizer="bobyqa",
                                                                                         optCtrl = list(maxfun = 100000)), data = modern2, family=poisson
                                  )
                                          
  }
}

# ### REPORT SINGULARITY WARNING
## coefficent biologically reasonable
# see if leaving out random effect changes coefficient (estimate)

aictab(models)
#mod2 wins
summary(models[[2]])
summary(models[[37]])
summary(models[[32]])
mod2 <- glmer(Count_ray_flourets ~ BIO18 + (1 | Locality), data=modern2, family=poisson)
# make sure all fits work out the same--you can ignore singularity warning
allFit <- allFit(mod2)

dispersion_glmer(mod2)
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


#disks  ########

#exclude potentially problematic specimens?

modern <- subset(subset, year>2000)
modern <- modern %>% drop_na(Count_disk_flourets)

modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)
modern2$BIO12_scaled <- scale(modern2$BIO12)
modern2$CHILI_scaled <- scale(modern2$CHILI)
modern2$BIO18_scaled <- scale(modern2$BIO18)
modern2$BIO15_scaled <- scale(modern2$BIO15)
modern2$fog_scaled <- scale(modern2$fog)



vars2 <- c("BIO18_scaled", "fog_scaled", "BIO12_scaled", "BIO15_scaled", "CHILI_scaled") #, 
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
aictab(models2)

#mod 32 wins 
summary(models2[[36]])
summary(models2[[45]])
summary(models2[[37]])
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


ggplot(modern, aes(x=Count_disk_flourets, y=fog)) + geom_point()

### burned ###

burns <- filter %>%
  filter(burned=="YES")

vars5 <- c("BIO18", "Locality", "number", "diff", "CHILI") #
models5 <- list()

# All additive combinations
for(i in 1:length(vars5)) {
  combos <- combn(vars5, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " + "))
    models5[[length(models5) + 1]] <- glm(as.formula(formula_str), data = burns, family=poisson)
  }
}

# Add interaction models -- 22nd
for(i in 2:length(vars5)) {
  combos <- combn(vars5, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " * "))
    models5[[length(models5) + 1]] <- glm(as.formula(formula_str), data = burns, family=poisson)
  }
}
aictab(models5)

################## proportion #######

library(glmmTMB)

subset <- subset(newgt, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality, BIO12, BIO15, 
                                BIO18, fog, CHILI, recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd, 
                                diff, burned, number, cause))

modern <- subset(subset, year>2000)
modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)
modern2$total_florets <- with(modern2, Count_ray_flourets + Count_disk_flourets)
modern2$ray_prop <- with(modern2, Count_ray_flourets/total_florets)
modern2$BIO12_scaled <- scale(modern2$BIO12)
modern2$CHILI_scaled <- scale(modern2$CHILI)
modern2$BIO18_scaled <- scale(modern2$BIO18)
modern2$BIO15_scaled <- scale(modern2$BIO15)
modern2$fog_scaled <- scale(modern2$fog)
rayprop <- modern2 %>% drop_na(ray_prop)
rayprop <- modern2 %>% drop_na(BIO12)
rayprop <- rayprop[rayprop$Locality != "Dangermond",]
rayprop$Locality <- as.factor(rayprop$Locality)

vars <- c("BIO12_scaled", "CHILI_scaled", "BIO18_scaled", "BIO15_scaled", "fog_scaled")
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
summary(models[[3]])
summary(models[[5]])
summary(models[[2]])


hist(newgt$BIO18)

plot(residuals(beta) ~ fitted(beta))

modern2$ray_disk_ratio <- with(modern2, Count_ray_flourets/Count_disk_flourets)


vars <- c("BIO18", "fog", "BIO12") 
models <- list()

# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_disk_ratio  ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), data = modern2, family=Gamma(link='log'))
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_disk_ratio ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), data = modern2, Gamma(link='log'))
  }
}

aictab(models)
summary(models[[6]])

####

vars <- c("BIO12", "BIO18", "fog") 
models <- list()

#additive
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_prop ~", paste(combo, collapse = " + "))
    models[[length(models) + 1]] <- betareg(as.formula(formula_str), data = modern2)
  }
}

# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("ray_prop ~", paste(combo, collapse = " * "))
    models[[length(models) + 1]] <- betareg(as.formula(formula_str), data = modern2)
  
  }
}
aictab(models)


####### stats #####

library(gtsummary)

png("Figures/phenotype/diskmodelstats.png", width=700, height=400)
tbl_regression(mod10, exponentiate = TRUE)
dev.off()

png("Figures/phenotype/raymodelstats.png", width=700, height=400)
tbl_regression(mod2, exponentiate = TRUE)
dev.off()

############## visualize #########


library(coefplot2)
library(sjPlot)

coefplot2(mod2)

coefplot2(mod10)

plot_model(mod10, type = "int", terms = c("fog", "BIO18"), transform='log')


modern2$pred <- predict(mod10)
modern2$pred <- exp(predict(mod10)) 

png("Figures/phenotype/diskmodelobsvspred.png", width=700, height=400)
ggplot(modern2, aes(x = Count_disk_flourets)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, bins=30, binwidth=1) +
  geom_density(aes(x = pred), color = "blue") +
  labs(title = "Observed vs predicted distributions") +
  xlab('Count disk florets') +
  ylab('Proportion') +
  theme_cowplot(font_family='Palatino')
dev.off()

ggplot(modern2, aes(x = Count_disk_flourets, y=fog)) +
  geom_point(alpha = 0.5) +
  geom_density(aes(x = pred), color = "blue") +
  labs(title = "Observed vs predicted distributions")

modern2$predray <- predict(mod2, re.form=NULL)
modern2$predray <- exp(predict(mod2, re.form=NULL)) 
ggplot(modern2, aes(x = Count_ray_flourets)) +
  geom_histogram(aes(y = ..density..), alpha = .75, bins=10) +
  geom_density(aes(x = pred), color = "blue") +
  labs(title = "Observed vs predicted distributions") + 
  xlab('Count disk florets')+
  ylab('Proportion') +
  theme_cowplot(font_family='Palatino')


#### ####

mod10noloc <- glm(Count_disk_flourets ~ BIO18 * fog, data=modern2, family = poisson)

grid <- expand.grid(
  fog = seq(min(modern2$fog), max(modern2$fog), length = 100),
  BIO18 = seq(min(modern2$BIO18), max(modern2$BIO18), length = 100)
)
grid$pred <- predict(mod10noloc, newdata = grid, type = "response")

library(geometry) # for convex hull

envars <- drop_na(gt, fog)
envars <- drop_na(envars, BIO18)
# Get convex hull indices
hull_idx <- convhulln(as.matrix(envars[, c("fog", "BIO18")]))

# Function to test if points are inside the hull
in_hull <- function(points, hull_pts, hull) {
  inhulln(hull, as.matrix(points))
}

# Mask grid to convex hull
hull_mat <- as.matrix(envars[, c("fog", "BIO18")])
hull <- convhulln(hull_mat)

grid$in_hull <- inhulln(hull, as.matrix(grid[, c("fog", "BIO18")]))
grid$pred_masked <- ifelse(grid$in_hull, grid$pred, NA)

png("Figures/phenotype/fogbio18disks.png", width=1200, height=900)
diskcontour <- ggplot() +
  geom_contour_filled(grid, mapping=aes(fog, BIO18, z = pred_masked), na.rm=T) +
  geom_point(data = gt,
             aes(fog, BIO18),
             size = 5) +
  theme_cowplot(font_family = 'Palatino', font_size = 36) + theme(legend.title = element_blank(),
                                                                  axis.title=element_text(face='bold'),
                                                                  axis.text=element_text(size=34),
                                                                  legend.text=element_text(size=32),
                                                                  legend.key.size=unit(1, 'cm')) +
  coord_cartesian(expand=F) + xlab('Fog') + scale_fill_viridis(option='viridis', discrete=T)
dev.off()


beta <- betareg(ray_prop ~ BIO12 + fog, data=modern2)

grid <- expand.grid(
  fog = seq(min(modern2$fog), max(modern2$fog), length = 100),
  BIO12 = seq(min(modern2$BIO12), max(modern2$BIO12), length = 100)
)
grid$pred <- predict(beta, newdata = grid, type = "response")

envars <- drop_na(gt, BIO12)
envars <- drop_na(envars, fog)
# Get convex hull indices
hull_idx <- convhulln(as.matrix(envars[, c("BIO12", "fog")]))

# Function to test if points are inside the hull
in_hull <- function(points, hull_pts, hull) {
  inhulln(hull, as.matrix(points))
}

# Mask grid to convex hull
hull_mat <- as.matrix(envars[, c("BIO12", "fog")])
hull <- convhulln(hull_mat)

grid$in_hull <- inhulln(hull, as.matrix(grid[, c("BIO12", "fog")]))
grid$pred_masked <- ifelse(grid$in_hull, grid$pred, NA)

png("Figures/phenotype/fogbio12rayprop.png", width=1200, height=900)
# raypropcontour <- 
raypropcontour <- ggplot() +
  geom_contour_filled(data=grid, mapping=aes(fog, BIO12, z = pred_masked), na.rm=T) +
  geom_point(data = gt,
             aes(fog, BIO12),
             size = 5)  +
  theme_cowplot(font_family = 'Palatino') + theme(legend.title = element_blank(),
                                                                  axis.title=element_text(face='bold', size=36), 
                                                  axis.text=element_text(size=34),
                                                  legend.text=element_text(size=30),
                                                  legend.key.size=unit(1, 'cm')) +
  coord_cartesian(expand=F) + xlab('Fog') + 
  scale_fill_viridis(option='mako', discrete=T)
dev.off()

##### with new env vars ####################

#modern and historic specimens
library(AICcmodavg)
library(tidyr)
library(lme4)


subset <- subset(gtenv, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality, BIO12, BIO15, 
                                   BIO18, CHILI, recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd
))
#diff, alarmdate, burned, acreage, Height_in, Width_in, Height_width, Elevation._masl))

modern <- subset(subset, year>2000)
modern <- modern %>% drop_na(Count_ray_flourets)
modern <- modern %>% drop_na(BIO12)

modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)

#### rays ####

rays <- subset %>% drop_na(Count_ray_flourets)

rays$Locality <- as.factor(rays$Locality)
vars <- c("BIO12", "CHILI", "BIO18")
models <- list()
# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), data = rays, family=poisson)
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_ray_flourets ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), data = rays, family=poisson
    )
    
  }
}

# ### REPORT SINGULARITY WARNING
## coefficent biologically reasonable
# see if leaving out random effect changes coefficient (estimate)

aictab(models)
summary(models[[9]]) # ns
summary(models[[8]])
summary(models[[10]])

#### disks ####
modern <- subset(subset, year>2000)
modern <- modern %>% drop_na(Count_disk_flourets)

modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)

vars <- c("BIO12", "BIO18", "CHILI")
models <- list()
# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), data = modern2, family=poisson)
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Count_disk_flourets ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- glmer(as.formula(formula_str), data = disks, family=poisson
    )
    
  }
}

# ### REPORT SINGULARITY WARNING

## coefficent biologically reasonable
# see if leaving out random effect changes coefficient (estimate)

aictab(models)
summary(models[[7]])
summary(models[[4]])
summary(models[[2]])
