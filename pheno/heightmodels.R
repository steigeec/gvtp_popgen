#add histogram of plant heights/widths
#highly colinear, just exploring multidecadal avgs
library(psych)
library(tidyverse)
library(dplyr)

subset<-modern2[,c("fog_scaled", "BIO12_scaled", "CHILI_scaled", "BIO15_scaled")] 
par(mar = c(1, 1, 1, 1))
pairs.panels(subset, scale=T) 


filter <- newgt %>%
  filter(recordNumber_CollectorNumber <=8484 | recordNumber_CollectorNumber>=8501 )
filter <-filter[filter$recordNumber_CollectorNumber != "8481A",]
filter <-filter[filter$recordNumber_CollectorNumber != "8481B",]
filter <- filter %>%
  filter(recordNumber_CollectorNumber <=8501 | recordNumber_CollectorNumber >=8505 )

subset <- subset(filter, select = c(year, Count_disk_flourets, Count_ray_flourets, Locality, BIO12, BIO15, 
                                    BIO18, fog, CHILI, recordNumber_CollectorNumber, Lat._dd.dddd, Lon_dd.dddd, 
                                    diff, burned, number, cause, Height_in, Width_in, Height_width))
subset$BIO12_scaled <- scale(subset$BIO12)
subset$CHILI_scaled <- scale(subset$CHILI)
subset$BIO18_scaled <- scale(subset$BIO18)
subset$BIO15_scaled <- scale(subset$BIO15)
subset$fog_scaled <- scale(subset$fog)
modern <- subset(subset, year>2023)
modern <- modern %>% drop_na(Height_width)


modern$Locality <- as.factor(modern$Locality)
modern2 <- modern[modern$Locality != "CDFW",]
modern2 <- modern2[modern2$Locality != "Dangermond",]
modern2$Locality <- droplevels(modern2$Locality)


######RATIO#####

vars3 <- c("fog_scaled", "BIO12_scaled", "CHILI_scaled", "BIO15_scaled") 
models3 <- list()

# additive models ##
for(i in 1:length(vars3)) {
  combos <- combn(vars3, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Height_width ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models3[[length(models3) + 1]] <- glmer(as.formula(formula_str), data = modern2, family= Gamma(link = log))
  }
} 


# Add interaction models -- 22nd
for(i in 2:length(vars3)) {
  combos <-combn(vars3, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("Height_width ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models3[[length(models3) + 1]] <- glmer(as.formula(formula_str), data = modern2, family=Gamma(link = log), 
                                            control=glmerControl(optimizer="bobyqa",
                                                                 optCtrl = list(maxfun = 100000))
                                        )
  }
}



aictab(models3)
summary(models3[[23]])
summary(models3[[20]])
summary(models3[[18]])
#mod 8 wins
heiwid <- glmer(Height_width ~ fog * BIO18 * BIO12 + (1 | Locality), data=modern2, family=Gamma(link=log))

########### height ##########
library(lmerTest)
vars <- c("fog_scaled", "BIO12_scaled", "CHILI_scaled", "BIO15_scaled")  
models <- list()

# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("log(Height_in) ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- lmer(as.formula(formula_str), data = modern2)
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("log(Height_in) ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- lmer(as.formula(formula_str), data = modern2)
  }
}

aictab(models)
summary(mod26)
summary(models[[25]])
summary(models[[23]])



heightmod <-lmer(log(Height_in) ~ BIO18  + (1 | Locality), data=modern2)
sim_resid <- simulateResiduals(bestmod)
summary(bestmod)

########### width ##########
library(lmerTest)
vars <- c("BIO15_scaled", "fog_scaled", "BIO12_scaled", "CHILI_scaled") 
models <- list()

# additive models ##
for(i in 1:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("log(Width_in) ~", paste(combo, collapse = " + "),  " + (1 | Locality)")
    models[[length(models) + 1]] <- lmer(as.formula(formula_str), data = modern2)
  }
}


# Add interaction models -- 22nd
for(i in 2:length(vars)) {
  combos <- combn(vars, i, simplify = FALSE)
  for(combo in combos) {
    formula_str <- paste("log(Width_in) ~", paste(combo, collapse = " * "), " + (1 | Locality)")
    models[[length(models) + 1]] <- lmer(as.formula(formula_str), data = modern2)
  }
}


aictab(models)

summary(models[[26]])

bestmod <-lmer(log(Width_in) ~ BIO18 + BIO12 + (1 | Locality), data=modern2)
sim_resid <- simulateResiduals(bestmod)
summary(bestmod)
plot(sim_resid)
#hmmmmm

#### graphs histograms ###

filter <- newgt %>%
  filter(recordNumber_CollectorNumber <=8484 | recordNumber_CollectorNumber>=8501 )
filter <-filter[filter$recordNumber_CollectorNumber != "8481A",]
filter <-filter[filter$recordNumber_CollectorNumber != "8481B",]
filter <- filter %>%
  filter(recordNumber_CollectorNumber <=8501 | recordNumber_CollectorNumber >=8505 )
modern <- subset(filter, year>2023)
modern <- modern %>% drop_na(Height_width)

library(ggplot2)
library(cowplot)
heighthisto <- ggplot(modern, aes(x=Height_in, y = ..density..)) + 
  geom_histogram(color="black", fill="lightgoldenrod3", binwidth=1, linewidth=0.75) +
  geom_density(alpha=.2, fill="yellow", linewidth=0.75) + xlab("") +
  ylab("Density") + scale_x_continuous(n.breaks=15) + 
  theme_cowplot() +
  theme(text = element_text(family = "Palatino", size=36), axis.text=element_text(size=34),
        axis.title=element_text(face='bold'))
widthhisto <- ggplot(modern, aes(x=Width_in, y = ..density..)) + 
  geom_histogram(color="black", fill="lightgoldenrod3", binwidth=1, linewidth=0.75) +
  geom_density(alpha=.2, fill="yellow", linewidth=0.75) + xlab("") +
  ylab("Density") + scale_x_continuous(n.breaks=15) + 
  theme_cowplot() +
  theme(text = element_text(family = "Palatino", size=36), axis.text=element_text(size=34),
        axis.title=element_text(face='bold'))
heightwidthhisto <- ggplot(modern, aes(x=Height_width, y = ..density..)) + 
  geom_histogram(color="black", fill="lightgoldenrod3", binwidth=1, linewidth=0.75) +
  geom_density(alpha=.2, fill="yellow", linewidth=0.75) + xlab("") +
  ylab("Density") + scale_x_continuous(n.breaks=15) + 
  theme_cowplot() +
  theme(text = element_text(family = "Palatino", size=36), axis.text=element_text(size=34),
        axis.title=element_text(face='bold'))

library(patchwork)

png("Figures/phenotype/heightwidthhistograms.png", width=2600, height=700)
heighthisto + widthhisto + heightwidthhisto + plot_annotation(tag_levels = 'a') &
 theme(plot.tag = element_text(size = 50),
        plot.tag.position  = c(.07, .04))
dev.off() 