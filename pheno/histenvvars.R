library(tidyverse)
envs <-read.csv("Data/phenotype/emma_gee_results.csv", header=T)
envs <- envs[rowSums(!is.na(envs) & envs != "") > 0, ]

#filter all observations for only phenotype ones
numbs <- as.vector(gt$recordNumber_CollectorNumber)
envssub <- subset(envs, SampleID %in% numbs)

colnames(envssub)[colnames(envssub) == 'SampleID'] <- 'recordNumber_CollectorNumber'

envssubset<-envssub[,c('recordNumber_CollectorNumber', 'terraclimate_climate_water_deficit_mm_10y',
                       'modis_vegetation_indices_senescence_doy_1y',
                       'terraclimate_palmer_drought_severity_index_10y', 'terraclimate_palmer_drought_severity_index_30y',
                       'terraclimate_annual_precip_10y', 'terraclimate_precip_seasonality_10y',
                       'terraclimate_precip_warmest_quarter_10y',
                       'terraclimate_climate_water_deficit_mm_10y', 
                       'terraclimate_precip_seasonality_30y')]

pairs.panels(envssubset, scale=T) 

df_list <- list(gt, envssub)
gtenv <- df_list %>% reduce(full_join, by='recordNumber_CollectorNumber')
