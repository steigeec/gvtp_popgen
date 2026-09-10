
# read in inat rays
inatrays <- read.csv("Data/phenotype/inatrays.csv", header=T, fill=T)
inatrays <- subset(inatrays, select = c(Year, iNat_ray_count, latitude, longitude, 
                                        inat_Specimen, Locality))
colnames(inatrays)[colnames(inatrays) == 'longitude'] <- 'Lon_dd.dddd'
colnames(inatrays)[colnames(inatrays) == 'latitude'] <- 'Lat._dd.dddd'
colnames(inatrays)[colnames(inatrays) == 'Year'] <- 'year'
colnames(inatrays)[colnames(inatrays) == 'inat_Specimen'] <- 'recordNumber_CollectorNumber'
#not sure if i should do this next one
colnames(inatrays)[colnames(inatrays) == 'iNat_ray_count'] <- 'Count_ray_flourets'
inatrays$year <- as.integer(inatrays$year)
inatrays$recordNumber_CollectorNumber <- as.character(inatrays$recordNumber_CollectorNumber)
inatrays <- inatrays %>% drop_na(Count_ray_flourets)
 inatrays <- inatrays %>% drop_na(iNat_ray_count)
inatrays <- inatrays %>% drop_na(Lat._dd.dddd)

gt4<-read.table("Data/phenotype/gvtp_phenotype_13feb2026.tsv", sep='\t', header=T)

#combine gt and inat
gt <- bind_rows(gt, inatrays)
#merged



