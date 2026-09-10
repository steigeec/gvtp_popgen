setwd("~/Library/CloudStorage/GoogleDrive-sscheim@ucsc.edu/.shortcut-targets-by-id/1th4k9MOXb4NFZdnAlKncS7yQaC1BQnAT/GaviotaTarplant/Analysis")
library(terra)
library(sf)
library(janitor)
library(tidyr)
library(dplyr)
library(psych)
library(geodata)
library(raster)
library(paletteer)
library(viridis)
library(moments)
library(ggplot2)
library(cowplot)
gt <-read.table("Data/phenotype/gvtp_phenotype_13feb2026.tsv", sep='\t', header=T)
# removing blank rows and columns
gt <- gt[rowSums(is.na(gt)) < ncol(gt), colSums(is.na(gt)) < nrow(gt)]
# Get rid of potentially duplicated rows…
gt<-unique(gt)
nrow(gt)
gt <-gt[gt$recordNumber_CollectorNumber != "233000618",]

# Make a SpatVector point file for use in extraction
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- vect(gt[,c(3, 14, 15)], geom=c("Lon_dd.dddd", "Lat._dd.dddd"), crs=proj, keepgeom=FALSE)
bbox_extent <- ext(-120.7, -120.1, 34.4, 34.9)
box = c(xmin=-120.7, xmax=-120.1, ymin=34.4, ymax=34.9)

my_sf <- st_read("Data/phenotype/Fire/California_Fire_Perimeters/California_Fire_Perimeters_(all).shp", quiet = TRUE, stringsAsFactors = FALSE, as_tibble = TRUE)
# plot(st_geometry(my_sf), col = "#f2f2f2", bg = "skyblue", lwd = 0.25, border = 0)

# define projection for fire shapefile
my_sf<-st_transform(my_sf, crs=proj)

# fix invalid geometries
my_sf<-st_make_valid(my_sf)


cropped <- st_crop(my_sf
                    , xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)
# plot(st_geometry(cropped), col = "#f2f2f2", bg = "skyblue", lwd = 0.25, border = 0)

# Convert the SpatVector of points to an sf object
sf_points <- sf::st_as_sf(points)
# Verify the class of the new object
class(sf_points)
# slow # fire<- st_intersection(my_sf, sf_points)



#fire <- st_intersects(my_sf, sf_points)
#fire<-unique(as.data.frame(fire))
#names(fire)<-c("recordNumber_CollectorNumber", "fire")
#gtfire <- merge(gt, fire, by="recordNumber_CollectorNumber", all=T)

fire2 <- st_join(my_sf, sf_points)
fire2 <- fire2 %>% drop_na(recordNumber_CollectorNumber)

#prescribed burns
presc <- st_read("Data/phenotype/Fire/prescribed_burns/ds397.gdb")

# define projection for fire shapefile
presc<-st_transform(presc, crs=proj)

#fix error
presc <- st_cast(presc, "MULTIPOLYGON")
# fix invalid geometries
presc<-st_make_valid(presc)
sf_points <- sf::st_as_sf(points)
cropped <- st_crop(presc
                   , xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)
prescburn <- st_join(cropped, sf_points)
prescburn <- prescburn %>% drop_na(recordNumber_CollectorNumber)

#firedf <- data.frame(matrix(ncol = 4, nrow = 0))
#colnames(firedf) <- c('unique_number', 'fire_date', 'recent', 'decadeprior')
gt$diff <- NA
gt$alarmdate<-""
gt$burned <- ""
gt$acreage <- ""
gt$cause <- ""
gt$number <- ""

fire2 <- st_drop_geometry(fire2[, c("recordNumber_CollectorNumber", "ALARM_DATE", "GIS_ACRES", "CAUSE", "YEAR_")])
prescburn <- st_drop_geometry(prescburn[, c("recordNumber_CollectorNumber", "START_DATE", "GIS_ACRES", "YEAR_")])
prescburn$CAUSE <- "0"
names(prescburn)[2] <- "ALARM_DATE"
fire2 <- rbind(fire2, prescburn)


# for loop
for(SAMPLE in 1:length(unique(gt$recordNumber_CollectorNumber))) {
  # sample
  s <- gt$recordNumber_CollectorNumber[SAMPLE] 
  #fire year
  y <- as.numeric(unique(fire2$YEAR_[fire2$recordNumber_CollectorNumber==s]))
  #collection year
  c <- as.numeric(unique(gt$year[gt$recordNumber_CollectorNumber==s]))
  #difference
  diff <- y - c
  print(diff)
  #choose most recent fire before collection date
 firediff <- if (any(diff < 0)) (max(diff[diff < 0])) else NA
 print(firediff)
 gt$diff[gt$recordNumber_CollectorNumber==s] <- firediff
 
 # grab alarm date for fire
 index <- if (any(diff < 0)) which(diff==max(diff[diff < 0])) else NA
 fireyr <- y[index]
 ALARM_DATE <- fire2$ALARM_DATE[fire2$recordNumber_CollectorNumber==s & fire2$YEAR_==fireyr]
 #replace NAs
 index2 <- if(any(!is.na(ALARM_DATE))) which(ALARM_DATE>0) else NA
 alarm <- ifelse((!is.na(index2)), as.character(ALARM_DATE[index2]), NA)
 gt$alarmdate[gt$recordNumber_CollectorNumber==s] <- as.character(alarm)
 
 #create value for burned vs unburned
 burned <- ifelse(!is.na(gt$diff), "YES", "NO")
 gt$burned <- as.character(burned)
 
 #add acreage - proxy for intensity?
 a <- fire2$GIS_ACRES[fire2$recordNumber_CollectorNumber==s & fire2$ALARM_DATE==ALARM_DATE]
 index3 <- if(any(!is.na(a))) which(a>0) else NA
 acre <- ifelse((!is.na(index3)), as.character(a[index3]), NA)
 gt$acreage[gt$recordNumber_CollectorNumber==s] <- as.character(acre)
 
 #add cause
 cau <- fire2$CAUSE[fire2$recordNumber_CollectorNumber==s & fire2$ALARM_DATE==ALARM_DATE]
 if(identical(cau, character(0))){
   cau <- NA
 }
 gt$cause[gt$recordNumber_CollectorNumber==s] <- as.character(cau)
 
 #count times burned
 if(identical(diff, numeric(0))){
   count <- NA
 } else {
   count <- length(diff[diff<0])
 }
 gt$number[gt$recordNumber_CollectorNumber==s] <- count
 
 
}
#replace nas with -250
gt$diff <- replace_na(gt$diff, -250)

########### repeat w merged dataframe #########

points2 <- vect(merged[,c(5, 3, 4)], geom=c("Lon_dd.dddd", "Lat._dd.dddd"), crs=proj, keepgeom=FALSE)
sf_points2 <- sf::st_as_sf(points2)
# Verify the class of the new object
class(sf_points2)

fire3 <- st_join(my_sf, sf_points2)
fire3 <- fire3 %>% drop_na(recordNumber_CollectorNumber)

#prescribed burns
prescburn2 <- st_join(cropped, sf_points2)
prescburn2 <- prescburn2 %>% drop_na(recordNumber_CollectorNumber)

#firedf <- data.frame(matrix(ncol = 4, nrow = 0))
#colnames(firedf) <- c('unique_number', 'fire_date', 'recent', 'decadeprior')
merged$diff <- NA
merged$alarmdate<-""
merged$burned <- ""
merged$acreage <- ""
merged$cause <- ""
merged$number <- ""

fire3 <- st_drop_geometry(fire3[, c("recordNumber_CollectorNumber", "ALARM_DATE", "GIS_ACRES", "CAUSE", "YEAR_")])
prescburn2 <- st_drop_geometry(prescburn2[, c("recordNumber_CollectorNumber", "START_DATE", "GIS_ACRES", "YEAR_")])
prescburn2$CAUSE <- "0"
names(prescburn2)[2] <- "ALARM_DATE"
fire3 <- rbind(fire3, prescburn2)

# for loop
for(SAMPLE in 1:length(unique(merged$recordNumber_CollectorNumber))) {
  # sample
  s <- merged$recordNumber_CollectorNumber[SAMPLE] 
  #fire year
  y <- as.numeric(unique(fire3$YEAR_[fire3$recordNumber_CollectorNumber==s]))
  #collection year
  c <- as.numeric(unique(merged$year[merged$recordNumber_CollectorNumber==s]))
  #difference
  diff <- y - c
  print(diff)
  #choose most recent fire before collection date
  firediff <- if (any(diff < 0)) (max(diff[diff < 0])) else NA
  print(firediff)
  merged$diff[merged$recordNumber_CollectorNumber==s] <- firediff
  
  # grab alarm date for fire
  index <- if (any(diff < 0)) which(diff==max(diff[diff < 0])) else NA
  fireyr <- y[index]
  ALARM_DATE <- fire3$ALARM_DATE[fire3$recordNumber_CollectorNumber==s & fire3$YEAR_==fireyr]
  #replace NAs
  index2 <- if(any(!is.na(ALARM_DATE))) which(ALARM_DATE>0) else NA
  alarm <- ifelse((!is.na(index2)), as.character(ALARM_DATE[index2]), NA)
  merged$alarmdate[merged$recordNumber_CollectorNumber==s] <- as.character(alarm)
  
  #create value for burned vs unburned
  burned <- ifelse(!is.na(merged$diff), "YES", "NO")
  merged$burned <- as.character(burned)
  
  #add acreage - proxy for intensity?
  a <- fire3$GIS_ACRES[fire3$recordNumber_CollectorNumber==s & fire3$ALARM_DATE==ALARM_DATE]
  index3 <- if(any(!is.na(a))) which(a>0) else NA
  acre <- ifelse((!is.na(index3)), as.character(a[index3]), NA)
  merged$acreage[merged$recordNumber_CollectorNumber==s] <- as.character(acre)
  
  #add cause
  cau <- fire3$CAUSE[fire3$recordNumber_CollectorNumber==s & fire3$ALARM_DATE==ALARM_DATE]
  if(identical(cau, character(0))){
    cau <- NA
  }
  merged$cause[merged$recordNumber_CollectorNumber==s] <- as.character(cau)
  
  #count times burned
  if(identical(diff, numeric(0))){
    count <- NA
  } else {
    count <- length(diff[diff<0])
  }
  merged$number[merged$recordNumber_CollectorNumber==s] <- count
  
  
}


### map ########

cropped <- st_crop(my_sf
                   , xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)
cropped2 <- st_crop(presc
                    , xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)
graphing <- st_join(cropped, cropped2)

DEM<-rast("Data/Maps/GIS/elevation/elevation.tif")
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(DEM) <-proj
# Define the cropping extent (the extent of the gaviota tarplant study landscape)
crop_extent <- st_bbox(c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9), crs = st_crs(proj))
# crop the coast off the DEM raster
CA <- st_read("Data/Maps/GIS/state/gshhg-shp-2.3.7/GSHHS_shp/GVTP_clip.shp")
CA <- st_make_valid(CA)  # Ensure geometry validity
CA <- st_transform(CA, crs = crs(DEM))
crop_extent_sf <- st_as_sfc(st_bbox(crop_extent))
CA <- st_crop(CA, crop_extent)
CA<- st_simplify(CA)
# Convert to terra SpatVector
CA_vect <- vect(CA)
# Mask the DEM using CA
DEM <- mask(DEM, CA_vect)
DEM <- crop(DEM, extent(c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)))
DEM_df <- as.data.frame(DEM, xy = TRUE)
sf_points <- sf_points %>%
  filter(recordNumber_CollectorNumber <=8484 | recordNumber_CollectorNumber>=8501 )


png("Figures/phenotype/fireperimslegend.png", width=1800, height=1800)
ggplot() + geom_raster(data=DEM_df, aes(x=x, y=y), fill='grey90') +
  geom_sf(data=graphing, aes(fill=as.numeric(YEAR_.x)), lwd=0) +
  scale_fill_viridis_c(option = "inferno", limits = c(1914, 2024),
                       breaks = c(1914, 1969, 2024),
                       labels = c(1914, 1969, 2024)) +
  theme_cowplot() +  theme(text = element_text(family = "Palatino"),
                           legend.position = 'right',
                           legend.key.height=unit(5, "cm"),
                           legend.key.width = unit(2, 'cm'),
                           legend.box = "horizontal",
                           legend.ticks = element_blank(),
                           axis.text = element_text(size=40, vjust=-0.2),
                           axis.title = element_blank(),
                           plot.margin = margin(1.2,4,1.2,1.2, "cm"),
                           legend.text = element_text(size=40)
                           ) +
  labs(fill='') +
  annotation_scale(text_cex=3, text_family='Palatino', pad_y = unit(1, "cm"), pad_x = unit(1, "cm")) +
  annotation_north_arrow(location='tl', which_north='true', height = unit(4.5, 'cm'), width = unit(3, 'cm'),
                         pad_y = unit(1, "cm"), pad_x = unit(1, "cm"),
                         style=north_arrow_orienteering(text_family='Palatino', text_size=35)) +
  theme(panel.background = element_rect(fill = "grey80")) +
  geom_sf(data=sf_points, size=9, shape=19, color='black') + coord_sf(xlim = c(-120.7, -120.1),
                                                                      ylim = c(34.4, 34.9),
                                                                      expand = FALSE)
dev.off()


#inset map 
png("Figures/phenotype/inset.png", width=1800, height=1800)
ggplot(data = camap,
                mapping = aes(x = long, y = lat)) + geom_polygon(fill='grey0') +
  coord_sf(crs= 4326) + 
  geom_rect(
    xmin = -120.7,
    ymin = 34.4,
    xmax = -120.1,
    ymax = 34.9,
    fill = NA, 
    colour = "black",
    size = 0.8) +
  theme(panel.background = element_rect(fill = "grey80")) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  guides(fill = "none")
dev.off()

### exploratory plots #####
#replace nas with -250
gt$diff <- replace_na(gt$diff, -250)
burns <- subset(gt, Locality=='VAFBArguello' | Locality=='CDFW')



ggplot(burns, aes(x=burned, y=Count_disk_flourets, fill=burned)) +
  geom_boxplot() + scale_fill_manual(values=c("darkgoldenrod1", "darkorange3")) +
  cowplot::theme_cowplot(font_family='Palatino')

#exploratory plots
library(ggplot2)
  
one <- ggplot(burns, aes(x=diff, y=Count_ray_flourets)) +
  geom_point(size=2) + geom_smooth(method='lm', se=F, color="darkorange3") + 
  xlab("Years since burned") + ylab("Count ray florets") + 
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) +
  theme(axis.title=element_text(face='bold'))
two <- burns %>% dplyr::filter(diff > -250 ) %>% ggplot(., aes(x=diff, y=Count_ray_flourets)) +
  geom_point(size=2) + geom_smooth(method='lm', se=F, color="darkorange3") +
  xlab("Years since burned") + ylab("Count ray florets") +
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) + 
  scale_x_continuous(limits=c(-25,-12), n.breaks=10) + 
  theme(axis.title=element_text(face='bold'))
three <- ggplot(burns, aes(x=burned, y=Count_ray_flourets, fill=burned)) +
  geom_boxplot() + scale_fill_manual(values=c("darkgoldenrod1", "darkorange3")) +
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) + 
  theme(legend.position = "none") + xlab('Burned') + ylab('Count ray florets') +
  theme(axis.title=element_text(face='bold'))

png("Figures/phenotype/fireexploratoryplots2.png", width=1700, height=600)
cowplot::plot_grid(one, two, three, align="v", ncol=3,
                   labels = c('a', 'b', 'c'), 
                   label_fontfamily = "Palatino",
                   label_size = 40,
                   label_x = 0, label_y = 0, hjust = -2, vjust = -2)
dev.off()

one <- ggplot(burns, aes(x=diff, y=Count_disk_flourets)) +
  geom_point(size=2) + geom_smooth(method='lm', se=F, color="darkorange3") +
  xlab("Years since burned") + ylab("Count disk florets") +
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) +
  theme(axis.title=element_text(face='bold'))

two <- burns %>% dplyr::filter(diff > -250 ) %>% ggplot(., aes(x=diff, y=Count_disk_flourets)) +
  geom_point(size=2) + geom_smooth(method='lm', se=F, color="darkorange3") +
  xlab("Years since burned") + ylab("Count disk florets") +
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) + scale_x_continuous(limits=c(-25,-12), n.breaks=10) +
  theme(axis.title=element_text(face='bold'))

three <- ggplot(burns, aes(x=burned, y=Count_disk_flourets, fill=burned)) +
  geom_boxplot() + scale_fill_manual(values=c("darkgoldenrod1", "darkorange3")) +
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) +
  theme(legend.position = "none") + xlab('Burned') + ylab('Count disk florets') +
  theme(axis.title=element_text(face='bold'))

png("Figures/phenotype/fireexploratoryplots.png", width=1700, height=600)
cowplot::plot_grid(one, two, three, align="v", ncol=3,
                   labels = c('a', 'b', 'c'), 
                   label_fontfamily = "Palatino",
                   label_size = 40,
                   label_x = 0, label_y = 0, hjust = -2, vjust = -2)
dev.off()

burns$total_florets <- with(burns, Count_ray_flourets + Count_disk_flourets)
burns$ray_prop <- with(burns, Count_ray_flourets/total_florets)
one <- ggplot(burns, aes(x=diff, y=ray_prop)) +
  geom_point(size=2) + geom_smooth(method='lm', se=F, color="darkorange3") + 
  xlab("Years since burned") + ylab("Proportion rays/total florets") + 
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30)  +
  theme(axis.title=element_text(face='bold'))
two <-  burns %>% dplyr::filter(diff > -250 ) %>% ggplot(., aes(x=diff, y=Count_disk_flourets)) +
  geom_point(size=2) + geom_smooth(method='lm', se=F, color="darkorange3") +
  xlab("Years since burned") + ylab("Proportion rays/total florets") +
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) + 
  scale_x_continuous(limits=c(-25,-12), n.breaks=10)  +
  theme(axis.title=element_text(face='bold'))
three <- ggplot(burns, aes(x=burned, y=ray_prop, fill=burned)) +
  geom_boxplot() + scale_fill_manual(values=c("darkgoldenrod1", "darkorange3")) +
  cowplot::theme_cowplot(font_family='Palatino', font_size = 30) + 
  theme(legend.position = "none") + xlab('Burned') + ylab('Proportion rays/total florets')  +
  theme(axis.title=element_text(face='bold'))
png("Figures/phenotype/fireexploratoryplots3.png", width=1700, height=600)
cowplot::plot_grid(one, two, three, align="v", ncol=3,
                   labels = c('a', 'b', 'c'), 
                   label_fontfamily = "Palatino",
                   label_size = 40,
                   label_x = 0, label_y = 0, hjust = -2, vjust = -2)
dev.off()

burns$Locality <- as.factor(burns$Locality)
burns$Locality <- droplevels(burns$Locality)

rfsburn <- glmer(Count_ray_flourets ~ diff + (1 | Locality), data = burns, family = poisson)
summary(rfsburn)
burns2 <- burns %>% drop_na(Count_disk_flourets) 
dfsburn <- glmer(Count_disk_flourets ~ diff + (1 | Locality), data = burns2, family = poisson)
summary(dfsburn)
burns3 <- burns %>% filter(diff>-250) 
rfsburn2 <- glmer(Count_ray_flourets ~ diff + (1 | Locality), data = burns3, family = poisson)
summary(rfsburn2)
burns4 <- burns2 %>% filter(diff>-250) 
dfsburn2 <- glmer(Count_disk_flourets ~ diff + (1 | Locality), data = burns2, family = poisson)
summary(dfsburn2)
#significant


ggplot(my_sf) + geom_sf()
ggplot(cropped) + geom_sf()


burns %>% 
  wilcox_test(Count_disk_flourets ~ burned, paired=F) %>%
  add_significance()

burns %>% 
  wilcox_test(Count_ray_flourets ~ burned, paired=F) %>%
  add_significance()
