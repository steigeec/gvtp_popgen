subset <- gt %>%   filter(recordNumber_CollectorNumber > 8480 &  recordNumber_CollectorNumber < 8506)
subset <- subset %>%  filter(!is.na(subset$Count_disk_flourets))
subdisk <- mean(subset$Count_disk_flourets)
subray <- mean(subset$Count_ray_flourets)
gt2 <- subset(gt, !(recordNumber_CollectorNumber %in% c(8481, 8482, 8483, 8485, 8486, 8487, 8488, 8489, 8490, 
                                                        8491, 8492, 8493, 8494, 8495, 
                                                        8496, 8497, 8498, 8499, 8500, 8501, 8502, 8503, 8504, 
                                                        8505
                                                        )))
gt2 <- gt2 %>%  filter(!is.na(gt2$Count_disk_flourets))
gtdisk <- mean(gt2$Count_disk_flourets)
gtray <- mean(gt2$Count_ray_flourets)              

##### ratios?? ####

gt$ratio <- with(gt, Count_ray_flourets/Count_disk_flourets)
