strauss<-gt[gt$Locality=="Strauss"  & gt$year>2020 & gt$Lon_dd.dddd < -120.475,]
strauss$period <- ifelse(strauss$year < 2024, "2020-2023", "2024")
strauss$period <- factor(strauss$period, levels = c("2020-2023", "2024"))


ggplot(strauss, aes(x = period, y = Count_ray_flourets)) +
  geom_boxplot(fill = "gray80", color = "black", outlier.size = 2) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  labs(x = "", y = "Ray floret count") +
  theme_cowplot(font_size=20, font_family='Palatino') +
  theme(panel.grid.major=element_line(linetype=1, color='gray90'),
        panel.grid.minor.y=element_line(linetype=1, color='gray90'))
ggsave("ray_before_after_2023_font.png", height=1200, width=1200, bg='white',units="px", path="Figures/phenotype/turbines")


        