spectro_2023 <- read.csv("Data/Field_light_spectrum.csv", sep = ";")

library(tidyr)

spectro_2023 <- spectro_2023 %>% gather("wave_length", "irradiance", 6:476)

head(spectro_2023)

library(dplyr)

spectro_2023  <- spectro_2023 %>%                                 # Group data
  group_by(bloc, culture,  location, location_culture, wave_length) %>%
  dplyr::summarize( av_irradiance = mean(irradiance)) %>% 
  as.data.frame()

spectro_2023 <- spectro_2023 %>% separate_wider_delim(wave_length, "X", names = c(NA, "Wl")) #Removing the X before wavelength

spectro_2023$Wl <- as.numeric(spectro_2023$Wl)

spectro_2023$photonflux <- (spectro_2023$Wl * spectro_2023$av_irradiance * 0.001) / 0.11962657 #Converting from energy to photon flux

spectro_2023 <- select(spectro_2023, -av_irradiance)

spectro_2023<-spread(spectro_2023, Wl, photonflux)

library(ggplot2)

spectro_2023$sensor_location <- paste(spectro_2023$location, spectro_2023$location_culture, sep="_") 

##### R/FR ratio analysis  #####

spectro_2023$red <- rowSums(spectro_2023[ , c(245:345)], na.rm=TRUE) # Calculating Red

spectro_2023$Fred <- rowSums(spectro_2023[ , c(345:445)], na.rm=TRUE) # Calculating Far Red

spectro_2023$R_FR <- spectro_2023$red/spectro_2023$Fred # Calculating R:FR ratio
head(spectro_2023)

spectro_2023$R_FR <- as.numeric(spectro_2023$R_FR)


spectro_2023_AV  <- spectro_2023  %>%                                 # Group data
  group_by( culture,sensor_location ) %>%
  dplyr::summarize( av_R_FR = mean(R_FR), n=n(),R_FR_sd= sd(R_FR), 
                    R_FR_se =R_FR_sd / sqrt(n)) %>% 
  as.data.frame()

spectro_2023_AV$height <-ifelse(spectro_2023_AV$sensor_location =="under_OSR", "5", 
                                ifelse(spectro_2023_AV$sensor_location =="above_OSR", "40", "86"))

spectro_2023_AV$height <- as.numeric(spectro_2023_AV$height)

R.FR_GH <- ggplot(spectro_2023_AV, aes(av_R_FR,height, colour = culture)) + 
  geom_point()+
  geom_errorbarh(aes(xmax = av_R_FR + R_FR_se, xmin = av_R_FR - R_FR_se, height = .2))+ 
  theme_bw()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  scale_y_continuous(breaks = c(5, 40, 86))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6,0.8,1))+
  labs(title="Field",x="Ratio R/FR", y = "Height under canopy (cm)")+  theme( legend.position="none")

R.FR_GH
library(readr)
write_rds(R.FR_GH, "~/working/Aboveground_characterization/3-Output/R.FR_Field.rds")

# Stat analysis under OSR canopy
spectro_2023_aboveOSR <- subset(spectro_2023, spectro_2023$location == "above")
spectro_2023_aboveOSR <- subset(spectro_2023_aboveOSR, spectro_2023_aboveOSR$location_culture == "OSR")
head(spectro_2023_aboveOSR)

kruskal.test(R_FR ~  culture, data = spectro_2023_aboveOSR)

#### PAR Analysis #####

spectro_2023$PAR <- rowSums(spectro_2023[ , c(46:346)], na.rm=TRUE)

spectro_2023$height <-ifelse(spectro_2023$sensor_location =="under_OSR", "5", 
                             ifelse(spectro_2023$sensor_location =="above_OSR", "40", "86"))

spectro_2023_figure <- spectro_2023  %>%                                 # Group data
  group_by(culture, height) %>%
  dplyr::summarize(av_PAR_corr = mean(PAR), n=n(), PAR_corr_sd= sd(PAR), 
                   PAR_corr_se =PAR_corr_sd / sqrt(n)) %>% 
  as.data.frame()

spectro_2023_figure $height <- as.numeric(spectro_2023_figure $height)

R.FR_GH <- ggplot(spectro_2023_figure, aes(av_PAR_corr,height, colour = culture)) + 
  geom_point()+
  geom_errorbarh(aes(xmax = av_PAR_corr + PAR_corr_se, xmin = av_PAR_corr - PAR_corr_se, height = .2))+ 
  theme_bw()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  scale_y_continuous(breaks = c(5, 40, 86))+
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800))+
  labs(title="Field",x="Corrected transmitted PAR", y = "Height under canopy (cm)")+  theme( legend.position="none")

R.FR_GH

library(readr)
write_rds(R.FR_GH, "~/working/Aboveground_characterization/3-Output/PAR_Field.rds")

# Stat analysis under OSR canopy
spectro_2023_aboveOSR <- subset(spectro_2023, spectro_2023$location == "above")
spectro_2023_aboveOSR <- subset(spectro_2023_aboveOSR, spectro_2023_aboveOSR$location_culture == "OSR")
head(spectro_2023_aboveOSR)

kruskal.test(PAR ~  culture, data = spectro_2023_aboveOSR)
