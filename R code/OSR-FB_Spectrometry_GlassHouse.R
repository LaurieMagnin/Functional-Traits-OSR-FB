library(dplyr)
library(ggplot2)
library(tidyr)

data_PAR <- read.csv("Data/GlassHouse_light_PAR.csv", sep = ";")
head(data_PAR)

##### Producing the figure for PAR ######

data_PAR$PAR_down <- as.numeric(data_PAR$PAR_down)
data_PAR$PAR_percent <- as.numeric(data_PAR$PAR_percent )

data_PAR <- data_PAR %>%                                 # Group data
  group_by(variety, culture, bac, distance, repetition) %>%
  dplyr::summarize( PAR_top = mean(PAR_140), PAR_down = mean(PAR_down), PAR_percent = mean(PAR_percent) ) %>% 
  as.data.frame()

data_PAR_AV <- data_PAR %>%                                 # Group data
  group_by( culture, distance) %>%
  dplyr::summarize(PAR_percent_AV = mean(PAR_percent), n=n(), PAR_percent_sd= sd(PAR_percent), 
                   PAR_percent_se =PAR_percent_sd / sqrt(n) ) %>% 
  as.data.frame()

data_PAR_AV <-  filter(data_PAR_AV , data_PAR_AV $distance != '125')
data_PAR_AV <-  filter(data_PAR_AV , data_PAR_AV $distance != '120')
data_PAR_AV <-  filter(data_PAR_AV , data_PAR_AV $distance != 'NA')

data_PAR_AV $distance <- as.factor (data_PAR_AV $distance)
data_PAR_AV$PAR_corr <- data_PAR_AV$PAR_percent_AV*5
data_PAR_AV$PAR_corr_se <- data_PAR_AV$PAR_percent_se*5
head(data_PAR_AV)
 
PAR_GH <- ggplot(data_PAR_AV, aes(PAR_corr, distance, colour = culture)) + 
  geom_point()+
  geom_errorbarh(aes(xmax =PAR_corr + PAR_corr_se, xmin = PAR_corr - PAR_corr_se, height = .2))+ 
  theme_bw()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800))+
  labs(title="Glasshouse",x="Corrected transmitted PAR", y = "Height under canopy (cm)")+  theme( legend.position="none")

PAR_GH

write_rds(PAR_GH, "~/working/Aboveground_characterization/3-Output/PAR_GH.rds")

#### Stat analysis 

data_PAR_underOSR <- subset(data_PAR, data_PAR$distance == "5")
head(data_PAR_underOSR)

kruskal.test(PAR_percent~  culture, data = data_PAR_underOSR)

# Stat analysis under OSR canopy
data_PAR$dicu <- paste(data_PAR$distance, data_PAR$culture)
data_PAR_aboveOSR <- subset(data_PAR, data_PAR$dicu == "65 Mo"| data_PAR$dicu == "55 As" )
head(data_PAR_underOSR)

kruskal.test(PAR_percent~  culture, data = data_PAR_aboveOSR)

###################################################################################################################################################################
##### R:FR analysis ######

spectro <- read.csv("Data/2022.06.09_light_quality_under_canopy.csv", sep = ";")

# Reshape data 

spectro <- spectro %>% gather("wave_length", "irradiance", 7:477)

# Averaging the three measures for each point 

spectro <- spectro %>%                                 # Group data
  group_by(variety, culture, bac, location, hight, wave_length) %>%
  dplyr::summarize( av_irradiance = mean(irradiance)) %>% 
  as.data.frame()

# R:FR ratio calculation 

spectro<- spectro %>% separate_wider_delim(wave_length, "X", names = c(NA, "Wl")) #Removing the X before wavelength

spectro$Wl <- as.numeric(spectro$Wl)

spectro$photonflux <- (spectro$Wl * spectro$av_irradiance * 0.001) / 0.11962657 #Converting from energy to photon flux

spectro <- select(spectro, -av_irradiance)

spectro<-spread(spectro, Wl, photonflux)

spectro$red <- rowSums(spectro[ , c(246:346)], na.rm=TRUE)

spectro$Fred <- rowSums(spectro[ , c(346:446)], na.rm=TRUE)

spectro$R_FR <- spectro$red/spectro$Fred
head(spectro)

spectro$PAR <- rowSums(spectro[ , c(46:346)], na.rm=TRUE)

spectro$BLUE <- rowSums(spectro[ , c(46:146)], na.rm=TRUE)

spectro$R_FR <- as.numeric(spectro$R_FR)

spectro$R_FR <- as.numeric(spectro$R_FR)
spectro$hight <- as.factor(spectro$hight)

## Data from the 06.07.2023 ####

spectro_1 <- read.csv("Data/2022.06.07_Light_quality_under_canopy.csv", sep = ";")

# Reshape data 
library(tidyr)

spectro_1  <- spectro_1  %>% gather("wave_length", "irradiance", 7:477)

# R:FR ratio calculation 

spectro_1<- spectro_1 %>% separate_wider_delim(wave_length, "X", names = c(NA, "Wl")) #Removing the X before wavelength

spectro_1$av_irradiance<-spectro_1$irradiance
spectro_1 <- select(spectro_1, -irradiance)

spectro_1 $Wl <- as.numeric(spectro_1 $Wl)

spectro_1 $photonflux <- (spectro_1 $Wl * spectro_1 $av_irradiance * 0.001) / 0.11962657 #Converting from energy to photon flux

spectro_1  <- select(spectro_1 , -av_irradiance)

spectro_1 <-spread(spectro_1 , Wl, photonflux)

spectro_1$red <- rowSums(spectro_1[ , c(246:346)], na.rm=TRUE)

spectro_1$Fred <- rowSums(spectro_1[ , c(346:446)], na.rm=TRUE)

spectro_1$R_FR <- spectro_1$red/spectro_1$Fred
head(spectro_1)

spectro_1$PAR <- rowSums(spectro_1[ , c(46:346)], na.rm=TRUE)
spectro_1$BLUE <- rowSums(spectro_1[ , c(46:146)], na.rm=TRUE)

spectro_1$PAR <- as.numeric(spectro_1$PAR)
spectro_1$R_FR <- as.numeric(spectro_1$R_FR)
spectro_1$hight <- as.factor(spectro_1$hight)

## Merge data and results ####
spectro_1 <- select(spectro_1, -file)

spectro_total <- rbind(spectro, spectro_1)

# Averaging the measures for each levels in the canopy

#For the first set of data (07.06.2022)
spectro_1_AV <- spectro_1 %>%                                 # Group data
  group_by(culture, hight) %>%
  dplyr::summarize( av_R_FR = mean(R_FR), n=n(),R_FR_sd= sd(R_FR), 
                    R_FR_se =R_FR_sd / sqrt(n)) %>% 
  as.data.frame()

#For the second set of data (09.06.2022)

spectro_AV <- spectro %>%                                 # Group data
  group_by(culture, hight) %>%
  dplyr::summarize( av_R_FR = mean(R_FR), n=n(),R_FR_sd= sd(R_FR), 
                    R_FR_se =R_FR_sd / sqrt(n)) %>% 
  as.data.frame()

# Combining both of the sets

spectro_total_AV <- spectro_total %>%                                 # Group data
  group_by(culture, hight) %>%
  dplyr::summarize( av_R_FR = mean(R_FR), n=n(),R_FR_sd= sd(R_FR), 
                    R_FR_se =R_FR_sd / sqrt(n)) %>% 
  as.data.frame()

spectro_total_AV <-  filter(spectro_total_AV , spectro_total_AV$hight != '125')
spectro_total_AV <-  filter(spectro_total_AV , spectro_total_AV$hight != '135')

R.FR_GH <- ggplot(spectro_total_AV, aes(av_R_FR,hight, colour = culture)) + 
  geom_point()+
  geom_errorbarh(aes(xmax = av_R_FR + R_FR_se, xmin = av_R_FR - R_FR_se, height = .2))+ 
  theme_bw()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3))+
  labs(title="Glasshouse",x="Ratio R/FR", y = "Height under canopy (cm)")+  theme( legend.position="none")

R.FR_GH

write_rds(R.FR_GH, "~/working/Aboveground_characterization/3-Output/R.FR_GH.rds")

#### Stat analysis 

spectro_total_underOSR <- subset(spectro_total, spectro_total$hight == "5")
head(spectro_total_underOSR)

kruskal.test(R_FR~  culture, data = spectro_total_underOSR)

# Stat analysis under OSR canopy
spectro_total$dicu <- paste(spectro_total$hight, spectro_total$culture)
spectro_total_aboveOSR <- subset(spectro_total, spectro_total$dicu == "65 Mo"| spectro_total$dicu == "55 As" )
head(spectro_total_aboveOSR)

kruskal.test(R_FR~  culture, data = spectro_total_aboveOSR)
