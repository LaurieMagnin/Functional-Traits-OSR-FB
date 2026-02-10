
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

#Dataframing Glasshouse (GH) stem elongation data 

Notations_structures_variete <- read.csv("Data/GlassHouse_ecophysiological_traits.csv", sep=";")

Notations_structures_variete$variete<-as.factor(Notations_structures_variete$variete)
Notations_structures_variete$culture<-as.factor(Notations_structures_variete$culture)
Notations_structures_variete$H_total<-as.numeric(Notations_structures_variete$H_total)
Notations_structures_variete<- subset(Notations_structures_variete, Notations_structures_variete$no_feuille =="1")

DegreeDay <- read.csv("Data/Date_DegreeDay_GlassHouse.csv", sep=";")

head(DegreeDay)

DegreeDay<- DegreeDay %>%
  rename(date = Date)

DegreeDay <- DegreeDay %>%
  mutate(date = as.Date(date, format = "%Y.%m.%d"))  # Convert `date` to Date object

Notations_structures_variete <- Notations_structures_variete %>%
  mutate(date = as.Date(date, format = "%Y.%m.%d"))  # Ensure `date` is also a Date object

library(dplyr)
head(Notations_structures_variete)

Notations_structures_variete <- Notations_structures_variete %>%
  left_join(DegreeDay, by = "date")


elongation_average <- Notations_structures_variete %>%
  group_by(culture, variete, Cum_DD) %>%
  dplyr::summarize(total_av = mean(H_total),
                   n=n(), 
                   total_sd= sd(H_total), total_se = total_sd / sqrt(n))%>% 
  as.data.frame()

elongation_average$Cum_DD<- as.numeric(elongation_average$Cum_DD)
#elongation_average <-subset(elongation_average, elongation_average$Cum_DD < 1500)

# Stem Elongation dynamic for mambo in GH ####

elongation_average_mambo <- subset(elongation_average, elongation_average$variete == "mambo")

plot_mambo <- ggplot(elongation_average_mambo, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(700, 2100), breaks = c(750, 1000, 1250, 1500, 1750, 2000))+
  theme_classic()

plot_mambo_serre <-  plot_mambo + labs(title="GlassHouse", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() 

plot_mambo_serre 

write_rds(plot_bites_serre, "~/working/Aboveground_characterization/3-Output/Stem_GH_mambo.rds")

# Stem Elongation dynamic for feliciano in GH ####

elongation_average_feliciano <- subset(elongation_average, elongation_average$variete == "feliciano")

plot_feliciano <- ggplot(elongation_average_feliciano, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(700, 2100), breaks = c(750, 1000, 1250, 1500, 1750, 2000))+
  theme_classic()

plot_feliciano_serre <-  plot_feliciano + labs(title="GlassHouse", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() 

plot_feliciano_serre 

write_rds(plot_feliciano_serre, "~/working/Aboveground_characterization/3-Output/Stem_GH_feliciano.rds")

# Stem Elongation dynamic for angelico in GH ####

elongation_average_angelico <- subset(elongation_average, elongation_average$variete == "angelico")

plot_angelico <- ggplot(elongation_average_angelico, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(700, 2100), breaks = c(750, 1000, 1250, 1500, 1750, 2000))+
  theme_classic()

plot_angelico_serre <-  plot_angelico + labs(title="GlassHouse", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() 

plot_angelico_serre 

write_rds(plot_bites_serre, "~/working/Aboveground_characterization/3-Output/Stem_GH_angelico.rds")


##Dataframing for the 2022-2023 field experiments

Notations_structures_variete <- read.csv("Data/Field_ecophysiological_traits.csv", sep=";")
head(Notations_structures_variete)

Notations_structures_variete$h_tige<-as.numeric(Notations_structures_variete$h_tige)

DegreeDay <- read.csv("Data/Climate_Changins.csv", sep=";")

head(DegreeDay)

DegreeDay<- DegreeDay %>%
  rename(date = Day)

DegreeDay <- DegreeDay %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y"))  # Convert `date` to Date object

Notations_structures_variete <- Notations_structures_variete %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y"))  # Ensure `date` is also a Date object

library(dplyr)
head(Notations_structures_variete)

Notations_structures_variete <- Notations_structures_variete %>%
  left_join(DegreeDay, by = "date")

Notations_structures_variete_2023 <- subset(Notations_structures_variete, Notations_structures_variete$date < "2023-08-01")

elongation_average_2023 <- Notations_structures_variete_2023 %>%
  group_by(culture, Cum_DD, variete) %>%
  dplyr::summarize(total_av = mean(h_tige),
                   n=n(), 
                   total_sd= sd(h_tige), total_se = total_sd / sqrt(n))%>% 
  as.data.frame()

elongation_average_2023 <- subset(elongation_average_2023, elongation_average_2023$total_av != "NA" )
elongation_average_2023 <- subset(elongation_average_2023, elongation_average_2023$Cum_DD != "1665.05")

# Stem Elongation dynamic for mambo in 2022-2023 Field ####
elongation_average_2023_mambo <- subset(elongation_average_2023, elongation_average_2023$variete == "2" )

plot_mambo_2023 <- ggplot(elongation_average_2023_mambo, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(1500, 2100), breaks = c( 1500, 1750, 2000))+ 
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 20, 40,60))+
  theme_classic()

plot_mambo_2023 <-plot_mambo_2023 + labs(title="Field 2022-2023", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() +  theme( legend.position="none")

plot_mambo_2023

write_rds(plot_mambo_2023, "~/working/Aboveground_characterization/3-Output/Stem_F1_mambo.rds")

# Stem Elongation dynamic for feliciano in 2022-2023 Field ####
elongation_average_2023_feliciano <- subset(elongation_average_2023, elongation_average_2023$variete == "3" )

plot_feliciano_2023 <- ggplot(elongation_average_2023_feliciano, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(1500, 2100), breaks = c( 1500, 1750, 2000))+ 
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 20, 40,60))+
  theme_classic()

plot_feliciano_2023 <-plot_feliciano_2023 + labs(title="Field 2022-2023", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() +  theme( legend.position="none")

plot_feliciano_2023

write_rds(plot_feliciano_2023, "~/working/Aboveground_characterization/3-Output/Stem_F1_feliciano.rds")

# Stem Elongation dynamic for angelico in 2022-2023 Field ####
elongation_average_2023_angelico <- subset(elongation_average_2023, elongation_average_2023$variete == "1" )

plot_angelico_2023 <- ggplot(elongation_average_2023_angelico, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(1500, 2100), breaks = c( 1500, 1750, 2000))+ 
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 20, 40,60))+
  theme_classic()

plot_angelico_2023 <-plot_angelico_2023 + labs(title="Field 2022-2023", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() +  theme( legend.position="none")

plot_angelico_2023

write_rds(plot_angelico_2023, "~/working/Aboveground_characterization/3-Output/Stem_F1_angelico.rds")

# Data framing for the 2023-2024 field experiment

Notations_structures_variete_2024 <- subset(Notations_structures_variete, Notations_structures_variete$date > "2023-08-01")

elongation_average_2024 <- Notations_structures_variete_2024 %>%
  group_by(culture, Cum_DD, variete) %>%
  dplyr::summarize(total_av = mean(h_tige),
                   n=n(), 
                   total_sd= sd(h_tige), total_se = total_sd / sqrt(n))%>% 
  as.data.frame()

elongation_average_2024 <- subset(elongation_average_2024, elongation_average_2024$total_av != "NA" )
elongation_average_2024 <- subset(elongation_average_2024, elongation_average_2024$Cum_DD != 1867.90 )
elongation_average_2024 <- subset(elongation_average_2024, elongation_average_2024$Cum_DD != "1983.55" )

# Stem Elongation dynamic for mambo in 2023-2024 Field ####
elongation_average_2024_mambo <- subset(elongation_average_2024, elongation_average_2024$variete == "2" )

plot_mambo_2024 <- ggplot(elongation_average_2024_mambo, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(1400, 2350), breaks = c( 1500, 1750, 2000, 2250))+ 
  scale_y_continuous(limits = c(0, 130), breaks = c(0, 20, 40, 60, 80, 100, 120))+
  theme_classic()

plot_mambo_2024 <-plot_mambo_2024 + labs(title="Field 2023-2024", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() +  theme( legend.position="none")

plot_mambo_2024

write_rds(plot_mambo_2024, "~/working/Aboveground_characterization/3-Output/Stem_F2_mambo.rds")

# Stem Elongation dynamic for feliciano in 2023-2024 Field ####
elongation_average_2024_feliciano <- subset(elongation_average_2024, elongation_average_2024$variete == "3" )

plot_feliciano_2024 <- ggplot(elongation_average_2024_feliciano, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(1400, 2350), breaks = c( 1500, 1750, 2000, 2250))+ 
  scale_y_continuous(limits = c(0, 140), breaks = c(0, 20, 40, 60, 80, 100, 120))+
  theme_classic()

plot_feliciano_2024 <-plot_feliciano_2024 + labs(title="Field 2023-2024", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() +  theme( legend.position="none")

plot_feliciano_2024

write_rds(plot_feliciano_2024, "~/working/Aboveground_characterization/3-Output/Stem_F2_feliciano.rds")

# Stem Elongation dynamic for angelico in 2023-2024 Field ####
elongation_average_2024_angelico <- subset(elongation_average_2024, elongation_average_2024$variete == "1" )

plot_angelico_2024 <- ggplot(elongation_average_2024_angelico, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  scale_x_continuous(limits = c(1400, 2350), breaks = c( 1500, 1750, 2000, 2250))+ 
  scale_y_continuous(limits = c(0, 140), breaks = c(0, 20, 40, 60, 80, 100, 120))+
  theme_classic()

plot_angelico_2024 <-plot_angelico_2024 + labs(title="Field 2023-2024", x="Cumulated degree days (°C)", y = "Stem elongation (cm)")+
  theme_classic() +  theme( legend.position="none")
plot_angelico_2024

write_rds(plot_angelico_2024, "~/working/Aboveground_characterization/3-Output/Stem_F2_angelico.rds")
