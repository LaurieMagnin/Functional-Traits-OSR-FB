library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

Donnees_Physiologiques <- read.csv("Data/GlassHouse_ecophysiological_traits.csv", sep=";")

Donnees_Physiologiques<- filter(Donnees_Physiologiques, no_feuille !="cot1")
Donnees_Physiologiques<- filter(Donnees_Physiologiques, no_feuille !="cot2")
Donnees_Physiologiques<- filter(Donnees_Physiologiques, no_feuille !="na")


Donnees_Physiologiques$no_feuille <- as.numeric(Donnees_Physiologiques$no_feuille)

Donnees_Physiologiques$side <- ifelse(Donnees_Physiologiques$bloc == "nordest", "est", ifelse(Donnees_Physiologiques$bloc == "sudest","est", "ouest"))

Donnees_Physiologiques$variety <- ifelse(Donnees_Physiologiques$variete == "angelico", "1_angelico", ifelse(Donnees_Physiologiques$variete == "feliciano","3_feliciano","2_mambo"))


feuilles <- Donnees_Physiologiques %>%                                 # Group data
  group_by(variety, culture, bac, plante, date, side) %>%
  dplyr::summarize(nb_feuilles = max(no_feuille)) %>% 
  as.data.frame()

DegreeDay <- read.csv("~/working/Aboveground_characterization/1-Data/Date_DegreeDay_GlassHouse.csv", sep=";")

head(DegreeDay)

DegreeDay<- DegreeDay %>%
  rename(date = Date)

DegreeDay <- DegreeDay %>%
  mutate(date = as.Date(date, format = "%Y.%m.%d"))  # Convert `date` to Date object

DegreeDay$microclimate <- DegreeDay$DD-DegreeDay$DD*0.02

DegreeDay <- DegreeDay %>%
  mutate(Cum_microclimate = cumsum(microclimate))

feuilles <- feuilles %>%
  mutate(date = as.Date(date, format = "%Y.%m.%d"))  # Ensure `date` is also a Date object


library(dplyr)
head(feuilles)

feuilles <- feuilles %>%
  left_join(DegreeDay, by = "date")

head(feuilles)

feuilles_average <- feuilles %>%
  group_by(culture, Cum_DD, variety) %>%
  dplyr::summarize(total_av = mean(nb_feuilles),
                   n=n(), 
                   total_sd= sd(nb_feuilles), total_se = total_sd / sqrt(n))%>% 
  as.data.frame()

feuilles_average$Cum_DD<- as.numeric(feuilles_average$Cum_DD)
feuilles_average <-subset(feuilles_average, feuilles_average$Cum_DD < 1500)

## Dynamic of leaf number mambo ####

feuilles_average_mambo <-subset(feuilles_average, feuilles_average$variety == "2_mambo")

plot <- ggplot(feuilles_average_mambo, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  theme_classic()+
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20))

plot_serre <- plot + labs( x="Cumulated degree days (°C)", y = "Cumulated number of leaves")+
  theme_classic() 
plot_serre

write_rds(plot_serre, "~/working/Aboveground_characterization/3-Output/Leaves_GH_mambo.rds")


## Dynamic of leaf number feliciano ####

feuilles_average_feliciano <-subset(feuilles_average, feuilles_average$variety == "3_feliciano")

plot_feliciano <- ggplot(feuilles_average_feliciano, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  theme_classic()+
  scale_y_continuous(limits = c(0, 22), breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20))

plot_feliciano <- plot_feliciano + labs(title="Cumulated number of leaves during the glasshouse experiment", x="Cumulated degree days (°C)", y = "Cumulated number of leaves")+
  theme_classic() 

plot_feliciano

write_rds(plot_feliciano, "~/working/Aboveground_characterization/3-Output/Leaves_GH_feliciano.rds")

## Dynamic of leaf number angelico ####

feuilles_average_angelico <-subset(feuilles_average, feuilles_average$variety == "1_angelico")

plot_angelico <- ggplot(feuilles_average_angelico, aes(x=Cum_DD, y=total_av, color = culture )) +
  geom_point()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_line(aes(group = culture)) + 
  geom_errorbar(aes(ymax = total_av+total_se, ymin = total_av-total_se, width = .02))+ 
  theme_classic()+
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20))

plot_serre_angelico <- plot_angelico + labs(title="Cumulated number of leaves during the glasshouse experiment", x="Cumulated degree days (°C)", y = "Cumulated number of leaves")+
  theme_classic() 

plot_serre_angelico

write_rds(plot_serre_angelico, "~/working/Aboveground_characterization/3-Output/Leaves_GH_angelico.rds")
