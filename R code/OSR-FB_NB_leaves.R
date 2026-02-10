##Analysis fot the glasshouse experiment all varieties ####
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(glmmTMB)

Donnees_Physiologiques <- read.csv("Data/GlassHouse_ecophysiological_traits.csv", sep=";")

Donnees_Physiologiques<- filter(Donnees_Physiologiques, no_feuille !="cot1")
Donnees_Physiologiques<- filter(Donnees_Physiologiques, no_feuille !="cot2")
Donnees_Physiologiques<- filter(Donnees_Physiologiques, no_feuille !="na")


Donnees_Physiologiques$no_feuille <- as.numeric(Donnees_Physiologiques$no_feuille)

Donnees_Physiologiques$side <- ifelse(Donnees_Physiologiques$bloc == "nordest", "est", ifelse(Donnees_Physiologiques$bloc == "sudest","est", "ouest"))

Donnees_Physiologiques$variety <- ifelse(Donnees_Physiologiques$variete == "angelico", "1", ifelse(Donnees_Physiologiques$variete == "feliciano","3","2"))


feuilles <- Donnees_Physiologiques %>%                                 # Group data
  group_by(variety, culture, bac, plante, date, side) %>%
  dplyr::summarize(nb_feuilles = max(no_feuille)) %>% 
  as.data.frame()


Feuilles_BFwinter <- subset(feuilles, date == " 2022.04.03")

hist(Feuilles_BFwinter$nb_feuilles, col='steelblue', main='Original')

head(Feuilles_BFwinter)

Feuilles_BFwinter$pot <- paste(Feuilles_BFwinter$variety, Feuilles_BFwinter$culture, Feuilles_BFwinter$bac)
head(Feuilles_BFwinter)

mymodel<- glmmTMB(nb_feuilles ~ variety*culture + (1|pot), data= Feuilles_BFwinter, family =compois)
summary(mymodel)

library(DHARMa)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Notations_structures_variete$variete)

library(rstatix)
Anova(mymodel, type=3)

library(emmeans)

emmeans(mymodel, list(pairwise ~ variety), adjust = "tukey")

emm<-emmeans(mymodel, ~ variety*culture, type = 'response')

library('data.table')

emm<- as.data.table(emm)
emm

emm$exp <- paste("Glasshouse")

plot<- ggplot(emm, aes(x=culture, y=response, fill = variety )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#E69F00','#999999','#69b3a2'))+
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8))

plot +labs(title="Number of leaves before wintering (Glasshouse) ", x="Modality", y = "EMMean number of leaves")+
  theme_classic() 

##Analysis for the field experiments

Notations_structures_variete <- read.csv("Data/Field_ecophysiological_traits.csv", sep=";")

Notations_structures_variete$variete<-as.factor(Notations_structures_variete$variete)
Notations_structures_variete$bloc<-as.factor(Notations_structures_variete$bloc)
Notations_structures_variete$nb_feuille<-as.numeric(Notations_structures_variete$nb_feuille)

##### Analysis for the number of leave in 2022

Notations_structures_variete <-subset(Notations_structures_variete, date=='02.11.2022')
Notations_structures_variete <-subset(Notations_structures_variete, nb_feuille!=19)
Notations_structures_variete <-subset(Notations_structures_variete, variete!=4)

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(glmmTMB)

ggplot(Notations_structures_variete , aes(x=culture, y=nb_feuille, fill=variete))  + geom_boxplot()

ggplot(Notations_structures_variete , aes(x=variete, y=nb_feuille))  + geom_boxplot()

ggplot(Notations_structures_variete , aes(x=culture, y=nb_feuille))  + geom_boxplot()

ggplot(Notations_structures_variete , aes(x=bloc, y=nb_feuille))  + geom_boxplot()

hist(Notations_structures_variete$nb_feuille, col='steelblue', main='Original')

head(Notations_structures_variete)

mymodel<- glmmTMB(nb_feuille ~ variete*culture + (1|bloc), data= Notations_structures_variete, family =compois)
summary(mymodel)

library(DHARMa)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Notations_structures_variete$variete)

library(rstatix)
Anova(mymodel, type=3)

library(emmeans)

emm_2023<-emmeans(mymodel, ~ variete*culture, type = 'response')

library('data.table')

emm_2023<- as.data.table(emm_2023)

emm_2023

emm_2023 <- emm_2023 %>% rename(variety = variete)

emm_2023$exp <- paste("Field 2022-2023")

plot <- ggplot(emm_2023, aes(x=culture, y=response, fill = variety )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#E69F00','#999999','#69b3a2','seagreen'))+
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8))

plot +labs(title="Number of green leaves before winter (02.11.2022) N=45", x="Modality", y = "EMMean number of leaves")+
  theme_classic() 


##### Analysis for the number of leave in  end of september 2023

Notations_structures_variete <- read.csv("Data/Field_ecophysiological_traits.csv", sep=";")

Notations_structures_variete$variete<-as.factor(Notations_structures_variete$variete)
Notations_structures_variete$bloc<-as.factor(Notations_structures_variete$bloc)
Notations_structures_variete$nb_feuille<-as.numeric(Notations_structures_variete$nb_feuille)

Notations_structures_variete <-subset(Notations_structures_variete, date=='26.09.2023')

Notations_structures_variete$select_plant <- ifelse(Notations_structures_variete$plante < 6, "yes",
                                                    ifelse(Notations_structures_variete$plante < 9, "no", 
                                                           ifelse(Notations_structures_variete$plante < 14, "yes", 
                                                                  ifelse(Notations_structures_variete$plante < 20, "no", "yes"))))


Notations_structures_variete <-subset(Notations_structures_variete, Notations_structures_variete$select_plant == "yes")

hist(Notations_structures_variete$nb_feuille, col='steelblue', main='Original')

head(Notations_structures_variete)
library(glmmTMB())

mymodel<- glmmTMB(nb_feuille ~ variete*culture + (1|bloc), data= Notations_structures_variete, family =compois)
summary(mymodel)

library(DHARMa)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Notations_structures_variete$variete)

library(rstatix)
Anova(mymodel, type=3)

library(emmeans)

emm_2024<-emmeans(mymodel, ~ variete*culture, type = 'response')

library('data.table')

emm_2024<- as.data.table(emm_2024)

emm_2024 <- emm_2024 %>% rename(variety = variete)

emm_2024$exp <- paste("Field 2023-2024")

plot <- ggplot(emm_2024, aes(x=culture, y=response, fill = variety )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#E69F00','#999999','#69b3a2','seagreen'))+
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8))

plot +labs(title="Number of leaves at P. chrysocephala migration", x="Modality", y = "EMMean number of leaves")+
  theme_classic() 


##### Analysis for the number of leave in  end of november 2023

Notations_structures_variete <- read.csv("Data/Field_ecophysiological_traits.csv", sep=";")

Notations_structures_variete$variete<-as.factor(Notations_structures_variete$variete)
Notations_structures_variete$bloc<-as.factor(Notations_structures_variete$bloc)
Notations_structures_variete$nb_feuille<-as.numeric(Notations_structures_variete$nb_feuille)

Notations_structures_variete <-subset(Notations_structures_variete, date=='21.11.2023')

Notations_structures_variete$select_plant <- ifelse(Notations_structures_variete$plante < 6, "yes",
                                                    ifelse(Notations_structures_variete$plante < 9, "no", 
                                                           ifelse(Notations_structures_variete$plante < 14, "yes", 
                                                                  ifelse(Notations_structures_variete$plante < 20, "no", "yes"))))
                                                    

Notations_structures_variete <-subset(Notations_structures_variete, Notations_structures_variete$select_plant == "yes")

hist(Notations_structures_variete$nb_feuille, col='steelblue', main='Original')

head(Notations_structures_variete)

mymodel<- glmmTMB(nb_feuille ~ variete*culture + (1|bloc), data= Notations_structures_variete, family =compois)
summary(mymodel)

library(DHARMa)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Notations_structures_variete$variete)

library(rstatix)
Anova(mymodel, type=3)

library(emmeans)

emm_2024<-emmeans(mymodel, ~ variete*culture, type = 'response')

library('data.table')

emm_2024<- as.data.table(emm_2024)

emm_2024 <- emm_2024 %>% rename(variety = variete)

emm_2024$exp <- paste("Field 2023-2024")

emmeans(mymodel, list(pairwise ~ variete), adjust = "tukey")

plot <- ggplot(emm_2024 , aes(x=culture, y=response, fill = variety )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#E69F00','#999999','#69b3a2','seagreen'))+
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8))

plot +labs(title="Number of leaves at P. chrysocephala migration", x="Modality", y = "EMMean number of leaves")+
  theme_classic() 

### Figure by variety ####

emmean <- rbind(emm, emm_2023,emm_2024)

#Figure for mambo 
emmean_mambo <- subset(emmean, emmean$variety =="2")

emmean_mambo <- emmean_mambo %>%
  mutate(exp = factor(exp, levels = c("Glasshouse", "Field 2022-2023", "Field 2023-2024")))
  
plot_mambo <- ggplot(emmean_mambo, aes(x=exp, y=response, fill =culture  )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9))+ 
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8))

plot_mambo <- plot_bites +labs( x="Experiment", y = "EMMean number of leaves")+
  theme_classic() 
plot_mambo
write_rds(plot_mambo, "~/working/Aboveground_characterization/3-Output/Leaves_exp_mambo.rds")

#Figure for feliciano
emmean_feliciano <- subset(emmean, emmean$variety =="3")

emmean_feliciano <- emmean_feliciano %>%
  mutate(exp = factor(exp, levels = c("Glasshouse", "Field 2022-2023", "Field 2023-2024")))

plot_feliciano <- ggplot(emmean_feliciano, aes(x=exp, y=response, fill =culture  )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8))

plot_feliciano <- plot_bites +labs( x="Experiment", y = "EMMean number of leaves")+
  theme_classic() 

plot_bites

write_rds(plot_feliciano, "~/working/Aboveground_characterization/3-Output/Leaves_exp_feliciano.rds")

#Figure for angelico
emmean_angelico <- subset(emmean, emmean$variety =="1")

emmean_angelico <- emmean_angelico %>%
  mutate(exp = factor(exp, levels = c("Glasshouse", "Field 2022-2023", "Field 2023-2024")))

plot_angelico <- ggplot(emmean_angelico, aes(x=exp, y=response, fill =culture  )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8))

plot_angelico <- plot_bites +labs( x="Experiment", y = "EMMean number of leaves")+
  theme_classic() 

plot_angelico

write_rds(plot_angelico, "~/working/Aboveground_characterization/3-Output/Leaves_exp_angelico.rds")



