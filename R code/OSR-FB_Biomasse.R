library(lme4)
library(ggplot2)
library(rstatix)
library(broom)
library(emmeans)
library('data.table')

Berlese <- read.csv("Data/Crop_characteristics_2023.csv")
head (Berlese)

Berlese <- subset(Berlese, variety !='4')
Berlese $bloc.y<- as.factor(Berlese $bloc.y)
Berlese $variety <- as.factor(Berlese $variety)

Berlese$dry_weight<- Berlese$dry_weight/5

mymodel <- lmer(dry_weight ~ variety*culture.y + (1|bloc.y), data = Berlese)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Berlese$variety)

Anova(mymodel, type=3)

emm<-emmeans(mymodel, ~ variety*culture.y, type = 'response')

emm<-as.data.frame(emm)

emm <- emm %>% rename(culture = culture.y)

emm$exp <- paste("Field 2022-2023")

plot <- ggplot(emm, aes(x=culture, y=emmean, fill = variety )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#E69F00','#999999','#69b3a2'))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) 

plot +labs(title="Plant dry biomass 2023", x="Culture type", y = "OSR plant dry biomass (g)")+
  theme_classic() 


Berlese <- read.csv("Data/Crop_characteristics_2024.csv", sep=";")
head(Berlese)

Berlese $bloc<- as.factor(Berlese $bloc)
Berlese $variety <- as.factor(Berlese $variety)

Berlese$dry_weight<- Berlese$poids_sec/3

Berlese <- subset(Berlese, Berlese$poids_sec < 100)

mymodel <- lmer(poids_sec ~ variety*culture+ (1|bloc), data = Berlese)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Berlese$variety)

Anova(mymodel, type=3)

emmeans(mymodel, list(pairwise ~ variety*culture), adjust = "tukey")

emm_2024<-emmeans(mymodel, ~ variety*culture, type = 'response')

emm_2024<-as.data.frame(emm_2024)

emm_2024$exp <- paste("Field 2023-2024")

plot <- ggplot(emm_2024, aes(x=culture, y=emmean, fill = variety )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#E69F00','#999999','#69b3a2'))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) 

plot +labs(title="Plant dry biomass 2024", x="Culture type", y = "OSR plant dry biomass (g)")+
  theme_classic() 

emmean <- rbind(emm, emm_2024)
emmean

#Figure for mambo 
emmean<-as.data.frame(emmean)
emmean
emmean_mambo <- subset(emmean, variety == "2")

plot_mambo <- ggplot(emmean_mambo, aes(x=exp, y=emmean, fill =culture  )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) 

plot_mambo <- plot_mambo +labs( x="Experiment", y = "EMMean of plant biomass (g)")+
  theme_classic() 

plot_mambo

write_rds(plot_mambo, "~/working/Aboveground_characterization/3-Output/Biomass_exp_mambo.rds")

#Figure for feliciano 
emmean<-as.data.frame(emmean)
emmean
emmean_feliciano <- subset(emmean, variety == "3")

plot_feliciano <- ggplot(emmean_feliciano, aes(x=exp, y=emmean, fill =culture  )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) 

plot_feliciano <- plot_feliciano +labs( x="Experiment", y = "EMMean of plant biomass (g)")+
  theme_classic() 

plot_feliciano

write_rds(plot_feliciano, "~/working/Aboveground_characterization/3-Output/Biomass_exp_feliciano.rds")

#Figure for angelico 

emmean_angelico <- subset(emmean, variety == "1")

plot_angelico <- ggplot(emmean_angelico, aes(x=exp, y=emmean, fill =culture  )) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#5773CCFF','#FFB900FF'))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) 

plot_angelico <- plot_angelico +labs( x="Experiment", y = "EMMean of plant biomass (g)")+
  theme_classic() 

plot_angelico

write_rds(plot_angelico, "~/working/Aboveground_characterization/3-Output/Biomass_exp_angelico.rds")
