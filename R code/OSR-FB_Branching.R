library(glmmTMB)
library(emmeans)
library(ggplot2)
library(dplyr)
library(DHARMa)
library(rstatix)
library('data.table')

##### Glass house #####

Donnees_physiologiques <- read.csv("Data/GlassHouse_ecophysiological_traits.csv", sep=";")
head(Donnees_physiologiques)

Donnees_ramifications <- subset(Donnees_physiologiques,
                                  Donnees_physiologiques$date == " 2022.07.12")

Donnees_ramifications$Ramification <- as.numeric(Donnees_ramifications$Ramification)

ramifications <- Donnees_ramifications %>%                                 # Group data
  group_by(variete, culture, bac, plante, date, bloc) %>%
  dplyr::summarize(nb_ramif = max(Ramification)) %>% 
  as.data.frame()

hist(ramifications$nb_ramif)
ramifications_asso <-subset(ramifications, culture=='asso')
hist(ramifications_asso$nb_ramif)
ramifications_mono <-subset(ramifications, culture=='mono')
hist(ramifications_mono$nb_ramif)

ramifications <- subset(ramifications, ramifications$nb_ramif != "NA")

ramifications$pot <- paste(ramifications$variete, ramifications$culture, ramifications$bac)
head(ramifications)

mymodel<- glmmTMB(nb_ramif ~ variete*culture + (1|pot), data= ramifications, family =nbinom1)
summary(mymodel)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Notations_structures_variete$variete)

Anova(mymodel, type=3)

emm<-emmeans(mymodel, ~ variete*culture, type = 'response')

emm<- as.data.table(emm)

emm

emm$exp <- paste("Glasshouse")

emm_GH <- emm %>%
  rename(predicted_count = response, SE_count =SE ,df_count = df ,  asymp.LCL_count =asymp.LCL, asymp.UCL_count=asymp.UCL )

#### Field experiment 2022-2023 ####

Notations_structures_variete <- read.csv("Data/Field_ecophysiological_traits.csv", sep=";")
Notations_structures_variete$variete <- ifelse(Notations_structures_variete$variete == "1", "angelico", ifelse(Notations_structures_variete$variete == "2","mambo","feliciano"))

Notations_structures_variete <-subset(Notations_structures_variete, nb_ramification!='NA')

Notations_structures_variete <-subset(Notations_structures_variete, date=='26.06.2023')

hist(Notations_structures_variete$nb_ramification)
Notations_structures_variete_asso <-subset(Notations_structures_variete, culture=='asso')
hist(Notations_structures_variete_asso$nb_ramification)
Notations_structures_variete_mono <-subset(Notations_structures_variete, culture=='mono')
hist(Notations_structures_variete_mono$nb_ramification)

Notations_structures_variete <-subset(Notations_structures_variete, variete!='4')
Notations_structures_variete$variete <- as.factor(Notations_structures_variete$variete)

mymodel<- glmmTMB(nb_ramification ~ variete*culture + (1|bloc), data= Notations_structures_variete, family =nbinom1, zi = ~variete*culture)
summary(mymodel)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Notations_structures_variete$variete)

Anova(mymodel, type=3)

Anova(mymodel, component = "zi", type = 3)

# Get estimated marginal means for the count model (conditional)
emm_count <- emmeans(mymodel, ~ variete * culture, type = "response", component = "cond")

# Get estimated marginal means for the zero-inflation model
emm_zi <- emmeans(mymodel, ~ variete * culture, type = "response", component = "zi")

# Convert to data frames
df_count <- as.data.frame(emm_count)
df_zi <- as.data.frame(emm_zi)

# Merge the two datasets on variete and culture
df_combined_Field_2023 <- left_join(df_count, df_zi, by = c("variete", "culture"), suffix = c("_count", "_zi"))

# Rename for clarity
df_combined_Field_2023 <- df_combined_Field_2023 %>%
  rename(predicted_count = response_count, zero_prob = response_zi)

df_combined_Field_2023$exp <- paste("Field 2022-2023")
df_combined_Field_2023

df_combined <- rbind(emm_GH, df_combined_Field_2023, fill=TRUE )
df_combined

##### Field 2023-2024 #####

Notations_structures_variete <- read.csv("Data/Field_ecophysiological_traits.csv", sep=";")
Notations_structures_variete$variete <- ifelse(Notations_structures_variete$variete == "1", "angelico", ifelse(Notations_structures_variete$variete == "2","mambo","feliciano"))
Notations_structures_variete <-subset(Notations_structures_variete, nb_ramification!='NA')

Notations_structures_variete <-subset(Notations_structures_variete, date=='02.04.2024')

hist(Notations_structures_variete$nb_ramification)

Notations_structures_variete$variete <- as.factor(Notations_structures_variete$variete)

mymodel<- glmmTMB(nb_ramification ~ variete*culture + (1|bloc), data= Notations_structures_variete, family =nbinom1, , zi = ~variete*culture)
summary(mymodel)

n_sim <- 500

simulationOutput <- simulateResiduals(fittedModel = mymodel, n = n_sim)
plot(simulationOutput, asFactor = F)

testDispersion(simulationOutput)

plot(simulationOutput, form = Notations_structures_variete$variete)

Anova(mymodel, type=3)
Anova(mymodel, component = "zi", type = 3)

# Get estimated marginal means for the count model (conditional)
emm_count <- emmeans(mymodel, ~ variete * culture, type = "response", component = "cond")

# Get estimated marginal means for the zero-inflation model
emm_zi <- emmeans(mymodel, ~ variete * culture, type = "response", component = "zi")

# Convert to data frames
df_count <- as.data.frame(emm_count)
df_zi <- as.data.frame(emm_zi)

# Merge the two datasets on variete and culture
df_combined_Field_2024 <- left_join(df_count, df_zi, by = c("variete", "culture"), suffix = c("_count", "_zi"))

# Rename for clarity
df_combined_Field_2024 <- df_combined_Field_2024 %>%
  rename(predicted_count = response_count, zero_prob = response_zi)

df_combined_Field_2024$exp <- paste("Field 2023-2024")
df_combined_Field_2024

emmean <- rbind(df_combined, df_combined_Field_2024, fill=TRUE)
emmean

#Figure for mambo 
emmean<-as.data.frame(emmean)
emmean
emmean_mambo <- subset(emmean, variete == "mambo")

plot_mambo <- ggplot(emmean_mambo, aes(x = exp, group = culture)) +
  scale_x_discrete(limits = c("Glasshouse", "Field 2022-2023", "Field 2023-2024")) +
  # Bar plot for predicted count with proper dodge width
  geom_bar(aes(y = predicted_count, fill = culture), 
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  # Error bars for count model (bars)
  geom_errorbar(aes(ymin = predicted_count - SE_count, ymax = predicted_count + SE_count, 
                    group = culture), 
                width = 0.2, position = position_dodge(width = 0.9), color = "black") +
  
  scale_fill_manual(values = c('#5773CCFF', '#FFB900FF')) +
  # Zero-inflation probability points with error bars, properly positioned
  geom_point(aes(y = zero_prob * 10, color = culture),
             position = position_dodge(width = 0.9), size = 3) +
  # Error bars for zero-inflation points (matching point colors)
  geom_errorbar(aes(ymin = (zero_prob - SE_zi) * 10, 
                    ymax = (zero_prob + SE_zi) * 10, 
                    group = culture, color = culture),
                width = 0.2, position = position_dodge(width = 0.9)) +
  # Set left y-axis (0 to 10, step of 2) and scale right y-axis (0 to 1)
  scale_y_continuous( breaks = seq(0, 10, by = 2),
                     sec.axis = sec_axis(~ . / 10, name = "Zero Inflation Probability (Dots)")) +
  scale_color_manual(values = c('#2F4A99', '#C28800')) +
  
  labs(y = "EMMeans of number of branching with count model (Bars)",
       fill = "Culture",
       color = "Culture") +
  theme_classic()

plot_mambo

write_rds(plot_mambo, "~/working/Aboveground_characterization/3-Output/Branching_exp_mambo.rds")
read_rds("~/working/Aboveground_characterization/3-Output/Branching_exp_mambo.rds")

#Figure for feliciano 
emmean<-as.data.frame(emmean)
emmean
emmean_feliciano <- subset(emmean, variete == "feliciano")

plot_feliciano <- ggplot(emmean_feliciano, aes(x = exp, group = culture)) +
  scale_x_discrete(limits = c("Glasshouse", "Field 2022-2023", "Field 2023-2024")) +
  # Bar plot for predicted count with proper dodge width
  geom_bar(aes(y = predicted_count, fill = culture), 
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  # Error bars for count model (bars)
  geom_errorbar(aes(ymin = predicted_count - SE_count, ymax = predicted_count + SE_count, 
                    group = culture), 
                width = 0.2, position = position_dodge(width = 0.9), color = "black") +
  
  scale_fill_manual(values = c('#5773CCFF', '#FFB900FF')) +
  # Zero-inflation probability points with error bars, properly positioned
  geom_point(aes(y = zero_prob * 10, color = culture),
             position = position_dodge(width = 0.9), size = 3) +
  # Error bars for zero-inflation points (matching point colors)
  geom_errorbar(aes(ymin = (zero_prob - SE_zi) * 10, 
                    ymax = (zero_prob + SE_zi) * 10, 
                    group = culture, color = culture),
                width = 0.2, position = position_dodge(width = 0.9)) +
  # Set left y-axis (0 to 10, step of 2) and scale right y-axis (0 to 1)
  scale_y_continuous(breaks = seq(0, 10, by = 2),
                      sec.axis = sec_axis(~ . / 10, name = "Zero Inflation Probability (Dots)")) +
  
  scale_color_manual(values = c('#2F4A99', '#C28800')) +
  
  labs(y = "EMMeans of number of branching with count model (Bars)",
       fill = "Culture",
       color = "Culture") +
  theme_classic()

plot_feliciano
write_rds(plot_feliciano, "~/working/Aboveground_characterization/3-Output/Branching_exp_feliciano.rds")

#Figure for angelico 
emmean<-as.data.frame(emmean)
emmean
emmean_angelico <- subset(emmean, variete == "angelico")

plot_angelico <- ggplot(emmean_angelico, aes(x = exp, group = culture)) +
  scale_x_discrete(limits = c("Glasshouse", "Field 2022-2023", "Field 2023-2024")) +
  # Bar plot for predicted count with proper dodge width
  geom_bar(aes(y = predicted_count, fill = culture), 
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  # Error bars for count model (bars)
  geom_errorbar(aes(ymin = predicted_count - SE_count, ymax = predicted_count + SE_count, 
                    group = culture), 
                width = 0.2, position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = c('#5773CCFF', '#FFB900FF')) +
  # Zero-inflation probability points with error bars, properly positioned
  geom_point(aes(y = zero_prob * 10, color = culture),
             position = position_dodge(width = 0.9), size = 3) +
  # Error bars for zero-inflation points (matching point colors)
  geom_errorbar(aes(ymin = (zero_prob - SE_zi) * 10, 
                    ymax = (zero_prob + SE_zi) * 10, 
                    group = culture, color = culture),
                width = 0.2, position = position_dodge(width = 0.9)) +
  # Set left y-axis (0 to 10, step of 2) and scale right y-axis (0 to 1)
  scale_y_continuous(breaks = seq(0, 12, by = 2),
                     sec.axis = sec_axis(~ . / 12, name = "Zero Inflation Probability (Dots)")) +
  scale_color_manual(values = c('#2F4A99', '#C28800')) +
  
  labs(y = "EMMeans of number of branching with count model (Bars)",
       fill = "Culture",
       color = "Culture") +
  theme_classic()

plot_angelico
write_rds(plot_angelico, "~/working/Aboveground_characterization/3-Output/Branching_exp_angelico.rds")
