library (broom)
library(rstatix)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggstatsplot)
library(viridis) # or RColorBrewer

EVA_Gluco <- read.csv("Data/Glucosinolates.csv", sep=";")
head(EVA_Gluco)
EVA_Gluco <- subset(EVA_Gluco, EVA_Gluco$total > 60)
EVA_Gluco <- subset(EVA_Gluco, EVA_Gluco$variety != 'Ex')

EVA_Gluco$Progoitrin.isomer<-as.numeric(as.character(EVA_Gluco$Progoitrin.isomer))
EVA_Gluco$Glucoalyssin<-as.numeric(as.character(EVA_Gluco$ Glucoalyssin))
EVA_Gluco$Gluconapoleiferin<-as.numeric(as.character(EVA_Gluco$Gluconapoleiferin))
EVA_Gluco$Gluconapin<-as.numeric(as.character(EVA_Gluco$Gluconapin))
EVA_Gluco$Butyl.GS<-as.numeric(as.character(EVA_Gluco$Butyl.GS))
EVA_Gluco$Glucobrassicanapin<-as.numeric(as.character(EVA_Gluco$Glucobrassicanapin))
EVA_Gluco$Glucobrassicin<-as.numeric(as.character(EVA_Gluco$Glucobrassicin))
EVA_Gluco$Gluconasturtiin<-as.numeric(as.character(EVA_Gluco$Gluconasturtiin))
EVA_Gluco$Methoxyglucobrassicin<-as.numeric(as.character(EVA_Gluco$Methoxyglucobrassicin))
EVA_Gluco$Neoglucobrassicin<-as.numeric(as.character(EVA_Gluco$Neoglucobrassicin))
EVA_Gluco$Unknown.GLS..C16H20N2O11S2.<-as.numeric(as.character(EVA_Gluco$Unknown.GLS..C16H20N2O11S2.))
EVA_Gluco$Total<-as.numeric(as.character(EVA_Gluco$Total))
EVA_Gluco$Glucoraphanin<-as.numeric(as.character(EVA_Gluco$Glucoraphanin))

EVA_Gluco$date <- as.factor(EVA_Gluco$date)

##### Analysis of total glucosinolate content across experiments (Fig.5A)
head(EVA_Gluco)
EVA_Gluco <- subset(EVA_Gluco, EVA_Gluco$sample_type == "feuille")

EVA_Gluco$Sampling_campain <- as.factor(EVA_Gluco$Sampling_campain)
EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)

model <- lm(log(total) ~ Sampling_campain + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

Anova(model, type=3)

library(emmeans)

pwc <- 
  emmeans(model,pairwise~Sampling_campain)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ Sampling_campain, type = 'response')

emm<-as.data.frame(emm)
emm

emm$Sampling_campain <- ifelse(emm$Sampling_campain=="GH_BF", "Glasshouse before wintering", 
                               ifelse(emm$Sampling_campain=="GH_AF", "Glasshouse after wintering",
                                      ifelse(emm$Sampling_campain=="EVA_1_BF", "Field 2022-2023 autumn",
                                             ifelse(emm$Sampling_campain=="EVA_2_BF", "Field 2023-2024 autumn", "Field 2023-2024 spring" ))))

emm$Sampling_campain <- factor(emm$Sampling_campain, levels = c("Glasshouse before wintering", "Glasshouse after wintering", "Field 2022-2023 autumn", "Field 2023-2024 autumn", "Field 2023-2024 spring"))

plot_total <- ggplot(emm, aes(x =response , y = Sampling_campain, color = Sampling_campain)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),linewidth = 4, alpha=0.3) + 
  scale_color_manual(values = c("#339900","#66CC00", "#660006", "#000099", "#0000FF" ), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucosinolates Total  (ug/g)", y = "Variety") + 
  theme_classic() + coord_flip() +
  theme(axis.title.x = element_blank(),  # Remove X-axis title
        axis.text.x = element_text(angle = 30, hjust = 1))  # Rotate X labels

total <- plot_total +  theme( legend.position="none")

total

write_rds(total, "~/working/Aboveground_characterization/3-Output/total.rds")

### Proportion of glucosinolates type across sampling campaigns (Fig. 5B) ####

EVA_Gluco_long <- EVA_Gluco %>% 
  pivot_longer(
    cols = 'Glucoraphanin':'total', 
    names_to = "glucosinolates",
    values_to = "concentration"
  )

EVA_Gluco_long <- EVA_Gluco_long %>%
  group_by(Sample.ID) %>%
  mutate(proportion = concentration / sum(concentration)) %>%
  ungroup()

EVA_Gluco_long <- subset(EVA_Gluco_long,EVA_Gluco_long$glucosinolates != "total")

EVA_Gluco_long <- subset(EVA_Gluco_long, EVA_Gluco_long$Sampling_campain !="GH_seed")

EVA_Gluco_long$Sampling_campain <- ifelse(EVA_Gluco_long$Sampling_campain=="GH_BF", "Glasshouse before wintering", 
                               ifelse(EVA_Gluco_long$Sampling_campain=="GH_AF", "Glasshouse after wintering",
                                      ifelse(EVA_Gluco_long$Sampling_campain=="EVA_1_BF", "Field 2022-2023 autumn",
                                             ifelse(EVA_Gluco_long$Sampling_campain=="EVA_2_BF", "Field 2023-2024 autumn", "Field 2023-2024 spring" ))))

EVA_Gluco_long$Sampling_campain <- factor(EVA_Gluco_long$Sampling_campain, levels = c("Glasshouse before wintering", "Glasshouse after wintering", "Field 2022-2023 autumn", "Field 2023-2024 autumn", "Field 2023-2024 spring"))

EVA_Gluco_long<- EVA_Gluco_long %>%
  group_by(Sampling_campain, glucosinolates, ) %>%
  dplyr::summarize(proportion = sum(proportion))%>% 
  as.data.frame()

Proportion_gluco <- ggplot(EVA_Gluco_long, aes(fill = glucosinolates, y = proportion, x = Sampling_campain)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(axis.title.x = element_blank(),  # Remove X-axis title
        axis.text.x = element_text(angle = 30, hjust = 1))  # Rotate X labels

Proportion_gluco

write_rds(Proportion_gluco, "~/working/Aboveground_characterization/3-Output/Proportion_gluco.rds")

## PLS-DA performed on the standardized concentration of each of the GLS type (Fig. 5C)

data_log <- EVA_Gluco
data_log[, 9:20] <- log(data_log[, 9:20]+1)
data_log <- subset(data_log, select = -c(total))
data_log <- subset(data_log, data_log$sample_type =="feuille" )

for (i in 9:20) {
  data_log[, i] <- scales::rescale(as.numeric(data_log[, i]), to = c(0, 1))
}

library(mixOmics)

X <- data_log [,9:20]
X<- data.matrix(X, rownames.force = NA)
data_log$Sampling_campain <- ifelse(data_log$Sampling_campain=="GH_BF", "Glasshouse before wintering", 
                                    ifelse(data_log$Sampling_campain=="GH_AF", "Glasshouse after wintering",
                                           ifelse(data_log$Sampling_campain=="EVA_1_BF", "Field 2022-2023 autumn",
                                                  ifelse(data_log$Sampling_campain=="EVA_2_BF", "Field 2023-2024 autumn", "Field 2023-2024 spring" ))))
Y <- data_log$Sampling_campain

plsda.gluco <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.gluco, ind.names = FALSE, ellipse = TRUE, legend = TRUE)
plot(plsda.gluco)

linn.vip <- vip(plsda.gluco)

biplot(plsda.gluco, ind.names= FALSE)

plotIndiv(plsda.gluco, ind.names = FALSE, ellipse = TRUE, ellipse.level = 0.95)

PLS_DA_gluco <- biplot(plsda.gluco, ind.names = FALSE,col.per.group = c("#660006", "#000099", "#0000FF" ,"#66CC00", "#339900"))
PLS_DA_gluco

write_rds(PLS_DA_gluco, "~/working/Aboveground_characterization/3-Output/PLS_DA_gluco.rds")
