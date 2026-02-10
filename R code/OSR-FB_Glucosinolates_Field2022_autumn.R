library (broom)
library(rstatix)
library(ggplot2)

EVA_Gluco <- read.csv("Data/Glucosinolates.csv", sep=";")
head(EVA_Gluco)

EVA_Gluco$date <- as.factor(EVA_Gluco$date)

EVA_Gluco <- subset(EVA_Gluco, EVA_Gluco$total > 60)

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
EVA_Gluco$total<-as.numeric(as.character(EVA_Gluco$total))
EVA_Gluco$Glucoraphanin<-as.numeric(as.character(EVA_Gluco$Glucoraphanin))


EVA_Gluco <- subset(EVA_Gluco, EVA_Gluco$Sampling_campain == "EVA_1_BF")
EVA_Gluco <- subset(EVA_Gluco, EVA_Gluco$variety != "Ex")

EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)
EVA_Gluco$bloc <- as.factor(EVA_Gluco$bloc)

head(EVA_Gluco)

# Analysis culture and variety on Total gluco ####
head(EVA_Gluco)

model <- lm(total ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

anova_table <- EVA_Gluco  %>% anova_test(total ~ bloc + variety*culture )
anova_table
anova_results <- data.frame(variety_p = numeric(), culture_p = numeric(), CxV_p = numeric())

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~variety*culture)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety*culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

plot_total <- ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucosinolates Total  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

total <- plot_total +  theme( legend.position="none")

total

# Analysis culture and variety on Total Glucoraphanin ####

head(EVA_Gluco)

EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)

model <- lm(log(Glucoraphanin+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Glucoraphanin <- log(EVA_Gluco$Glucoraphanin+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Glucoraphanin ~ bloc + variety*culture )
anova_table

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~culture)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-emmeans(model, ~ culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-emmeans(model, ~ culture*variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange", name = "", labels = c("Intercropping", "Monocropping"))) + 
  labs(x = "EMMeans of Glucoraphanin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

plot_Glucoraphanin <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucoraphanin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Glucoraphanin <- plot_Glucoraphanin +  theme( legend.position="none")

Glucoraphanin

# Progoitrin.isomer####

head(EVA_Gluco)

EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)

model <- lm(log(Progoitrin.isomer+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Progoitrin.isomer <- log(EVA_Gluco$Progoitrin.isomer+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Progoitrin.isomer ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

emm<-emmeans(model, ~ culture*variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

plot_Progoitrin.isomer<- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Progoitrin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Progoitrin <- plot_Progoitrin.isomer +  theme( legend.position="none")

Progoitrin

# Glucoalyssin ####

head(EVA_Gluco)

EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)

model <- lm(Glucoalyssin ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

anova_table <- EVA_Gluco  %>% anova_test(Glucoalyssin ~ bloc + variety*culture )
anova_table

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

emm<-emmeans(model, ~ culture*variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

plot_Glucoalyssin <- ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucoalyssin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Glucoalyssin <- plot_Glucoalyssin +  theme( legend.position="none")

Glucoalyssin

# Gluconapoleiferin ####

head(EVA_Gluco)

EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)

model <- lm(log(Gluconapoleiferin+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Gluconapoleiferin <- log(EVA_Gluco$Gluconapoleiferin+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Gluconapoleiferin ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~culture)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-emmeans(model, ~ culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-emmeans(model, ~ culture*variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange", name = "", labels = c("Intercropping", "Monocropping"))) + 
  labs(x = "EMMeans of Glucoraphanin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

plot_Gluconapoleiferin <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucoraphanin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Gluconapoleiferin <- plot_Gluconapoleiferin +  theme( legend.position="none")

Gluconapoleiferin

#Gluconapin####

head(EVA_Gluco)

EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)

model <- lm(Gluconapin ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

anova_table <- EVA_Gluco  %>% anova_test(Gluconapin ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

emm<-emmeans(model, ~ culture*variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

plot_Gluconapin<- ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Gluconapinn  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Gluconapin <- plot_Gluconapin +  theme( legend.position="none")

Gluconapin

# Analysis culture and variety on Total Butyl####
head(EVA_Gluco)

EVA_Gluco$variety <- as.factor(EVA_Gluco$variety)
EVA_Gluco_butyl <- subset(EVA_Gluco, EVA_Gluco$Sample.ID != "230112_IH_GLS_sample_4-44")
EVA_Gluco_butyl <- subset(EVA_Gluco_butyl, EVA_Gluco_butyl$Sample.ID != "230112_IH_GLS_sample_4-71")

model <- lm(log(Butyl.GS+1) ~ bloc + variety*culture  , data=EVA_Gluco_butyl ) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

anova_table <- EVA_Gluco  %>% anova_test(Butyl.GS ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~variety)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety*culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange", name = "", labels = c("Intercropping", "Monocropping"))) + 
  labs(x = "Butyl glucosinolate  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()


plot_Butyl <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Butyl glucosinolate  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Butyl <- plot_Butyl+  theme( legend.position="none")

Butyl

# Analysis culture and variety on Total Glucobrassicanapin ####

model <- lm(log(Glucobrassicanapin+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Glucobrassicanapin <- log(EVA_Gluco$Glucobrassicanapin+1)

anova_table <- EVA_Gluco  %>% anova_test(Glucobrassicanapin ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~variety)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety*culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange" )) + 
  labs(x = "EMMeans of Glucobrassicanapin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

plot_Glucobrassicanapin <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucobrassicanapin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Glucobrassicanapin <- plot_Glucobrassicanapin +  theme( legend.position="none")

Glucobrassicanapin

# Analysis culture and variety on Total Glucobrassicin ####
head(EVA_Gluco)

model <- lm(log(Glucobrassicin+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Glucobrassicin <- log(EVA_Gluco$Glucobrassicin+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Glucobrassicin ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~variety)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety*culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange")) + 
  labs(x = "EMMeans of Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()


plot_Glucobrassicin <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Glucobrassicin <- plot_Glucobrassicin +  theme( legend.position="none")

Glucobrassicin

#Gluconasturtiin####

head(EVA_Gluco)

model <- lm(log(Gluconasturtiin+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Gluconasturtiin <- log(EVA_Gluco$Gluconasturtiin+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Gluconasturtiin ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

emm<-emmeans(model, ~ culture*variety, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

plot_Gluconasturtiin <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Gluconasturtiin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Gluconasturtiin <- plot_Gluconasturtiin +  theme( legend.position="none")
Gluconasturtiin 

# Methoxyglucobrassicin ####
head(EVA_Gluco)

model <- lm(log(Methoxyglucobrassicin+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Methoxyglucobrassicin <- log(EVA_Gluco$Methoxyglucobrassicin+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Methoxyglucobrassicin ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~variety)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety*culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange")) + 
  labs(x = "EMMeans of Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()


plot_Methoxyglucobrassicin <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Methoxyglucobrassicin <- plot_Methoxyglucobrassicin +  theme( legend.position="none")

Methoxyglucobrassicin

# Neoglucobrassicin ####

head(EVA_Gluco)

model <- lm(sqrt(Neoglucobrassicin+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Neoglucobrassicin <- sqrt(EVA_Gluco$Neoglucobrassicin+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Neoglucobrassicin ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~variety)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety*culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange")) + 
  labs(x = "EMMeans of Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()


plot_Neoglucobrassicin <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

Neoglucobrassicin <- plot_Neoglucobrassicin +  theme( legend.position="none")

Neoglucobrassicin

#Unknown.GLS ####

head(EVA_Gluco)

model <- lm(log(Unknown.GLS..C16H20N2O11S2.+1) ~ bloc + variety*culture  , data=EVA_Gluco) #calculer le model, la covariable passe en premier 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid)
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ variety*culture) #OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #OK

EVA_Gluco$log_Unknown.GLS <- log(EVA_Gluco$Unknown.GLS..C16H20N2O11S2.+1)

anova_table <- EVA_Gluco  %>% anova_test(log_Unknown.GLS ~ bloc + variety*culture )

# Function to extract p-values and append to results table
variety_p <- anova_table$p[anova_table$Effect == "variety"]
culture_p <- anova_table$p[anova_table$Effect == "culture"]
CxV_p <- anova_table$p[anova_table$Effect == "variety:culture"]

# Append to the global results table
anova_results <- rbind(anova_results, data.frame(variety_p, culture_p, CxV_p ))

library(emmeans)

pwc <- 
  emmeans(model,pairwise~culture)
pwc
pairwise <- pwc$contrasts

emm<-emmeans(model, ~ variety*culture, type = 'response')

plot(emm, by = NULL, comparisons = FALSE, adjust = "tukey", horizontal = FALSE)

emm<-as.data.frame(emm)
emm

ggplot(emm, aes(x =emmean , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 2, alpha=0.5) + 
  scale_color_manual(values = c("blue","orange")) + 
  labs(x = "EMMeans of Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()


plot_log_Unknown.GLS <- ggplot(emm, aes(x =response , y = variety, color = culture)) + 
  geom_point(position = position_dodge(width = 0.5), size=3, alpha = 1) + 
  geom_linerange(aes(xmin = lower.CL, xmax = upper.CL),position = position_dodge(width = 0.5),size = 4, alpha=0.3) + 
  scale_color_manual(values = c("blue","orange"), name = "", labels = c("Intercropping", "Monocropping")) + 
  labs(x = "Glucobrassicin  (ug/g)", y = "variety") + 
  theme_classic() + coord_flip()

log_Unknown.GLS <- plot_log_Unknown.GLS +  theme( legend.position="none")

log_Unknown.GLS

# FIGURE ####

write_csv(anova_results, "~/working/Aboveground_characterization/3-Output/anova_results_EVA_2022.csv" )

library(ggpubr)

figure <- ggarrange( total, Butyl, Glucoraphanin, Glucobrassicin, Glucobrassicanapin, 
                    Glucoalyssin, Gluconasturtiin, Gluconapin, Progoitrin,
                    labels = c("A", "B", "C", "D", "E", "F", "G","H", "I", "J"), 
                    ncol = 2, nrow = 5, 
                    common.legend = TRUE, legend = "bottom", 
                    widths = 2)
figure
