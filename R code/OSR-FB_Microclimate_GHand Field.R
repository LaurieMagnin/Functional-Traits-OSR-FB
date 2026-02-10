library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(rstatix)
library('data.table')
library(emmeans)

### Analysis in temperature difference with air ############
#### GlassHouse####
# Temperature campbell est side 
temp_BF_est<- read.csv("Data/GlassHouse_Temperature_est.csv", sep=",")

temp_BF_est <- subset(temp_BF_est, temp_BF_est$modality != "OUT")
temp_BF_est <- subset(temp_BF_est, temp_BF_est$date != "2022-04-05")

temp_BF_est <- subset(temp_BF_est, temp_BF_est$date == "2022-03-31" | 
                        temp_BF_est$date == "2022-04-01" | 
                        temp_BF_est$date == "2022-04-02" | temp_BF_est$date == "2022-04-03" | temp_BF_est$date == "2022-04-04")

temp_BF_est  <- temp_BF_est  %>%                                 # Group data
  group_by(culture, date, sonde, modality) %>%
  dplyr::summarize(cumulated_D = mean(value), n=n()) %>% 
  as.data.frame()

temp_BF_est <- subset(temp_BF_est,temp_BF_est$cumulated_D>0)

#Data transformation in relative cumulated temperature 

temp_BF_est_ext  <- temp_BF_est  %>%                                 # Group data
  group_by(culture, date, modality) %>%
  dplyr::summarize(cumulated_D = mean(cumulated_D)) %>% 
  as.data.frame()

temp_BF_est_ext  <- subset(temp_BF_est_ext,temp_BF_est$modality == "Ext")
temp_BF_est_ext  <- subset(temp_BF_est_ext, select = c(2:4))
temp_BF_est_ext  <- temp_BF_est_ext %>% slice(1:5)
temp_BF_est_ext  <- spread(temp_BF_est_ext, key= date, value= cumulated_D)

temp_BF_est_spread  <- spread(temp_BF_est, key= date, value= cumulated_D)
temp_BF_est_spread$day1 <- temp_BF_est_spread$"2022-03-31"-temp_BF_est_ext$"2022-03-31" 
temp_BF_est_spread$day2 <- temp_BF_est_spread$"2022-04-01"- temp_BF_est_ext$"2022-04-01" 
temp_BF_est_spread$day3 <- temp_BF_est_spread$"2022-04-02"- temp_BF_est_ext$"2022-04-02" 
temp_BF_est_spread$day4 <- temp_BF_est_spread$"2022-04-03"-temp_BF_est_ext$"2022-04-03" 
temp_BF_est_spread$day5 <- temp_BF_est_spread$"2022-04-04"- temp_BF_est_ext$"2022-04-04" 

# temp_BF_est_corr <- gather(temp_BF_est_spread, "date", "value", 5:9) 
temp_BF_est_corr <- gather(temp_BF_est_spread, "day", "T_cum_corr", 10:14) 
temp_BF_est_corr  <- subset(temp_BF_est_corr, select = -c(5:9))
temp_BF_est_corr <- subset(temp_BF_est_corr, temp_BF_est_corr$T_cum_corr != "NA")
temp_BF_est_corr <- subset(temp_BF_est_corr, temp_BF_est_corr$modality != "Ext")

# Temperature campbell west side 

temp_BF_west<- read.csv("Data/GlassHouse_Temperature_west.csv", sep=",")

temp_BF_west <- subset(temp_BF_west, temp_BF_west$modality != "Other")
temp_BF_west <- subset(temp_BF_west, temp_BF_west$date == "2022-03-31" | 
                         temp_BF_west$date == "2022-04-01" | 
                         temp_BF_west$date == "2022-04-02" | temp_BF_west$date == "2022-04-03" | temp_BF_west$date == "2022-04-04")

#Exploration

temp_BF_west  <- temp_BF_west  %>%                                 # Group data
  group_by(culture, date, sonde, modality) %>%
  dplyr::summarize(cumulated_D = mean(value), n=n()) %>% 
  as.data.frame()


temp_BF_west <- subset(temp_BF_west,temp_BF_west$cumulated_D>0)
temp_BF_west <- subset(temp_BF_west,temp_BF_west$sonde != "T_13")

#Data transformation in relative cumulated temperature 

temp_BF_west_ext  <- temp_BF_west  %>%                                 # Group data
  group_by(culture, date, modality) %>%
  dplyr::summarize(cumulated_D = mean(cumulated_D)) %>% 
  as.data.frame()

temp_BF_west_ext  <- subset(temp_BF_west_ext,temp_BF_west$modality == "Ext")
temp_BF_west_ext  <- subset(temp_BF_west_ext, select = c(2:4))
temp_BF_west_ext  <- temp_BF_west_ext %>% slice(1:5)
temp_BF_west_ext  <- spread(temp_BF_west_ext, key= date, value= cumulated_D)

temp_BF_west_spread  <- spread(temp_BF_west, key= date, value= cumulated_D)
temp_BF_west_spread$day1 <- temp_BF_west_spread$"2022-03-31"- temp_BF_west_ext$"2022-03-31" 
temp_BF_west_spread$day2 <- temp_BF_west_spread$"2022-04-01"-  temp_BF_west_ext$"2022-04-01" 
temp_BF_west_spread$day3 <- temp_BF_west_spread$"2022-04-02"- temp_BF_west_ext$"2022-04-02" 
temp_BF_west_spread$day4 <- temp_BF_west_spread$"2022-04-03"- temp_BF_west_ext$"2022-04-03" 
temp_BF_west_spread$day5 <- temp_BF_west_spread$"2022-04-04"- temp_BF_west_ext$"2022-04-04" 

# temp_BF_west_corr <- gather(temp_BF_west_spread, "date", "value", 5:9) 
temp_BF_west_corr <- gather(temp_BF_west_spread, "day", "T_cum_corr", 10:14) 
temp_BF_west_corr  <- subset(temp_BF_west_corr, select = -c(5:9))
temp_BF_west_corr <- subset(temp_BF_west_corr, temp_BF_west_corr$T_cum_corr != "NA")
temp_BF_west_corr <- subset(temp_BF_west_corr, temp_BF_west_corr$modality != "Ext")

# Merging the datasets of the sides of the glass house 

T_BF_campbell <- rbind(temp_BF_west_corr,temp_BF_est_corr )
head(T_BF_campbell)

# Stat analysis
model <- lm(T_cum_corr ~ day + culture , data=T_BF_campbell) 
#Normalité des résidus 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid) #  OK 
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ culture) # OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() # OK

T_BF_campbell$log_T_cum_corr <- log(T_BF_campbell$T_cum_corr +1)

T_BF_campbell  %>% anova_test(log_T_cum_corr ~ culture)

emm<-emmeans(model, ~ culture,, type = 'response')

emm<- as.data.table(emm)
emm

Temperature_GH <- ggplot(emm, aes(y=emmean, x= culture, color=culture )) + 
  geom_point()+
  geom_errorbar(aes(ymax = emmean+SE, ymin = emmean-SE, width = .02))+ 
  theme_classic()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  scale_x_discrete(labels = c("As" = "Intercropping", "Mo" = "Monocropping")) +
  labs(title="Glasshouse ", y = "EMMeans temperature difference (°C)")+  theme( legend.position="none")
Temperature_GH

library(readr)
saveRDS(Temperature_GH, "~/working/Aboveground_characterization/3-Output/DeltaT_GH.rds")

# Field 2022-2023 ####

#Data transformation in relative cumulated temperature 

T_2022 <- read.csv("Data/Field_2022_Temperature.csv", sep=",")
head(T_2022)

mclimate_BF  <- T_2022  %>%                                 # Group data
  group_by(day, sonde, culture) %>%
  dplyr::summarize(cumulated_D = mean(value)) %>% 
  as.data.frame()

mclimate_BF<- subset(mclimate_BF, mclimate_BF$cumulated_D != "NA")

mclimate_BF <- subset(mclimate_BF, mclimate_BF$day == "08.11.2022" | 
                        mclimate_BF$day == "09.11.2022" | mclimate_BF$day == "10.11.2022" |
                        mclimate_BF$day == "11.11.2022" | mclimate_BF$day == "12.11.2022" | 
                        mclimate_BF$day == "13.11.2022"| mclimate_BF$day == "14.11.2022" | mclimate_BF$day == "15.11.2022")


mclimate_BF_ext  <- mclimate_BF  %>%                                 # Group data
  group_by(culture, day) %>%
  dplyr::summarize(cumulated_D = mean(cumulated_D)) %>% 
  as.data.frame()

mclimate_BF_ext   <- subset(mclimate_BF_ext ,mclimate_BF_ext $culture == "ctrl")
mclimate_BF_ext  <- subset(mclimate_BF_ext, select = c(2:4))
mclimate_BF_ext  <- mclimate_BF_ext %>% slice(1:8)
mclimate_BF_ext  <- spread(mclimate_BF_ext, key= day, value= cumulated_D)

mclimate_BF_spread  <- spread(mclimate_BF, key= day, value= cumulated_D)
mclimate_BF_spread  <- subset(mclimate_BF_spread, mclimate_BF_spread$culture != "OUT")
mclimate_BF_spread$day1 <- mclimate_BF_spread$"08.11.2022"- mclimate_BF_ext$"08.11.2022" 
mclimate_BF_spread$day2 <- mclimate_BF_spread$"09.11.2022"-mclimate_BF_ext$"09.11.2022" 
mclimate_BF_spread$day3 <- mclimate_BF_spread$"10.11.2022"-mclimate_BF_ext$"10.11.2022" 
mclimate_BF_spread$day4 <- mclimate_BF_spread$"11.11.2022"-mclimate_BF_ext$"11.11.2022" 
mclimate_BF_spread$day5 <- mclimate_BF_spread$"12.11.2022"-mclimate_BF_ext$"12.11.2022" 
mclimate_BF_spread$day6 <- mclimate_BF_spread$"13.11.2022"-mclimate_BF_ext$"13.11.2022" 
mclimate_BF_spread$day7 <- mclimate_BF_spread$"14.11.2022"-mclimate_BF_ext$"14.11.2022" 
mclimate_BF_spread$day8 <- mclimate_BF_spread$"15.11.2022"- mclimate_BF_ext$"15.11.2022" 


mclimate_BF_spread  <- subset(mclimate_BF_spread, select = -c(3:10))
mclimate_BF_corr <- gather(mclimate_BF_spread, "day", "T_cum_corr", 3:10) 

mclimate_BF_corr <-subset(mclimate_BF_corr, mclimate_BF_corr$culture != "ctrl")
mclimate_BF_corr <-subset(mclimate_BF_corr, mclimate_BF_corr$sonde != "T_13")
head(mclimate_BF_corr)

# Stat analysis
mclimate_BF_corr <- subset(mclimate_BF_corr, mclimate_BF_corr$T_cum_corr >-0.3414801)

model <- lm(log(T_cum_corr+1) ~  culture , data=mclimate_BF_corr)
#Normalité des résidus 
model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid) #  OK 
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ culture) # OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #not OK

mclimate_BF_corr  %>% anova_test(T_cum_corr ~ culture)

emm<-emmeans(model, ~ culture,, type = 'response')

emm<- as.data.table(emm)
emm

Temperature_2023 <- ggplot(emm, aes(y = response, x = culture, color = culture)) + 
  geom_point() +
  geom_errorbar(aes(ymax = response + SE, ymin = response - SE, width = 0.02)) + 
  theme_classic() +
  scale_x_discrete(labels = c("asso" = "Intercropping", "mono" = "Monocropping")) +
  scale_color_manual(values = c('#5773CCFF', '#FFB900FF')) +
  labs(title = "Field 2022-2023", y = "EMMeans temperature difference (°C)") +  
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-0.05, 0.25), breaks = seq(0, 0.3, by = 0.1))

Temperature_2023

saveRDS(Temperature_2023, "~/working/Aboveground_characterization/3-Output/DeltaT_2022.rds")

# Field 2023-2024
#Data transformation in relative cumulated temperature 

T_2023 <- read.csv("Data/Field_2023_Temperature.csv", sep=",")
head(T_2023)

mclimate_BF  <- T_2023  %>%                                 # Group data
  group_by(day, sonde, culture) %>%
  dplyr::summarize(cumulated_D = mean(value)) %>% 
  as.data.frame()

mclimate_BF<- subset(mclimate_BF, mclimate_BF$cumulated_D > 0)

mclimate_BF <- subset(mclimate_BF, mclimate_BF$day == "08.11.2023 00:00" | 
                        mclimate_BF$day == "09.11.2023 00:00" | mclimate_BF$day == "10.11.2023 00:00" |
                        mclimate_BF$day == "11.11.2023 00:00" | mclimate_BF$day == "12.11.2023 00:00" | 
                        mclimate_BF$day == "13.11.2023 00:00"| mclimate_BF$day == "14.11.2023 00:00" | mclimate_BF$day == "15.11.2023 00:00")


mclimate_BF_ext  <- mclimate_BF  %>%                                 # Group data
  group_by(culture, day) %>%
  dplyr::summarize(cumulated_D = mean(cumulated_D)) %>% 
  as.data.frame()

mclimate_BF_ext   <- subset(mclimate_BF_ext ,mclimate_BF_ext$culture == "ctrl")
mclimate_BF_ext  <- subset(mclimate_BF_ext, select = c(2:4))
mclimate_BF_ext  <- mclimate_BF_ext %>% slice(1:8)
mclimate_BF_ext  <- spread(mclimate_BF_ext, key= day, value= cumulated_D)

mclimate_BF_spread  <- spread(mclimate_BF, key= day, value= cumulated_D)
mclimate_BF_spread  <- subset(mclimate_BF_spread, mclimate_BF_spread$culture != "OUT")
mclimate_BF_spread$day1 <- mclimate_BF_spread$"08.11.2023 00:00"-mclimate_BF_ext$"08.11.2023 00:00" 
mclimate_BF_spread$day2 <- mclimate_BF_spread$"09.11.2023 00:00"- mclimate_BF_ext$"09.11.2023 00:00" 
mclimate_BF_spread$day3 <- mclimate_BF_spread$"10.11.2023 00:00"- mclimate_BF_ext$"10.11.2023 00:00" 
mclimate_BF_spread$day4 <- mclimate_BF_spread$"11.11.2023 00:00"- mclimate_BF_ext$"11.11.2023 00:00" 
mclimate_BF_spread$day5 <- mclimate_BF_spread$"12.11.2023 00:00"- mclimate_BF_ext$"12.11.2023 00:00" 
mclimate_BF_spread$day6 <- mclimate_BF_spread$"13.11.2023 00:00"- mclimate_BF_ext$"13.11.2023 00:00" 
mclimate_BF_spread$day7 <- mclimate_BF_spread$"14.11.2023 00:00"- mclimate_BF_ext$"14.11.2023 00:00" 
mclimate_BF_spread$day8 <- mclimate_BF_spread$"15.11.2023 00:00"- mclimate_BF_ext$"15.11.2023 00:00" 


mclimate_BF_spread  <- subset(mclimate_BF_spread, select = -c(3:10))
mclimate_BF_corr <- gather(mclimate_BF_spread, "day", "T_cum_corr", 3:10) 

mclimate_BF_corr <-subset(mclimate_BF_corr, mclimate_BF_corr$culture != "ctrl")
mclimate_BF_corr <-subset(mclimate_BF_corr, mclimate_BF_corr$T_cum_corr != "NA")
mclimate_BF_corr <-subset(mclimate_BF_corr, mclimate_BF_corr$sonde != "T_13")

# Stat analysis
model <- lm(T_cum_corr ~ culture , data=mclimate_BF_corr) 
#Normalité des résidus 

model.metrics <- augment(model) 
head(model.metrics)

shapiro_test(model.metrics$.resid) #  OK 
qqnorm(model.metrics$.resid)
qqline(model.metrics$.resid)

#Homogénéité des variances 

model.metrics %>% levene_test(.resid ~ culture) # OK

# Valeurs aberrantes 

model.metrics %>% filter(abs(.std.resid) > 3) %>% as.data.frame() #not OK

mclimate_BF_corr  %>% anova_test(T_cum_corr ~  culture)

emm<-emmeans(model, ~ culture,, type = 'response')

emm<- as.data.table(emm)
emm

Temperature_2024 <- ggplot(emm, aes(y=emmean, x= culture, color=culture )) + 
  geom_point()+
  geom_errorbar(aes(ymax = emmean+SE, ymin = emmean-SE, width = .02))+ 
  theme_classic()+
  scale_color_manual(values=c('#5773CCFF','#FFB900FF'))+
  labs(title="Field 2023-2024 ", y = "EMMeans temperature difference (°C)")+
  scale_y_continuous(limits = c(-0.30, 0))+
  scale_x_discrete(labels = c("asso" = "Intercropping", "mono" = "Monocropping"))+  theme( legend.position="none")

Temperature_2024

saveRDS(Temperature_2024, "~/working/Aboveground_characterization/3-Output/DeltaT_2023.rds")
