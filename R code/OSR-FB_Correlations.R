library(dplyr)
library ('corrplot')
library(RColorBrewer)
library(Hmisc)
library("ggpubr")

#Correlation matrix (Fig. S2)

data <- read.csv2("Data/Data_Field_Summary.csv")
head(data)

data <- data[, !(colnames(data) %in% c("M", "year", "variety","culture", "bloc", "D_FB_february", "shotholes_NB_median"))]

sapply(data, is.numeric)

data$Collar_D_feb <- as.numeric(data$Collar_D_feb)

col_order <- c("Rdt.net", "Leaves_NB", "Collar_D", "biomasse", "Stem_elongation","Collar_D_feb","D_OSR_february", 
               "shotholes_NB_mean","larva.plant", "oviposition_NB", "Pollenbeetle_NB", 
               "Stem_height", "Ramification_NB", "PMG", "seeds.surface")
data <- data[, col_order]

data <- data %>% 
  rename(
    Yield = Rdt.net,
    Biomasse = biomasse, 
    Density_OSR_feb = D_OSR_february,
    Shotholes_cotyledons_NB = shotholes_NB_mean, 
    Larval.infestation = larva.plant, 
    Thousand_kernel_weight = PMG, 
    Seeds_per_smeter = seeds.surface
    
  )

correlation_coef <-cor(data, method= c("pearson"))

correlation <- rcorr(as.matrix(data))
correlation

colpractice<- colorRampPalette(c("orange", "white", "blue"))(10)

corrplot(correlation_coef, type = "upper", order = "hclust", tl.col = "black", tl.srt =45)

corrplot(correlation$r, type = "upper", p.mat = correlation$P, sig.level= 0.05,tl.col = "black", insig= "blank", tl.srt =45, addCoef.col = 'black', col=colpractice)

correlation_P <- correlation$P
correlation_R <- correlation$r

# Zoom IN, correlation plot (Fig. 4)####
data <- read.csv2("~/working/Aboveground_characterization/1-Data/Data_Field_Summary.csv")
head(data)

data$year <- as.factor(data$year)

Corr_plot_ramif <- ggscatter(data, x = "Ramification_NB", y = "Stem_height",
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          color = "year",            # Color by groups 
          shape = "culture", 
          xlab = "Average number of branching", ylab = "Average stem height (cm)"
)+
  scale_color_manual(values = c('black','#999999'))+
  scale_fill_manual(values = c('black','#999999'))+
  stat_cor(aes(color = year), label.x = 2) 

Corr_plot_ramif

Corr_plot_rdt <-ggscatter(data, x = "Rdt.net", y = "Stem_height",
          add = "reg.line",                         
          conf.int = TRUE,                         
          color = "year",            
          shape = "culture", 
          xlab = " Yield (dt/ha)", ylab = "Average stem height (cm)")+
  scale_color_manual(values = c('black','#999999'))+
  scale_fill_manual(values = c('black','#999999'))+
  stat_cor(aes(color = year), label.x = 2) 

Corr_plot_rdt

Corr_plot_ramif <- ggscatter(data, x = "Ramification_NB", y = "Rdt.net",
                             add = "reg.line",                         # Add regression line
                             conf.int = TRUE,                          # Add confidence interval
                             color = "year",            # Color by groups 
                             shape = "culture", 
                             xlab = "Average number of branching", ylab = "Average stem height (cm)"
)+
  scale_color_manual(values = c('black','#999999'))+
  scale_fill_manual(values = c('black','#999999'))+
  stat_cor(aes(color = year), label.x = 2) 

Corr_plot_ramif


ggarrange(Corr_plot_ramif,Corr_plot_rdt, 
  nrow = 2, common.legend = TRUE
) 
