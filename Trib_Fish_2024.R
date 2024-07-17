

#### Looking at trib fish data from preliminary shocking in July 2024

# Loading any relevant packages
library(car)
library(multcompView)
library(multcomp)
library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(emmeans)
library(stringr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(knitr)
library(readxl)
library(outliers)
library(psych)
library(vtable)
library(lattice)
library(magrittr)
library(mgcv)
library(grid)


trib_data = read_xlsx("Trib Fish Data.xlsx")
view(trib_data)

trib_data_BRN_RBT = trib_data[trib_data$Species %in% c("BRN", "RBT"),]

# Data Summary 

trib_data_sum <- trib_data_BRN_RBT %>% 
  group_by(Waterbody, Species) %>%  
  summarize(sample_size = n(),
            avg_mass_g = mean(`Mass_(g)`, na.rm=TRUE),
            SE_mass_g = sd(`Mass_(g)`, na.rm=TRUE)/sqrt(n()),
            avg_length_mm = mean(`Length_(mm)`, na.rm=TRUE),
            SE_length_mm = sd(`Length_(mm)`, na.rm=TRUE)/sqrt(n())
  )

sample_sizes <- trib_data_BRN_RBT %>%
  group_by(Waterbody, Species) %>%
  summarise(n = n()) %>%
  ungroup()

ggplot(trib_data_BRN_RBT, aes(x = Waterbody, y = `Length_(mm)`, color = Species)) + 
  geom_point(position = position_dodge(1), alpha = 0.4) +
  stat_summary(aes(color = Species, width = 0.2), fun.data="mean_se",  fun.args = list(mult=1),geom = "errorbar", size = 1,  position = position_dodge(1)) +
  stat_summary(aes(group = Species, color = Species, shape = Location), geom = "point",fun.args = list(mult=1), fun = mean, shape =16, size = 4 , position = position_dodge(1)) + 
  labs(y = "Length (mm)") +
  geom_text(data = sample_sizes, aes(x = Waterbody, y = max(trib_data_BRN_RBT$`Length_(mm)`) + 10, label = n), 
            position = position_dodge(1), vjust = 10)


