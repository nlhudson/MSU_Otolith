
#### Otolith Data from FWP 2024 Spring Electrofishing 

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

benny_data = read_xlsx("MFWP.Electrofishing.2024.xlsx")
view(benny_data)

benny_data_ll_rb = benny_data[benny_data$Species %in% c("LL", "RB"),]

# Data Summary 

fwp_2024_summary_full <- benny_data %>% 
  group_by(Waterbody, Species) %>%  
  summarize(sample_size = n(),
            avg_mass_g = mean(Weight, na.rm=TRUE),
            SE_mass_g = sd(Weight, na.rm=TRUE)/sqrt(n()),
            avg_length_mm = mean(Length, na.rm=TRUE),
            SE_length_mm = sd(Weight, na.rm=TRUE)/sqrt(n())
            )

fwp_2024_summary_ll_rb <- benny_data_ll_rb %>% 
  group_by(Waterbody, Species) %>%  
  summarize(sample_size = n(),
            avg_mass_g = mean(Weight, na.rm=TRUE),
            SE_mass_g = sd(Weight, na.rm=TRUE)/sqrt(n()),
            avg_length_mm = mean(Length, na.rm=TRUE),
            SE_length_mm = sd(Weight, na.rm=TRUE)/sqrt(n())
  )



## mass by river (450g = 1lb)

# all species
ggplot(benny_data, aes(x = Waterbody, y = Weight, color = Species)) + 
  geom_point(position = position_dodge(1), alpha = 0.2) +
  stat_summary(aes(color = Species, width = 0.2), fun.data="mean_se",  fun.args = list(mult=1),geom = "errorbar", size = 1,  position = position_dodge(1)) +
  stat_summary(aes(group = Species, color = Species, shape = Location), geom = "point",fun.args = list(mult=1), fun = mean, shape =16, size = 5 , position = position_dodge(1)) + 
  labs(y = "Mass (g)")

# brown and rainbows
ggplot(benny_data_ll_rb, aes(x = Waterbody, y = Weight, color = Species)) + 
  geom_point(position = position_dodge(1), alpha = 0.2) +
  stat_summary(aes(color = Species, width = 0.2), fun.data="mean_se",  fun.args = list(mult=1),geom = "errorbar", size = 1,  position = position_dodge(1)) +
  stat_summary(aes(group = Species, color = Species, shape = Location), geom = "point",fun.args = list(mult=1), fun = mean, shape =16, size = 5 , position = position_dodge(1)) + 
  labs(y = "Mass (g)")

## length by river (25mm = 1in ; 300mm = 12in)

# all species
ggplot(benny_data, aes(x = Waterbody, y = Length, color = Species)) + 
  geom_point(position = position_dodge(1), alpha = 0.2) +
  stat_summary(aes(color = Species, width = 0.2), fun.data="mean_se",  fun.args = list(mult=1),geom = "errorbar", size = 1,  position = position_dodge(1)) +
  stat_summary(aes(group = Species, color = Species, shape = Location), geom = "point",fun.args = list(mult=1), fun = mean, shape =16, size = 5 , position = position_dodge(1)) + 
  labs(y = "Length (mm)")

# brown and rainbows
ggplot(benny_data_ll_rb, aes(x = Waterbody, y = Length, color = Species)) + 
  geom_point(position = position_dodge(1), alpha = 0.2) +
  stat_summary(aes(color = Species, width = 0.2), fun.data="mean_se",  fun.args = list(mult=1),geom = "errorbar", size = 1,  position = position_dodge(1)) +
  stat_summary(aes(group = Species, color = Species, shape = Location), geom = "point",fun.args = list(mult=1), fun = mean, shape =16, size = 5 , position = position_dodge(1)) + 
  labs(y = "Length (mm)")

## Mass length relationship by river

# all species
ggplot(benny_data, aes(x = Length, y = Weight, color = Species)) + 
  geom_point(alpha = 0.6) +
facet_grid(Waterbody ~ Species) +
  labs(x = "Length (mm)", y = "Weight (g)")

# brown and rainbows
ggplot(benny_data_ll_rb, aes(x = Length, y = Weight, color = Species)) + 
  geom_point(alpha = 0.6) +
  facet_grid(Waterbody ~ Species) +
  labs(x = "Length (mm)", y = "Weight (g)")

## Linear model for mass length relationship by species 

# all species
mass_length = lm(log(Weight) ~ log(Length) + Species, data = benny_data)
summary(mass_length)

ggplot(benny_data, aes(x = log(Length), y = log(Weight), color = Species)) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method='lm', color = "black") +
  facet_grid(Waterbody ~ Species) +
  labs(x = "Log Length (mm) ", y = "Log Weight (g)")

# brown and rainbows
mass_length_ll_rb = lm(log(Weight) ~ log(Length) + Species, data = benny_data_ll_rb)
summary(mass_length_ll_rb)

ggplot(benny_data, aes(x = log(Length), y = log(Weight), color = Species)) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method='lm', color = "black") +
  facet_grid(Waterbody ~ Species) +
  labs(x = "Log Length (mm) ", y = "Log Weight (g)")





