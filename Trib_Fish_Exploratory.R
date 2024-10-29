

#### Looking at trib fish data from 2024

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
library(writexl)
library(outliers)
library(psych)
library(vtable)
library(lattice)
library(magrittr)
library(mgcv)
library(grid)


trib_data = read.csv("2024 Trib Fish Data.csv")
trib_data$Length_mm = as.numeric(trib_data$Length_mm)
trib_data$Mass_g = as.numeric(trib_data$Mass_g)

trib_data_BRN_RBT = trib_data[trib_data$Species %in% c("BRN", "RBT"),]

trib_data_BRN = trib_data[trib_data$Species %in% c("BRN"),]

trib_data_RBT = trib_data[trib_data$Species %in% c("RBT"),]

trib_data_BRN_pass1 = trib_data_BRN[trib_data_BRN$Method %in% "Pass 1",] # only pass 1 data 

trib_data_RBT_pass1 = trib_data_RBT[trib_data_RBT$Method %in% "Pass 1",] # only pass 1 data 

# Average Length and Weight (not useful until I fill in the missing weights)
trib_data_sum <- trib_data_BRN_RBT %>% 
  group_by(Waterbody, Species) %>%  
  summarize(sample_size = n(),
            avg_mass_g = mean(Mass_g, na.rm=TRUE),
            SE_mass_g = sd(Mass_g, na.rm=TRUE)/sqrt(n()),
            avg_length_mm = mean(Length_mm, na.rm=TRUE),
            SE_length_mm = sd(Length_mm, na.rm=TRUE)/sqrt(n())
  )

# Sample sizes for different tribs, species, age classes - useful for 1st pass CPUE

# for BRN less than 140 is YOY 
trib_sample_sizes_BRN <- trib_data_BRN_pass1 %>%
  mutate(Age_Class = ifelse(Length_mm >= 140, "Other", "Juvenile")) %>%
  group_by(Basin, Waterbody, Species, Age_Class) %>%
  summarise(n = n()) %>%
  ungroup()

# for RBT less than 90 is YOY 
trib_sample_sizes_RBT <- trib_data_RBT_pass1 %>%
  mutate(Age_Class = ifelse(Length_mm >= 90, "Other", "Juvenile")) %>%
  group_by(Basin, Waterbody, Species, Age_Class) %>%
  summarise(n = n()) %>%
  ungroup()

trib_sample_sizes_BRN_RBT <- bind_rows(
  trib_sample_sizes_BRN, 
  trib_sample_sizes_RBT
)

write_xlsx(trib_sample_sizes_BRN_RBT, "pass1_sample_sizes.xlsx") 
# above writes xlsx for 1st pass data, then I added reach length and calculated CPUE and renamed the file "Trib_Pass1_CPUE.xlsx"

# Plots for CPUE between sites 
Pass1_CPUE = read_xlsx("Trib_Pass1_CPUE.xlsx")

Pass1_CPUE_BH = Pass1_CPUE[Pass1_CPUE$Basin == "Big Hole",]
Pass1_CPUE_BV = Pass1_CPUE[Pass1_CPUE$Basin == "Beaverhead",]
Pass1_CPUE_RB = Pass1_CPUE[Pass1_CPUE$Basin == "Ruby",]

# jittered plots to show CPUE by age class and species 
BH_trib_pass1_cpue=
  ggplot(Pass1_CPUE_BH, aes(x = Waterbody, y = CPUE, fill = Age_Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "CPUE (fish/meter)") +
  theme(text = (text=element_text(size=20))) +
  scale_y_continuous(breaks = seq(0,2.0,0.2), limits = c(0, 2.0)) +
  facet_wrap(~Species)

BV_trib_pass1_cpue=
ggplot(Pass1_CPUE_BV, aes(x = Waterbody, y = CPUE, fill = Age_Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "CPUE (fish/meter)") +
  scale_y_continuous(breaks = seq(0,1.0,0.2), limits = c(0, 1.0)) +
  theme(text = (text=element_text(size=20))) +
  facet_wrap(~Species)

RB_trib_pass1_cpue=
ggplot(Pass1_CPUE_RB, aes(x = Waterbody, y = CPUE, fill = Age_Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "CPUE (fish/meter)") +
  theme(text = (text=element_text(size=20))) +
  scale_y_continuous(breaks = seq(0,2.0,0.2), , limits = c(0, 1.6)) +
  facet_wrap(~Species)

ggsave(BH_trib_pass1_cpue, file = "BH_trib_pass1_cpue.jpg", dpi = 600)
ggsave(BV_trib_pass1_cpue, file = "BV_trib_pass1_cpue.jpg", dpi = 600)
ggsave(RB_trib_pass1_cpue, file = "RB_trib_pass1_cpue.jpg", dpi = 600)


## Exploratory plots
ggplot(trib_data_BRN_RBT, aes(x = Waterbody, y = Length_mm, color = Species)) + 
  geom_point(position = position_dodge(1), alpha = 0.4) +
  stat_summary(aes(color = Species, width = 0.2), fun.data="mean_se",  fun.args = list(mult=1),geom = "errorbar", size = 1,  position = position_dodge(1)) +
  stat_summary(aes(group = Species, color = Species, shape = Location), geom = "point",fun.args = list(mult=1), fun = mean, shape =16, size = 4 , position = position_dodge(1)) + 
  labs(y = "Length (mm)") +
  geom_text(data = sample_sizes, aes(x = Waterbody, y = max(trib_data_BRN_RBT$`Length (mm)`) + 10, label = n), 
            position = position_dodge(1), vjust = 10)

ggplot(trib_data_BRN_RBT, aes(x = log(Length_mm), y = log(Mass_g), color = Waterbody)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)


# normal linear model for length predicting mass by trib
trib_fish_lm = lm(log(Mass_g) ~ log(Length_mm)*Waterbody, data = trib_data_BRN)
summary(trib_fish_lm)

# linear mixed effects model for length predicting mass by trib, trib as random effect 
trib_fish_lmer = lmer(log(Mass_g) ~ log(Length_mm) + (1|Waterbody), data = trib_data_BRN)
summary(trib_fish_lmer)

## Looking at just trapper, cherry, Mill, and upper moose
trib_data_BRN_select = trib_data_BRN[trib_data_BRN$Waterbody %in% c("Trapper - Lower", "Cherry Creek", "Mill", "Upper Moose Creek"),]

ggplot(trib_data_BRN_select, aes(x = Length_mm, color = Waterbody)) +
  geom_histogram(fill="white", alpha=0.5, position="dodge", size = 3)

ggplot(trib_data_BRN_select, aes(x = Length_mm)) +
  geom_histogram(color = "black", fill ="white") +
  scale_x_continuous(breaks = seq(0,350,25)) +
  scale_y_continuous(breaks = seq(0,150,25)) +
  facet_wrap(~Waterbody)

ggplot(trib_data_BRN_select, aes(x = log(Length_mm), y = log(Mass_g), color = Waterbody)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

ggplot(trib_data_BRN_select[trib_data_BRN_select$Waterbody == "Trapper - Lower", ], aes(x = Length_mm, y = Mass_g, color = Waterbody)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

summary(lm(log(Mass_g) ~ log(Length_mm)*Waterbody, data = trib_data_BRN_select))


