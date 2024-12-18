
#### Redd data analysis

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

redd_dat = read.csv("Redd_data.csv")

redd_dat$Redds_by_increment = as.numeric(redd_dat$Redds_by_increment)
redd_dat$Total_redds = as.numeric(redd_dat$Total_redds)
redd_dat$Active_Spawners = as.numeric(redd_dat$Active_Spawners)
redd_dat$Reader = as.factor(redd_dat$Reader)
redd_dat$Site = as.factor(redd_dat$Site)

redd_dat_1 = redd_dat[redd_dat$Reader %in% c("1"),]

redd_dat_sum = redd_dat %>%
  group_by(Site, Reader) %>%
  summarize(total_reach_length = max(Reach_Distance)/1000,
            total_redds = sum(Redds_by_increment)) %>%
  mutate(redds_per_km = total_redds/total_reach_length)

redds = 
ggplot(redd_dat_sum, aes(x=factor(Site, level=c('Cherry', 'Trapper', 'Moose', 'Canyon', 'Camp','Deep', 
                                               'Cox', 'Blacktail','Grasshopper',
                                               'Wisconsin', 'Mill', 'Silver','Alder')), 
                                        y = redds_per_km, fill = as.factor(Reader))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_grey(base_size = 20) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = -45, hjust = -0.05)) +
  labs(y = "Redds/km", fill = "Reader")
  
ggsave(redds, file = "redds.jpg", dpi = 600)

#### Recruit per Spawner ####

depletion_dat = read.csv("depletion_dat.csv") # depletion data

# blacktain had no yoy, camp was not sampleable, deep had no depletion estimate

redd_dat_sum_1 = redd_dat_sum[redd_dat_sum$Reader %in% c("1"),] # removing blacktail
redd_dat_sum_1 = redd_dat_sum_1[-2,]
redd_dat_sum_1 = redd_dat_sum_1[-5,]

depletion_dat = depletion_dat[-2,] # removing camp



view(redd_dat_sum_1)
view(depletion_dat)

sr = merge(redd_dat_sum_1, depletion_dat, by = "Site")

view(sr)

sr_plot = 
ggplot(sr, aes(x=redds_per_km, y = Np, color = Site)) +
  geom_point(size = 3) +
  theme(text=element_text(size = 20)) +
  labs(y = "Recruits", x = "Spawners")

ggsave(sr_plot, file = "sr_plot.jpg", dpi = 600)





