
### 2024 Otolith collection length/weight distributions

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

BH_otos = read.csv("BH_Otos_24.csv")
BV_otos = read.csv("BV_Otos_24.csv")
RB_otos = read.csv("RB_Otos_24.csv")

# check numeric
BH_otos$Length_mm = as.numeric(BH_otos$Length_mm)
BV_otos$Length_mm = as.numeric(BV_otos$Length_mm)
RB_otos$Length_mm = as.numeric(RB_otos$Length_mm)

BH_otos_BRN = BH_otos[BH_otos$Species %in% c("BRN"),]
BH_otos_RBT = BH_otos[BH_otos$Species %in% c("RBT"),]

 # summary tables for counts 
BH_otos <- BH_otos %>%
  mutate(Length_Range = cut(
    Length_mm,
    breaks = c(0, 305, 405, 5000),
    labels = c("0-305", "306-405", "406-5000"),
    right = TRUE
  ))

BH_sum <- BH_otos %>%
  group_by(Species,Year) %>%
  summarize(Catch = n(), .groups = "drop")

BV_otos <- BV_otos %>%
  mutate(Length_Range = cut(
    Length_mm,
    breaks = c(0, 305, 405, 5000),
    labels = c("0-305", "306-405", "406-5000"),
    right = TRUE
  ))

BV_sum <- BV_otos %>%
  group_by(Species, Year) %>%
  summarize(Catch = n(), .groups = "drop")

RB_otos <- RB_otos %>%
  mutate(Length_Range = cut(
    Length_mm,
    breaks = c(0, 305, 405, 5000),
    labels = c("0-305", "306-405", "406-5000"),
    right = TRUE
  ))

RB_sum <- RB_otos %>%
  group_by(Species, Year) %>%
  summarize(Catch = n(), .groups = "drop")

view(RB_sum)

# Length histogram for BH BRN - match the bin width with the break width for better visualization

BH_BRN_hist = 
  ggplot(BH_otos_BRN, aes(x = Length_mm, color = Year)) +
  geom_histogram(color = "black", fill = "white", binwidth = 20, boundary = 0, closed = "left") + 
  scale_x_continuous(breaks = seq(0,600,20)) +
  scale_y_continuous(breaks = seq(0,20,2)) +
  geom_vline(xintercept=305, color = "red") +
  geom_vline(xintercept=406, color = "red") +
  theme(text = element_text(size = 20))

BH_BRN_hist2 = 
  ggplot(BH_otos_BRN, aes(x = Length_mm, fill = as.factor(Year))) +
    geom_histogram(position = "dodge", binwidth = 20, boundary = 0, closed = "left", color = "black") + 
    scale_x_continuous(breaks = seq(0, 600, 50)) +
    scale_y_continuous(breaks = seq(0, 20, 2)) +
    geom_vline(xintercept = 305, color = "red") +
    geom_vline(xintercept = 406, color = "red") +  
    theme(text = element_text(size = 20)) +
    labs(fill = "Year")

ggsave(BH_BRN_hist2, file = "BH_BRN_hist.jpg", dpi = 600)

# Length histogram for BH RBT

BH_RBT_hist = 
  ggplot(BH_otos_RBT, aes(x = Length_mm)) +
  geom_histogram(color = "black", fill = "white", binwidth = 20, boundary = 0, closed = "left") + 
  scale_x_continuous(breaks = seq(0,600,50)) +
  scale_y_continuous(breaks = seq(0,20,2)) +
  geom_vline(xintercept=305, color = "red") +
  geom_vline(xintercept=406, color = "red") +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(0.1,1,0.1,0.1), "cm")) +
  

ggsave(BH_RBT_hist2, file = "BH_RBT_hist.jpg", dpi = 600)


# Length histogram for BV BRN
BV_BRN_hist = 
  ggplot(BV_otos, aes(x = Length_mm)) +
  geom_histogram(color = "black", fill = "white", binwidth = 20, boundary = 0, closed = "left") + 
  scale_x_continuous(breaks = seq(0,600,20)) +
  scale_y_continuous(breaks = seq(0,20,2)) +
  geom_vline(xintercept=305, color = "red") +
  geom_vline(xintercept=406, color = "red") +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(0.1,1,0.1,0.1), "cm"))

BV_BRN_hist2 = 
  ggplot(BV_otos, aes(x = Length_mm, fill = as.factor(Year))) +
  geom_histogram(position = "dodge", binwidth = 20, boundary = 0, closed = "left", color = "black") + 
  scale_x_continuous(breaks = seq(0, 600, 50)) +
  scale_y_continuous(breaks = seq(0, 20, 2)) +
  geom_vline(xintercept = 305, color = "red") +
  geom_vline(xintercept = 406, color = "red") + # Apply color palette
  theme(text = element_text(size = 20)) +
  labs(fill = "Year")

ggsave(BV_BRN_hist2, file = "BV_BRN_hist.jpg", dpi = 600)


# Length histogram for RB BRN
RB_BRN_hist = 
ggplot(RB_otos, aes(x = Length_mm)) +
  geom_histogram(color = "black", fill = "white", binwidth = 20, boundary = 0, closed = "left") + 
  scale_x_continuous(breaks = seq(0,600,20)) +
  scale_y_continuous(breaks = seq(0,5,1)) +
  geom_vline(xintercept=305, color = "red") +
  geom_vline(xintercept=406, color = "red") +
  theme(text = element_text(size = 20))

RB_BRN_hist2 = 
  ggplot(RB_otos, aes(x = Length_mm, fill = as.factor(Year))) +
  geom_histogram(position = "dodge", binwidth = 20, boundary = 0, closed = "left", color = "black") + 
  scale_x_continuous(breaks = seq(0, 600, 50)) +
  scale_y_continuous(breaks = seq(0, 20, 2)) +
  geom_vline(xintercept = 305, color = "red") +
  geom_vline(xintercept = 406, color = "red") + # Apply color palette
  theme(text = element_text(size = 20)) +
  labs(fill = "Year")

ggsave(RB_BRN_hist2, file = "RB_BRN_hist.jpg", dpi = 600)


