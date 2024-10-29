
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

BH_otos = read_xlsx("BH_Otos_24.xlsx")
BV_otos = read_xlsx("BV_Otos_24.xlsx")
RB_otos = read_xlsx("RB_Otos_24.xlsx")

# check numeric
class(BH_otos$Length_mm)
class(BV_otos$Length_mm)
class(RB_otos$Length_mm)

BH_otos_BRN = BH_otos[BH_otos$Species %in% c("BRN"),]
BH_otos_RBT = BH_otos[BH_otos$Species %in% c("RBT"),]

# Length histogram for BH BRN - match the bin width with the break width for better visualization

BH_BRN_hist = 
  ggplot(BH_otos_BRN, aes(x = Length_mm)) +
  geom_histogram(color = "black", fill = "white", binwidth = 20, boundary = 0, closed = "left") + 
  scale_x_continuous(breaks = seq(0,600,20)) +
  scale_y_continuous(breaks = seq(0,20,2)) +
  geom_vline(xintercept=305, color = "red") +
  geom_vline(xintercept=406, color = "red") +
  theme(text = element_text(size = 20))

ggsave(BH_BRN_hist, file = "BH_BRN_hist.jpg", dpi = 600)

# Length histogram for BH RBT

BH_RBT_hist = 
  ggplot(BH_otos_RBT, aes(x = Length_mm)) +
  geom_histogram(color = "black", fill = "white", binwidth = 20, boundary = 0, closed = "left") + 
  scale_x_continuous(breaks = seq(0,600,20)) +
  scale_y_continuous(breaks = seq(0,20,2)) +
  geom_vline(xintercept=305, color = "red") +
  geom_vline(xintercept=406, color = "red") +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(0.1,1,0.1,0.1), "cm"))

ggsave(BH_RBT_hist, file = "BH_RBT_hist.jpg", dpi = 600)


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

ggsave(BV_BRN_hist, file = "BV_BRN_hist.jpg", dpi = 600)


# Length histogram for RB BRN
RB_BRN_hist = 
ggplot(RB_otos, aes(x = Length_mm)) +
  geom_histogram(color = "black", fill = "white", binwidth = 20, boundary = 0, closed = "left") + 
  scale_x_continuous(breaks = seq(0,600,20)) +
  scale_y_continuous(breaks = seq(0,5,1)) +
  geom_vline(xintercept=305, color = "red") +
  geom_vline(xintercept=406, color = "red") +
  theme(text = element_text(size = 20))

ggsave(RB_BRN_hist, file = "RB_BRN_hist.jpg", dpi = 600)


