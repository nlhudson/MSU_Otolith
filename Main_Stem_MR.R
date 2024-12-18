

### Main Stem-Mark Recap Analysis

library(car)
library(multcompView)
library(multcomp)
library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(knitr)
library(readxl)
library(writexl)
library(magrittr)
library(mgcv)
library(grid)
library(FSA)

main_stem_mr_full = read.csv("Main_Stem_MR.csv")
main_stem_mr_full$Length_mm = as.numeric(main_stem_mr_full$Length_mm)
main_stem_mr_full$Mass_g = as.numeric(main_stem_mr_full$Mass_g)
main_stem_mr_full$Site = as.factor(main_stem_mr_full$Site)
main_stem_mr_full$Basin = as.factor(main_stem_mr_full$Basin)
main_stem_mr_full$Pass = as.factor(main_stem_mr_full$Pass)

mr_dat = main_stem_mr_full %>%
  group_by(Basin, Site, Pass, Species, Reach.Length_mm) %>%
  summarize(Catch = n(),
  Recaps = sum(Recap == "R" & Pass == 2, na.rm = TRUE)) 

mr_results <- mr_dat %>%
  group_by(Basin, Site, Species, Reach.Length_mm) %>%
  summarize(
    m = sum(Catch[Pass == 1], na.rm = TRUE), # Total Catch in Pass 1
    n = sum(Catch[Pass == 2], na.rm = TRUE), # Total Catch in Pass 2
    r = sum(Recaps[Pass == 2], na.rm = TRUE) # Recaps in Pass 2
  ) %>%
  filter(r>0) %>%
  group_by(Basin, Site, Species) %>%
  mutate(N = (mrClosed(m,n,r)$N)*2, 
         uci = (confint(mrClosed(m,n,r))[,2])*2,
         lci = (confint(mrClosed(m,n,r))[,1])*2) %>%
  mutate(yoy_per_km = N/(Reach.Length_mm/1000), # correcting N for reach length
         uci_std = uci/(Reach.Length_mm/1000), # correcting uci and lci for reach length 
         lci_std = lci/(Reach.Length_mm/1000))

mr_results_bh = mr_results[mr_results$Basin %in% c("Big Hole"),]
mr_results_bv = mr_results[mr_results$Basin %in% c("Beaverhead"),]
mr_results_rb = mr_results[mr_results$Basin %in% c("Ruby"),]

mr = 
  ggplot(mr_results, aes(x=factor(Site, levels = c('Browns Bridge', 'Salmon Fly', 'George Grant', 'Jerry Creek',
                                                   'Ten Mile', 'Poindexter', 'Pipe Organ', 'Hildreth',
                                                   'Silver Springs', 'Barnosky', 'Alder', 'Vigilante')), y = yoy_per_km, color = Species)) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = lci_std, ymax = uci_std, width = 0.2, linetype = "95% CI",), 
                position = position_dodge(0.5), size = 1) +
  theme_grey(base_size = 20) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = -45, hjust = -0.05)) +
  scale_linetype_manual(values = c("95% CI" = "solid"), name = "Error Bars") +
  labs(y = "YOY / km", x = "Site")

ggsave(mr, file = "mr.jpg", dpi = 600)


mr_bh = 
ggplot(mr_results_bh, aes(x=factor(Site, levels = c('Browns Bridge', 'Salmon Fly', 'George Grant', 'Jerry Creek')), y = N, color = Species)) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = lci, ymax = uci, width = 0.2, linetype = "95% CI",), 
                position = position_dodge(0.5), size = 1) +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05)) +
  scale_linetype_manual(values = c("95% CI" = "solid"), name = "Error Bars") +
  labs(y = "Np", x = "Site")

mr_bv = 
ggplot(mr_results_bv, aes(x=factor(Site, levels = c('Ten Mile', 'Poindexter', 'Pipe Organ', 'Hildreth')), y = N, color = Species)) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = lci, ymax = uci, width = 0.2, linetype = "95% CI",), 
                position = position_dodge(0.5), size = 1) +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05)) +
  scale_linetype_manual(values = c("95% CI" = "solid"), name = "Error Bars") +
  labs(y = "Np", x = "Site")

mr_rb = 
  ggplot(mr_results_bv, aes(x=factor(Site, levels = c('Silver Springs', 'Barnosky', 'Alder', 'Vigilante')), y = N, color = Species)) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = lci, ymax = uci, width = 0.2, linetype = "95% CI",), 
                position = position_dodge(0.5), size = 1) +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05)) +
  scale_linetype_manual(values = c("95% CI" = "solid"), name = "Error Bars") +
  labs(y = "Np", x = "Site")
  

