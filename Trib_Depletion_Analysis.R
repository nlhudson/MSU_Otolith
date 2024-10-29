
###### Tributary Depletion Data Analysis

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
library(FSA)

###### Loading Data ######

trib_fish_data_full = read.csv("2024 Trib Fish Data.csv")
trib_fish_data_full$Length_mm = as.numeric(trib_data$Length_mm)
trib_fish_data_full$Mass_g = as.numeric(trib_data$Mass_g)

depletion_data = trib_fish_data_full[trib_fish_data_full$Method %in% c("Pass 1", "Pass 2", "Pass 3"),]
depletion_data_onlytrout = depletion_data[depletion_data$Species %in% c("BRN", "RBT", "BRK"),]

reach_lengths <- c("Alder" = 200, "Blacktail" = 250, "Camp" = 200, "Canyon" = 200, "Cherry" = 250, "Cox Slough" = 200,
                   "Deep" = 230, "Grasshopper" = 250, "Mill" = 200, "Moose" = 200, "Silver Springs" = 200,
                   "Stone" = 200, "Trapper" = 250, "Wisconsin" = 250)

by_pass = depletion_data_onlytrout %>%
  group_by(Method, Site, Species) %>%
  summarize(Catch = n()) %>%
  mutate(Reach_Length_m = reach_lengths[Site])


#### Leslie ####

leslie_dat1 = by_pass %>%
  group_by(Site, Species, Method, Catch, Reach_Length_m) %>%
  reframe(cpe = Catch/Reach_Length_m) # calculate catch per effort vector

leslie_dat2 = leslie_dat1 %>%
  group_by(Site, Species) %>%
  reframe(K = cumsum(Catch)-Catch) # calculating cumulative catch (K)

leslie_dat1$K = leslie_dat2$K # combining K with cpe onto same dataframe

# summary figure
ggplot(leslie_dat1, aes(x=K,y=cpe, color = Species)) +
  geom_point() + 
  geom_smooth(method='lm') +
  facet_wrap(~Site)

# Alder
leslie_alder = leslie_dat1[leslie_dat1$Site %in% "Alder",] # BRN only
leslie_model_alder <- with(leslie_alder,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_alder)
confint(leslie_model_alder)
plot(leslie_model_alder)

# Blacktail
leslie_blacktail_BRN = subset(leslie_dat1, Site == "Blacktail" & Species == "BRN") # BRN
leslie_model_blacktail_BRN <- with(leslie_blacktail_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_blacktail_BRN)
confint(leslie_model_blacktail_BRN)
plot(leslie_model_blacktail_BRN)

leslie_blacktail_BRK = subset(leslie_dat1, Site == "Blacktail" & Species == "BRK") # BRK
leslie_model_blacktail_BRK <- with(leslie_blacktail_BRK,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_blacktail_BRK)
confint(leslie_model_blacktail_BRK)
plot(leslie_model_blacktail_BRK)

# Camp
leslie_camp_BRN = subset(leslie_dat1, Site == "Camp" & Species == "BRN") # BRN
leslie_model_camp_BRN <- with(leslie_camp_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_camp_BRN)
confint(leslie_model_camp_BRN)
plot(leslie_model_camp_BRN)

leslie_camp_BRK = subset(leslie_dat1, Site == "Camp" & Species == "BRK") # BRK
leslie_model_camp_BRK <- with(leslie_camp_BRK,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_camp_BRK)
confint(leslie_model_camp_BRK)
plot(leslie_model_camp_BRK)

# Canyon
leslie_canyon_BRN = subset(leslie_dat1, Site == "Canyon" & Species == "BRN") # BRN
leslie_model_canyon_BRN <- with(leslie_canyon_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_canyon_BRN)
confint(leslie_model_canyon_BRN)
plot(leslie_model_canyon_BRN)

leslie_canyon_RBT = subset(leslie_dat1, Site == "Canyon" & Species == "RBT") # RBT
leslie_model_canyon_RBT <- with(leslie_canyon_RBT,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_canyon_RBT)
confint(leslie_model_canyon_RBT)
plot(leslie_model_canyon_RBT)

leslie_canyon_BRK = subset(leslie_dat1, Site == "Canyon" & Species == "BRK") # BRK
leslie_model_canyon_BRK <- with(leslie_canyon_BRK,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_canyon_BRK)
confint(leslie_model_canyon_BRK)
plot(leslie_model_canyon_BRK)

# Cherry
leslie_cherry_BRN = leslie_dat1[leslie_dat1$Site %in% "Cherry",] # BRN only
leslie_model_cherry_BRN <- with(leslie_cherry_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_cherry_BRN)
confint(leslie_model_cherry_BRN)
plot(leslie_model_cherry_BRN)

# Cox
leslie_cox_BRN = leslie_dat1[leslie_dat1$Site %in% "Cox Slough",] # BRN only
leslie_model_cox_BRN <- with(leslie_cox,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_cox_BRN)
confint(leslie_model_cox_BRN)
plot(leslie_model_cox_BRN)

# Deep
leslie_deep_BRN = subset(leslie_dat1, Site == "Deep" & Species == "BRN") # BRN
leslie_model_deep_BRN <- with(leslie_deep_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE)) 
summary(leslie_model_deep_BRN)
confint(leslie_model_deep_BRN) # only 2 passe with catch, cant work
plot(leslie_model_deep_BRN)

leslie_deep_RBT = subset(leslie_dat1, Site == "Deep" & Species == "RBT") # RBT
leslie_model_deep_RBT <- with(leslie_deep_RBT,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_deep_RBT)
confint(leslie_model_deep_RBT)
plot(leslie_model_deep_RBT)

leslie_deep_BRK = subset(leslie_dat1, Site == "Deep" & Species == "BRK") # BRK
leslie_model_deep_BRK <- with(leslie_deep_BRK,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_deep_BRK)
confint(leslie_model_deep_BRK) # only 2 passe with catch, cant work
plot(leslie_model_deep_BRK)

# Grasshopper
leslie_grasshopper_BRN = subset(leslie_dat1, Site == "Grasshopper" & Species == "BRN") # BRN
leslie_model_grasshopper_BRN <- with(leslie_grasshopper_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_grasshopper_BRN)
confint(leslie_model_grasshopper_BRN)
plot(leslie_model_grasshopper_BRN)

leslie_grasshopper_BRK = subset(leslie_dat1, Site == "Grasshopper" & Species == "BRK") # BRK
leslie_model_grasshopper_BRK <- with(leslie_grasshopper_BRK,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_grasshopper_BRK)
confint(leslie_model_grasshopper_BRK)
plot(leslie_model_grasshopper_BRK)

# Mill
leslie_mill_BRN = leslie_dat1[leslie_dat1$Site %in% "Mill",] # BRN only
leslie_model_mill_BRN <- with(leslie_mill_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_mill_BRN)
confint(leslie_model_mill_BRN)
plot(leslie_model_mill_BRN)

# Moose
leslie_moose_BRN = subset(leslie_dat1, Site == "Moose" & Species == "BRN") # BRN
leslie_model_moose_BRN <- with(leslie_moose_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_moose_BRN)
confint(leslie_model_moose_BRN)
plot(leslie_model_moose_BRN)

leslie_moose_RBT = subset(leslie_dat1, Site == "Moose" & Species == "RBT") # RBT
leslie_model_moose_RBT <- with(leslie_moose_RBT,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_moose_RBT)
confint(leslie_model_moose_RBT)
plot(leslie_model_moose_RBT)

# Silver Springs
leslie_silver_BRN = subset(leslie_dat1, Site == "Silver Springs" & Species == "BRN") # BRN
leslie_model_silver_BRN <- with(leslie_silver_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_silver_BRN)
confint(leslie_model_silver_BRN)
plot(leslie_model_silver_BRN)

leslie_silver_RBT = subset(leslie_dat1, Site == "Silver Springs" & Species == "RBT") # RBT
leslie_model_silver_RBT <- with(leslie_silver_RBT,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_silver_RBT)
confint(leslie_model_silver_RBT)
plot(leslie_model_silver_RBT)

# Stone
leslie_stone_BRN = leslie_dat1[leslie_dat1$Site %in% "Stone",] # BRN only
leslie_model_stone_BRN <- with(leslie_stone_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE)) # only 2 passes, doesnt work
summary(leslie_model_stone)
confint(leslie_model_stone)
plot(leslie_model_stone)

# Trapper
leslie_trapper_BRN = subset(leslie_dat1, Site == "Trapper" & Species == "BRN") # BRN only
leslie_model_trapper_BRN <- with(leslie_trapper_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_trapper_BRN)
confint(leslie_model_trapper_BRN)
plot(leslie_model_trapper_BRN)

# Wisconsin
leslie_wisconsin_BRN = leslie_dat1[leslie_dat1$Site %in% "Wisconsin",] # BRN only
leslie_model_wisconsin_BRN <- with(leslie_wisconsin_BRN,depletion(Catch,Reach_Length_m, Ricker.mod=TRUE))
summary(leslie_model_wisconsin_BRN)
confint(leslie_model_wisconsin_BRN)
plot(leslie_model_wisconsin_BRN)

### Leslie needs a miniumum of 3 passes to work, improves with linear depletion and more passes 
### More linear depletion reduces error
### Now try K Pass

#### K Pass ####

# Alder 
k_alder_BRN = subset(by_pass, Site == "Alder" & Species == "BRN") # BRN
k_alder_est_BRN = removal(k_alder_BRN$Catch)
summary(k_alder_BRN_est)
confint(k_alder_BRN_est) # lower se than leslie

# Blacktail 
k_blacktail_BRN = subset(by_pass, Site == "Blacktail" & Species == "BRN") # BRN
k_blacktail_BRN_est = removal(k_blacktail_BRN$Catch)
summary(k_blacktail_BRN_est)
confint(k_blacktail_BRN_est) # lower se than leslie

k_blacktail_BRK = subset(by_pass, Site == "Blacktail" & Species == "BRK") # BRK - cannot calculate with only 1 pass of catch data
k_blacktail_BRK_est = removal(k_blacktail_BRK$Catch)
summary(k_blacktail_BRK_est)
confint(k_blacktail_BRK_est) # cannot calculate with only 1 pass of catch data

# Camp
k_camp_BRN = subset(by_pass, Site == "Camp" & Species == "BRN") # BRN
k_camp_BRN_est = removal(k_camp_BRN$Catch)
summary(k_camp_BRN_est)
confint(k_camp_BRN_est) # lower se than leslie

k_camp_BRK = subset(by_pass, Site == "Camp" & Species == "BRK") # RBT
k_camp_BRK_est = removal(k_camp_BRK$Catch)
summary(k_camp_BRK_est)
confint(k_camp_BRK_est) # lower se than leslie

# Canyon
k_canyon_BRN = subset(by_pass, Site == "Canyon" & Species == "BRN") # BRN
k_canyon_BRN_est = removal(k_canyon_BRN$Catch)
summary(k_canyon_BRN_est)
confint(k_canyon_BRN_est) # lower se than leslie

k_canyon_RBT = subset(by_pass, Site == "Canyon" & Species == "RBT") # RBT
k_canyon_RBT_est = removal(k_canyon_RBT$Catch)
summary(k_canyon_RBT_est)
confint(k_canyon_RBT_est) # higher se than leslie

k_canyon_BRK = subset(by_pass, Site == "Canyon" & Species == "BRK") #BRK
k_canyon_BRK_est = removal(k_canyon_BRK$Catch)
summary(k_canyon_BRK_est)
confint(k_canyon_BRK_est) # lower se than leslie

# Cherry
k_cherry_BRN = subset(by_pass, Site == "Cherry") # BRN only
k_cherry_BRN_est  = removal(k_cherry_BRN$Catch)
summary(k_cherry_BRN_est)
confint(k_cherry_BRN_est) # lower se than leslie

# Cox
k_cox_BRN = subset(by_pass, Site == "Cox Slough" & Species == "BRN") # BRN only
k_cox_BRN_est = removal(k_cox_BRN$Catch)
summary(k_cox_BRN_est)
confint(k_cox_BRN_est) 

# Deep
k_deep_BRN = subset(by_pass, Site == "Deep" & Species == "BRN") # BRN
k_deep_BRN_est = removal(k_deep_BRN$Catch)
summary(k_deep_BRN_est) # leslie does not work, this is best
confint(k_deep_BRN_est) 

k_deep_RBT = subset(by_pass, Site == "Deep" & Species == "RBT") # RBT
k_deep_RBT_est = removal(k_deep_RBT$Catch)
summary(k_deep_RBT_est) # higher se than leslie
confint(k_deep_RBT_est) 

k_deep_BRK = subset(by_pass, Site == "Deep" & Species == "BRK") # BRK
k_deep_BRK_est = removal(k_deep_BRK$Catch)
summary(k_deep_BRK_est) # leslie does not work, this is best
confint(k_deep_BRK_est) 

# Grasshopper
k_grasshopper_BRN = subset(by_pass, Site == "Grasshopper" & Species == "BRN") # BRN
k_grasshopper_BRN_est = removal(k_grasshopper_BRN$Catch)
summary(k_grasshopper_BRN_est) # leslie method better due to very good regression fit
confint(k_grasshopper_BRN_est) 

k_grasshopper_BRK = subset(by_pass, Site == "Grasshopper" & Species == "BRK") # BRK
k_grasshopper_BRK_est = removal(k_grasshopper_BRK$Catch)
summary(k_grasshopper_BRK_est) # higher se than leslie
confint(k_grasshopper_BRK_est) 

# Mill
k_mill= subset(by_pass, Site == "Mill" & Species == "BRN") # BRN only
k_mill_BRN_est = removal(k_mill$Catch)
summary(k_mill_BRN_est) # higher se than leslie
confint(k_mill_BRN_est) 

# Moose
k_moose_BRN = subset(by_pass, Site == "Moose" & Species == "BRN") # BRN
k_moose_BRN_est = removal(k_moose_BRN$Catch)
summary(k_moose_BRN_est)
confint(k_moose_BRN_est) # lower se than leslie

k_moose_RBT = subset(by_pass, Site == "Moose" & Species == "RBT") # RBT
k_moose_RBT_est = removal(k_moose_RBT$Catch)
summary(k_moose_RBT_est)
confint(k_moose_RBT_est) # lower se than leslie

# Silver Springs
k_silver_BRN = subset(by_pass, Site == "Silver Springs" & Species == "BRN") # BRN
k_silver_BRN_est = removal(k_silver_BRN$Catch)
summary(k_silver_BRN_est)
confint(k_silver_BRN_est) # higher se than leslie

k_silver_RBT = subset(by_pass, Site == "Silver Springs" & Species == "RBT") # RBT
k_silver_RBT_est = removal(k_silver_RBT$Catch)
summary(k_silver_RBT_est)
confint(k_silver_RBT_est) # higher se than leslie

# Stone
k_stone= subset(by_pass, Site == "Stone" & Species == "BRN") # BRN only
k_stone_BRN_est = removal(k_stone$Catch)
summary(k_stone_BRN_est)
confint(k_stone_BRN_est) # leslie does not work, this is better

# Trapper
k_trapper= subset(by_pass, Site == "Trapper" & Species == "BRN") # BRN only 
k_trapper_BRN_est = removal(k_trapper$Catch)
summary(k_trapper_BRN_est) 
confint(k_trapper_BRN_est) # higher se than leslie

# Wisconsin
k_wisconsin= subset(by_pass, Site == "Wisconsin" & Species == "BRN") # BRN only 
k_wisconsin_BRN_est = removal(k_wisconsin$Catch)
summary(k_wisconsin_BRN_est)
confint(k_wisconsin_BRN_est) # higher se than leslie


#### Summary ####

# K pass has lower standard error in half of the estimates compared to leslie method
# highly linear depletions favor the leslie method, such as in Grasshopper BRN
# it appears that the K pass results in smaller confidence intervals however, need to plot and find out 

leslie_estimate = c(summary(leslie_model_alder)[1], summary(leslie_model_blacktail_BRN)[1],
                    summary(leslie_model_camp_BRN)[1], summary(leslie_model_camp_BRK)[1], summary(leslie_model_canyon_BRN)[1],
                    summary(leslie_model_canyon_RBT)[1],summary(leslie_model_canyon_BRK)[1],summary(leslie_model_cherry_BRN)[1],
                    summary(leslie_model_cox_BRN)[1],summary(leslie_model_deep_RBT)[1],
                    summary(leslie_model_grasshopper_BRN)[1],summary(leslie_model_grasshopper_BRK)[1],
                    summary(leslie_model_mill_BRN)[1], summary(leslie_model_moose_BRN)[1],summary(leslie_model_moose_RBT)[1],
                    summary(leslie_model_silver_BRN)[1],summary(leslie_model_silver_RBT)[1],summary(leslie_model_trapper_BRN)[1],
                    summary(leslie_model_wisconsin_BRN)[1])

leslie_se = c(summary(leslie_model_alder)[3], summary(leslie_model_blacktail_BRN)[3],
              summary(leslie_model_camp_BRN)[3], summary(leslie_model_camp_BRK)[3], summary(leslie_model_canyon_BRN)[3],
              summary(leslie_model_canyon_RBT)[3],summary(leslie_model_canyon_BRK)[3],summary(leslie_model_cherry_BRN)[3],
              summary(leslie_model_cox_BRN)[3],summary(leslie_model_deep_RBT)[3],
              summary(leslie_model_grasshopper_BRN)[3],summary(leslie_model_grasshopper_BRK)[3],
              summary(leslie_model_mill_BRN)[3], summary(leslie_model_moose_BRN)[3],summary(leslie_model_moose_RBT)[3],
              summary(leslie_model_silver_BRN)[3],summary(leslie_model_silver_RBT)[3],summary(leslie_model_trapper_BRN)[3],
              summary(leslie_model_wisconsin_BRN)[3])

leslie_lci = c(confint(leslie_model_alder)[1], confint(leslie_model_blacktail_BRN)[1],
                    confint(leslie_model_camp_BRN)[1], confint(leslie_model_camp_BRK)[1], confint(leslie_model_canyon_BRN)[1],
                    confint(leslie_model_canyon_RBT)[1],confint(leslie_model_canyon_BRK)[1],confint(leslie_model_cherry_BRN)[1],
                    confint(leslie_model_cox_BRN)[1],confint(leslie_model_deep_RBT)[1],
                    confint(leslie_model_grasshopper_BRN)[1],confint(leslie_model_grasshopper_BRK)[1],
                    confint(leslie_model_mill_BRN)[1], confint(leslie_model_moose_BRN)[1],confint(leslie_model_moose_RBT)[1],
                    confint(leslie_model_silver_BRN)[1],confint(leslie_model_silver_RBT)[1],confint(leslie_model_trapper_BRN)[1],
                    confint(leslie_model_wisconsin_BRN)[1])

leslie_uci = c(confint(leslie_model_alder)[3], confint(leslie_model_blacktail_BRN)[3],
              confint(leslie_model_camp_BRN)[3], confint(leslie_model_camp_BRK)[3], confint(leslie_model_canyon_BRN)[3],
              confint(leslie_model_canyon_RBT)[3],confint(leslie_model_canyon_BRK)[3],confint(leslie_model_cherry_BRN)[3],
              confint(leslie_model_cox_BRN)[3],confint(leslie_model_deep_RBT)[3],
              confint(leslie_model_grasshopper_BRN)[3],confint(leslie_model_grasshopper_BRK)[3],
              confint(leslie_model_mill_BRN)[3], confint(leslie_model_moose_BRN)[3],confint(leslie_model_moose_RBT)[3],
              confint(leslie_model_silver_BRN)[3],confint(leslie_model_silver_RBT)[3],confint(leslie_model_trapper_BRN)[3],
              confint(leslie_model_wisconsin_BRN)[3])

k_estimate = c(summary(k_alder_BRN_est)[1], summary(k_blacktail_BRN_est)[1],
               summary(k_camp_BRN_est)[1], summary(k_camp_BRK_est)[1], summary(k_canyon_BRN_est)[1],
               summary(k_canyon_RBT_est)[1], summary(k_canyon_BRK_est)[1], summary(k_cherry_BRN_est)[1],
               summary(k_cox_BRN_est)[1], summary(k_deep_BRN_est)[1],
               summary(k_grasshopper_BRN_est)[1], summary(k_grasshopper_BRK_est)[1],
               summary(k_mill_BRN_est)[1], summary(k_moose_BRN_est)[1], summary(k_moose_RBT_est)[1],
               summary(k_silver_BRN_est)[1], summary(k_silver_RBT_est)[1], summary(k_trapper_BRN_est)[1],
               summary(k_wisconsin_BRN_est)[1])

k_lci = c(confint(k_alder_BRN_est)[1], confint(k_blacktail_BRN_est)[1],
               confint(k_camp_BRN_est)[1], confint(k_camp_BRK_est)[1], confint(k_canyon_BRN_est)[1],
               confint(k_canyon_RBT_est)[1], confint(k_canyon_BRK_est)[1], confint(k_cherry_BRN_est)[1],
               confint(k_cox_BRN_est)[1], confint(k_deep_BRN_est)[1],
               confint(k_grasshopper_BRN_est)[1], confint(k_grasshopper_BRK_est)[1],
               confint(k_mill_BRN_est)[1], confint(k_moose_BRN_est)[1], confint(k_moose_RBT_est)[1],
               confint(k_silver_BRN_est)[1], confint(k_silver_RBT_est)[1], confint(k_trapper_BRN_est)[1],
               confint(k_wisconsin_BRN_est)[1])

k_uci = c(confint(k_alder_BRN_est)[3], confint(k_blacktail_BRN_est)[3],
          confint(k_camp_BRN_est)[3], confint(k_camp_BRK_est)[3], confint(k_canyon_BRN_est)[3],
          confint(k_canyon_RBT_est)[3], confint(k_canyon_BRK_est)[3], confint(k_cherry_BRN_est)[3],
          confint(k_cox_BRN_est)[3], confint(k_deep_BRN_est)[3],
          confint(k_grasshopper_BRN_est)[3], confint(k_grasshopper_BRK_est)[3],
          confint(k_mill_BRN_est)[3], confint(k_moose_BRN_est)[3], confint(k_moose_RBT_est)[3],
          confint(k_silver_BRN_est)[3], confint(k_silver_RBT_est)[3], confint(k_trapper_BRN_est)[3],
          confint(k_wisconsin_BRN_est)[3])


Site = c("Alder", "Blacktail", "Camp", "Camp", "Canyon","Canyon","Canyon","Cherry","Cox","Deep","Grasshopper","Grasshopper",
         "Mill", "Moose", "Moose", "Silver", "Silver", "Trapper","Wisconsin")

Species = c("BRN", "BRN","BRN","BRK","BRN","RBT","BRK","BRN","BRN","RBT","BRN","BRK","BRN","BRN","RBT",
            "BRN","RBT","BRN","BRN")

df = data.frame(cbind(Site, Species, leslie_estimate, leslie_se, leslie_lci, leslie_uci, k_estimate, k_lci, k_uci))
df$leslie_estimate = as.numeric(df$leslie_estimate)
df$leslie_se = as.numeric(df$leslie_se)
df$leslie_lci = as.numeric(df$leslie_lci)
df$leslie_uci = as.numeric(df$leslie_uci)
df$k_estimate = as.numeric(df$k_estimate)
df$k_lci = as.numeric(df$k_lci)
df$k_uci = as.numeric(df$k_uci)
view(df)

lci = (df$leslie_uci - df$leslie_lci)
kci = (df$k_uci - df$k_lci)
view(lci - kci) # the K-pass method produces smaller confidence intervals at every site besides grasshopper

ggplot(df, aes(x=Site, y = leslie_estimate, color = Species)) +
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbar(aes(ymin = leslie_lci, ymax = leslie_uci, width = 0.5, linetype = "95% CI",), 
                position = position_dodge(0.5), size = 0.8) +
  scale_y_continuous(breaks = seq(-205, 2200, by = 200)) +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05)) +
  scale_linetype_manual(values = c("95% CI" = "solid"), name = "Error Bars") +
  labs (y = "Leslie Abundance")

ggplot(df, aes(x=Site, y = k_estimate, color = Species)) +
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbar(aes(ymin = k_lci, ymax = k_uci, width = 0.5, linetype = "95% CI",), 
                position = position_dodge(0.5), size = 0.8) +
  scale_y_continuous(breaks = seq(0, 2000, by = 200)) +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05)) +
  scale_linetype_manual(values = c("95% CI" = "solid"), name = "Error Bars") +
  labs (y = "K-Pass Abundance")


