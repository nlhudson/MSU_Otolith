
#### Von Bertalanffy Growth Modeling 

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
library(FSAdata)
library(nlstools)

#### FSA Example ####
## loading in FSA croaker data 
data(Croaker2)
crm <- subset(Croaker2,sex=="M")
view(crm)

svTypical <- vbStarts(tl~age,data=crm) # determining starting values for Linf, K, and t0
unlist(svTypical) 

vbTypical <- tl~Linf*(1-exp(-K*(age-t0))) # defining VB equation
fitTypical <- nls(vbTypical,data=crm,start=svTypical) # fitting the model to the data

overview(fitTypical) # overview of model fit to croaker data

crm$predicted_tl <- predict(fitTypical, newdata = crm) # extracting predicted model fit

ggplot(crm, aes(x = age, y = tl)) +
  geom_point(color = "blue", size = 2) + # Scatterplot of the original data
  geom_line(aes(y = predicted_tl), color = "red", size = 1) + # Fitted curve using extracted values
  labs(x = "Age", y = "Total Length (mm)", title = "Von Bertalanffy Model Fit") +
  theme_minimal()

bootTypical <- nlsBoot(fitTypical,niter=1000) # bootstrapping to obtain better parameter estimates and CIs
confint(bootTypical,plot=TRUE)

## Plotting Von B curve with 95% CIs based on bootstrapped parameter estimates
# 1. set age limits
age_seq <- 1:10

# 2. Use bootstrapped results to compute confidence intervals
boot_preds <- apply(bootTypical$coefboot, 1, function(params) {
  Linf <- params["Linf"]
  K <- params["K"]
  t0 <- params["t0"]
  Linf * (1 - exp(-K * (age_seq - t0)))
})

# 3. Compute mean and 95% confidence intervals
mean_preds <- rowMeans(boot_preds)
ci_lower <- apply(boot_preds, 1, quantile, probs = 0.025)
ci_upper <- apply(boot_preds, 1, quantile, probs = 0.975)

# 4. Create a data frame for plotting
plot_data <- data.frame(
  age = age_seq,
  mean_preds = mean_preds,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

# 5. Plot size-at-age data, fitted curve, and confidence intervals
ggplot(crm, aes(x = age, y = tl)) +
  geom_point(color = "blue", size = 2) +                           # Scatterplot of observed data
  geom_line(data = plot_data, aes(x = age, y = mean_preds),        # Fitted VB curve
            color = "red", size = 1) +
  geom_ribbon(data = plot_data, aes(x = age, y = mean_preds, ymin = ci_lower, ymax = ci_upper), 
              alpha = 0.2, fill = "black") +                   # Shaded 95% confidence intervals
  labs(x = "Age", y = "Total Length (mm)", 
       title = "Von Bertalanffy Model Fit with 95% Confidence Intervals") +
  theme_minimal()


#### Gannon Data ####

agedat = read_xlsx("Otolith_Ages_Clean.xlsx")
view(agedat)
agedat_bh = agedat[agedat$Basin %in% "BH",]
agedat_bv = agedat[agedat$Basin %in% "BV",]

startval_bh <- vbStarts(Length~Age,data=agedat_bh) # determining starting values for Linf, K, and t0
unlist(startval_bh) 

startval_bv <- vbStarts(Length~Age,data=agedat_bv) # determining starting values for Linf, K, and t0
unlist(startval_bv) 

vbTypical <- Length~Linf*(1-exp(-K*(Age-t0))) # defining VB equation

fitTypical_bh <- nls(vbTypical,data=agedat_bh,start=startval_bh) # fitting the model to the data
fitTypical_bv <- nls(vbTypical,data=agedat_bv,start=startval_bv) # fitting the model to the data

overview(fitTypical_bh)
overview(fitTypical_bv)

agedat_bh$predicted_tl <- predict(fitTypical_bh, newdata = agedat_bh) # extracting predicted model fit
agedat_bv$predicted_tl <- predict(fitTypical_bv, newdata = agedat_bv) # extracting predicted model fit


ggplot(agedat_bh, aes(x = Age, y = Length)) +
  geom_point(color = "blue", size = 2) + # Scatterplot of the original data
  geom_line(aes(y = predicted_tl), color = "red", size = 1) + # Fitted curve using extracted values
  labs(x = "Age", y = "Total Length (mm)", title = "Von Bertalanffy Model Fit") +
  theme_minimal()

ggplot(agedat_bv, aes(x = Age, y = Length)) +
  geom_point(color = "blue", size = 2) + # Scatterplot of the original data
  geom_line(aes(y = predicted_tl), color = "red", size = 1) + # Fitted curve using extracted values
  labs(x = "Age", y = "Total Length (mm)", title = "Von Bertalanffy Model Fit") +
  theme_minimal()

bootTypical_bh <- nlsBoot(fitTypical_bh,niter=1000) # bootstrapping to obtain better parameter estimates and CIs
confint(bootTypical_bh,plot=TRUE)

bootTypical_bv <- nlsBoot(fitTypical_bv,niter=1000) # bootstrapping to obtain better parameter estimates and CIs
confint(bootTypical_bv,plot=TRUE)

## Plotting Von B curve with 95% CIs based on bootstrapped parameter estimates
# 1. set age limits
age_seq <- 1:12

# 2. Use bootstrapped results to compute confidence intervals 
boot_preds_bh <- apply(bootTypical_bh$coefboot, 1, function(params) {
  Linf <- params["Linf"]
  K <- params["K"]
  t0 <- params["t0"]
  Linf * (1 - exp(-K * (age_seq - t0)))
})

boot_preds_bv <- apply(bootTypical_bv$coefboot, 1, function(params) {
  Linf <- params["Linf"]
  K <- params["K"]
  t0 <- params["t0"]
  Linf * (1 - exp(-K * (age_seq - t0)))
})

# 3. Compute mean and 95% confidence intervals
mean_preds_bh <- rowMeans(boot_preds_bh)
ci_lower_bh <- apply(boot_preds_bh, 1, quantile, probs = 0.025)
ci_upper_bh <- apply(boot_preds_bh, 1, quantile, probs = 0.975)

mean_preds_bv <- rowMeans(boot_preds_bv)
ci_lower_bv <- apply(boot_preds_bv, 1, quantile, probs = 0.025)
ci_upper_bv <- apply(boot_preds_bv, 1, quantile, probs = 0.975)

# 4. Create a data frame for plotting
plot_data_bh <- data.frame(
  age = age_seq,
  mean_preds_bh = mean_preds_bh,
  ci_lower_bh = ci_lower_bh,
  ci_upper_bh = ci_upper_bh
)



plot_data_bv <- data.frame(
  age = age_seq,
  mean_preds_bv = mean_preds_bv,
  ci_lower_bv = ci_lower_bv,
  ci_upper_bv = ci_upper_bv
)

#### Model Fit Plotting ####

# 5. Plot size-at-age data, fitted curve, and confidence intervals
ggplot(agedat_bh, aes(x = Age, y = Length)) +
  geom_point(color = "blue", size = 2) +                           # Scatterplot of observed data
  geom_line(data = plot_data_bh, aes(x = age, y = mean_preds_bh),        # Fitted VB curve
            color = "red", size = 1) +
  geom_ribbon(data = plot_data_bh, aes(x = age, y = mean_preds_bh, ymin = ci_lower_bh, ymax = ci_upper_bh), 
              alpha = 0.2, fill = "black") +                   # Shaded 95% confidence intervals
  labs(x = "Age", y = "Total Length (mm)", 
       title = "BH Von Bertalanffy Model Fit with 95% Confidence Intervals") +
  theme_minimal()

ggplot(agedat_bv, aes(x = Age, y = Length)) +
  geom_point(color = "blue", size = 2) +                           # Scatterplot of observed data
  geom_line(data = plot_data_bv, aes(x = age, y = mean_preds_bv),        # Fitted VB curve
            color = "red", size = 1) +
  geom_ribbon(data = plot_data_bv, aes(x = age, y = mean_preds_bv, ymin = ci_lower_bv, ymax = ci_upper_bv), 
              alpha = 0.2, fill = "black") +                   # Shaded 95% confidence intervals
  labs(x = "Age", y = "Total Length (mm)", 
       title = "BV Von Bertalanffy Model Fit with 95% Confidence Intervals") +
  theme_minimal()

von_b = 
ggplot(agedat, aes(x = Age, y = Length)) +
  geom_point(data = agedat, aes(x= Age, y = Length, color = Basin), size = 2) +
  geom_line(data = plot_data_bh, aes(x = age, y = mean_preds_bh),        # Fitted VB curve
            color = "blue", size = 1) +
  geom_ribbon(data = plot_data_bh, aes(x = age, y = mean_preds_bh, ymin = ci_lower_bh, ymax = ci_upper_bh), 
              alpha = 0.2, fill = "blue") + 
  geom_line(data = plot_data_bv, aes(x = age, y = mean_preds_bv),        # Fitted VB curve
            color = "red", size = 1) +
  geom_ribbon(data = plot_data_bv, aes(x = age, y = mean_preds_bv, ymin = ci_lower_bv, ymax = ci_upper_bv), 
              alpha = 0.2, fill = "red") +
  scale_color_manual(values = c("blue", "red")) + # Shaded 95% confidence intervals
  theme(text = element_text(size = 20)) +
  labs(x = "Age", y = "Fork Length (mm)")

ggsave(von_b, file = "Von_B_plot.jpg", dpi = 600)


#### Parameter extraction plotting ####

bh_params <- data.frame(
  Parameter = c("Linf", "K", "t0"),
  bootTypical_bh$estiboot,
  CI_lower = confint(bootTypical_bh)[, 1],  # Lower bound of CIs
  CI_upper = confint(bootTypical_bh)[, 2],  # Upper bound of CIs
  Site = "BH"
)
 
bv_params <- data.frame(
  Parameter = c("Linf", "K", "t0"),
  bootTypical_bv$estiboot,
  CI_lower = confint(bootTypical_bv)[, 1],  # Lower bound of CIs
  CI_upper = confint(bootTypical_bv)[, 2],  # Upper bound of CIs
  Site = "BV"
)

param_data <- rbind(bh_params, bv_params)

# plotting 
ggplot(param_data, aes(x = Parameter, y = Estimate, color = Site)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Parameter, scales = "free_y") +
  labs(x = "Parameter", y = "Value", title = "Parameter Estimates with 95% CIs") +
  theme(text=element_text(size=15)) +
  scale_color_manual(values = c("BH" = "blue", "BV" = "red"))



  
