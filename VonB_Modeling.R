
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


# loading in FSA croaker data
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




