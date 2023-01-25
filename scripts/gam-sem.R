# GAM solution 
#under development
# idea is to get indep claims, then test them
# manually, then calculate fisher C, also manually

library(piecewiseSEM)
library(mgcv)

dat <- read.csv("data/Data_Control.csv")
dat$Stressor_app <- factor(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
dat$Block <- factor(dat$Block)

bs_gam = basisSet(list(
  lm(Shoot_density ~ Week*Block +
        Stressor_app + Week*Stressor_app + 
        Block, data = dat),
  lm(Avg_LSA ~ Week*Block  +
        Stressor_app + Week*Stressor_app + 
        Block, data = dat),
  lm(Crustacean_abundance ~ Avg_LSA + Shoot_density + Week*Block +
        Stressor_app + Week*Stressor_app + 
        Block, data = dat)
))
