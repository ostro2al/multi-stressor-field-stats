#Estimate treatment effects using CBs custom function
# for simulating CIs from GAMs
# CJ Brown
# 2023-03-25

library(tidyverse)
library(mgcv)
library(visreg)
library(cowplot)

#NB If you want to/can be bothered recreating your GAM graphs with 
# credible intervals (won't make much or any difference) then
# use code below but modify the creatoin of the prediction 
# dataframe with expand.grid() to get all combinations of 
# weeks and treatments. Then simulate predictions, calculate quantiles and finally 
# then plot as normal with ggplot2


#read csv
dat <- read.csv("data/Data_Control.csv")

#custon function
source("scripts/functions-sim-predictions.R")

#
# Parameters
#

k_val <- 4
dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))
dat$ln_Day0_density <- log(dat$Day0_density)
#log here rather than in the formula, makes using prediction
# below more intuitive

#
# Shoot density effect sizes 
#

m3.1 <- gam(Shoot_density ~ s(Week, by = Block, k = k_val) +
              Treatment + s(Week, by = Treatment, k = k_val) + 
              s(Block, bs = "re") + offset(ln_Day0_density),
            family = "poisson", data = dat)

#
# Create new dataframe for predictions 
#

newd <- with(dat, data.frame(Week = 6,
                   Block = 1,
                   Treatment = unique(Treatment),
                   ln_Day0_density = mean(ln_Day0_density))
)

newd$pred <- predict(m3.1, newdata = newd, type = "response")
newd


#row IDs for treatments (so we don't mix them up)
icontrol <- which(newd$Treatment == "Control")
istatic <- which(newd$Treatment == "StaticStatic")


#Predictions with credible intervals
# A list of formulas to evaluate sample by sample
# the function then takes quantiles over these outcomes
# use names if you want the output dataframes and 
# columns for CIs to be named. 

forms <- c(effects = ~exp(x), 
           control_mult = ~ exp(x - x[icontrol]),
          static_mult ~ exp(x - x[istatic]),
          prob_control = ~ x<x[icontrol],
           prob_static = ~ x<x[istatic])

gout <- simulate_gam_CIs(m3.1, 
                 newdata = newd,
                 forms = forms, 
                 random_var = "Block",
                 offset = mean(dat$ln_Day0_density),
                 probs = c(0.025, 0.5, 0.975),
                 nsims = 1000)

