
library(tidyverse)
library(mgcv)
library(visreg)
library(cowplot)

#v2 is handy to play around with teh code
# v3 uses a function

#NB If you want to/can be bothered recreating your GAM graphs with 
# credible intervals (won't make much or any difference) then
# use code below but modify the creatoin of the prediction 
# dataframe with expand.grid() to get all combinations of 
# weeks and treatments. Then simulate predictions, calculate quantiles and finally 
# then plot as normal with ggplot2


#read csv
dat <- read.csv("data/Data_Control.csv")
str(dat)

#
# Parameters
#

nsims <- 1000
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

#Predicting
Xp <- predict(m3.1, newdata = newd, type="lpmatrix") 

#Set block effects to zero
# this means our predictions remove
# block variation, we just focus on
# overall trends
# (they are 'conditional predictions')
# identify columns of the Xp matrix with 'block' 
# then set those columns to zero to ignore block effects

icol_block <- grepl("Block", colnames(Xp))
Xp[,icol_block] <- 0


#prediction of mean, same as newd$pred
#note the Xp matrix doesn't include the offset, so 
# we need to add that back on
exp(Xp %*% coef(m3.1) + newd$ln_Day0_density[1])

#
# Simulate random draws from the posterior to get intervals
#



br <- rmvn(nsims, coef(m3.1), m3.1$Vp) ## 1000 replicate param. vectors

#Setup some vectors/matrices to store random draws of
# predictions 
treatment_effects <- matrix(0,ncol = nsims, nrow = nrow(newd),nsims)
#treatment multiples of control
control_mult <- matrix(0, nrow = nrow(newd),
                       ncol = nsims)
#treatment multiples of static
control_mult <- static_mult <- matrix(0, nrow = nrow(newd),
                       ncol = nsims)
 
#row IDs for treatments (so we don't mix them up)
icontrol <- which(newd$Treatment == "Control")
istatic <- which(newd$Treatment == "StaticStatic")

for (i in 1:1000){ 
  pr <- Xp %*% br[i,] + newd$ln_Day0_density[1] ## replicate predictions
  #adding offset back in. This only matters for mean estimates
  #not important for the differences (sience the offset will cancel out)
  
  treatment_effects[,i] <- exp(pr)
  
  #Differences between treatments as multiples
  # use exp(x) - exp(y) if you want differences
  #in number shoots
  # or exp(x-y) if you want multiples
  
  control_mult[,i] <- exp(pr - pr[icontrol])
  static_mult[,i]  <- exp(pr - pr[istatic])
    
}

#
# Compute CIs
#

#treatment medians and CIs
treat_dat <- t(apply(treatment_effects, 1, quantile, probs = c(0.025, 0.5, 0.975))) %>%
  data.frame() %>%
  setNames(paste0('treatment_', names(.)))

#control multiples
control_mult_dat <- t(apply(control_mult, 1, quantile, probs = c(0.025, 0.5, 0.975))) %>%
  data.frame() %>%
  setNames(paste0('control_mult_', names(.)))


#static multiples
static_mult_dat <- t(apply(static_mult, 1, quantile, probs = c(0.025, 0.5, 0.975))) %>%
  data.frame() %>%
  setNames(paste0('static_mult_', names(.)))

#join them all
newd2 <- cbind(newd, treat_dat, control_mult_dat, 
               static_mult_dat)
newd2

#if you want probability that it is less then control then just sum
# up how many of the 1000 draws are less than a multiple of 1
newd2$control_mult_prob <- rowSums(control_mult<1)/1000
newd2$static_mult_prob <- rowSums(static_mult<1)/1000
  
  
  
  
  