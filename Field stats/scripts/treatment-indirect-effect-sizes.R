#Estimate indirect treatment effects 
# for simulating CIs from GAMs
# CJ Brown
# 2023-03-27

#Based on Imai et al. 2010 "A General Approach to Causal Mediation Analysis" 
# see also: https://cran.r-project.org/web/packages/mediation/vignettes/mediation.pdf
#
# Note one assumption is that there are no other causally related mediation effects
# therefore to estimate indirect effects we need to assume no shoot density effect
# (so its not quite model 3)
#

library(tidyverse)
library(mgcv)
library(visreg)
library(cowplot)
library(dagitty)

#read csv
dat <- read.csv("data/Data_Control.csv")

#custom function
source("scripts/functions-sim-predictions.R")

#
# Parameters
#
nsims <- 1000
k_val <- 4
dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))
dat$ln_Day0_density <- log(dat$Day0_density)


# 
# Crustacean mediation effect
# 
# Calculated using method of Imai et al. 2010

dat$ln_Day0_crustacean <- log(dat$Day0_crustacean+ 0.01)
m3_crust_directs <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
                          s(Week, by = Treatment, k = 4) + 
                          Avg_LSA +
                          s(Block, bs = "re") + offset(ln_Day0_crustacean),
                        family = nb(), data = dat)
m3.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
              Treatment + s(Week, by = Treatment, k = 4) + 
              s(Block, bs = "re") + offset(Day0_Avg_LSA),
            data = dat)


#
# Data to predict to 
#

newd <- with(dat, data.frame(Week = 6,
                             Block = 1,
                             Treatment = unique(Treatment),
                             ln_Day0_density = mean(ln_Day0_density),
                             ln_Day0_crustacean = mean(ln_Day0_crustacean),
                             Day0_Avg_LSA = mean(Day0_Avg_LSA))
)


# ----------------------
# Calculate mean indirect effect to check results 
# ----------------------

#Use lp matrix method to predict desnity for the average block
lp_LSA <- predict(m3.2, newdata = newd, type="lpmatrix")
icol_block <- grepl("Block", colnames(lp_LSA))
lp_LSA[,icol_block] <- 0
LSA_mean <- lp_LSA %*% coef(m3.2) + newd$Day0_Avg_LSA
LSA_mean
# Should be the same if we turn off the line that sets
#block effects to zero: 
predict(m3.2, newdata = newd)

newd2 <- expand.grid(Avg_LSA = as.numeric(LSA_mean),
                     Treatment = newd$Treatment)
newd2$Treatment_applied_to_LSA <- rep(newd$Treatment, 4)
newd2$Week <- 6
newd2$Block <- 1
newd2$ln_Day0_crustacean = mean(dat$ln_Day0_crustacean)

newd2$Avg_LSA <- as.numeric(LSA_mean)
Xp_crusties <- predict(m3_crust_directs, 
                       newdata = newd2, type="lpmatrix") 
icol_block <- grepl("Block", colnames(Xp_crusties))
Xp_crusties[,icol_block] <- 0
lncrustie_mean <- Xp_crusties %*% coef(m3_crust_directs) + newd2$ln_Day0_crustacean
newd2$lncrustie_mean <- as.numeric(lncrustie_mean)

# predict(m3_crust_directs, 
        # newdata = newd_temp, type = "response")
#Now apply equation 1 of Imai et al. 2010 

newd3 <- newd2 %>%
  group_by(Treatment) %>%
  mutate(CME = exp(lncrustie_mean - lncrustie_mean[Treatment_applied_to_LSA == "Control"])) %>%
  ungroup() %>%
  group_by(Treatment_applied_to_LSA) %>%
  mutate(DE = exp(lncrustie_mean - lncrustie_mean[Treatment == "Control"])) %>%
  ungroup()
data.frame(newd3)

#ACME: average across different direct effects
newd3 %>%
  group_by(Treatment_applied_to_LSA) %>%
  summarize(mean(CME))

#Average direct effect ADE
newd3 %>%
  group_by(Treatment) %>%
  summarize(mean(DE))

#Note these values are biased, so we should use values below


# ------------------
#Predict LSA
# ------------------
Xp_LSA <- predict(m3.2, newdata = newd, type="lpmatrix") 

#Set block effects to zero
icol_block <- grepl("Block", colnames(Xp_LSA))
Xp_LSA[,icol_block] <- 0

# Simulate random draws from the posterior to get intervals

br <- rmvn(nsims, coef(m3.2), m3.2$Vp) ## 1000 replicate param. vectors

#Setup some vectors/matrices to store random draws of
# predictions 
#needs to be a 3D array
treatment_LSA <- treatment_effects <- direct_control_mult <- direct_static_mult <- 
  indirect_control_mult <- indirect_static_mult <- 
  array(rep(0, nrow(newd)* nsims* nrow(newd)), c(nrow(newd), nrow(newd), nsims))
#Rows are treatment effects on crusties,
# columns are treatment effects on LSA,
# pages are simulations 

#row IDs for treatments (so we don't mix them up)
icontrol <- which(newd$Treatment == "Control")
istatic <- which(newd$Treatment == "StaticStatic")

for (i in 1:1000){ 
  # replicate prediction
  # for LSA
  pr_LSA <- Xp_LSA %*% br[i,] + newd$Day0_Avg_LSA[1] 
  
  for (iLSA in 1:nrow(newd)){
    newd_temp <- newd
    #for this one we need to calculate direct effect
    # under each treatment, conditional on LSA fixed 
    # at its value for every treatment
    # Loop over treatments to recalculate crusties direct effect
    # for each LSA indirect effect
    #
    newd_temp <- newd
    newd_temp$Avg_LSA <- as.numeric(pr_LSA)[iLSA]
    Xp_crusties <- predict(m3_crust_directs, 
                           newdata = newd_temp, type="lpmatrix") 
    icol_block2 <- grepl("Block", colnames(Xp_crusties))
    Xp_crusties[,icol_block2] <- 0
    br_crusties <- rmvn(1, coef(m3_crust_directs), m3_crust_directs$Vp)
    pr_crusties <- Xp_crusties %*% br_crusties[1,] + newd_temp$ln_Day0_density[1]
    
    #Differences between treatments as multiples
    treatment_effects[, iLSA,i] <- pr_crusties
    treatment_LSA[, iLSA,i] <- pr_LSA
  }
  
  #Would be more efficient to do this post-loop,
  # but easier to conceptualise this way
  
  #For indirect effect we need to do
  # outcome under each treatment (rows) given mediator at each treatment minus
  # outcome under each treatment given mediator at control (cols)
  #then take average across rows
  indirect_control_mult[,,i] <- treatment_effects[,,i] - 
    matrix(rep(treatment_effects[,icontrol,i], nrow(newd)), nrow(newd), byrow = FALSE)
  indirect_static_mult[,,i]  <-  treatment_effects[,,i] - 
    matrix(rep(treatment_effects[,istatic,i], nrow(newd)), nrow(newd), byrow = FALSE)
  
  #For direct effects we need to do: 
  # outcome under each treatment (rows) given mediator at each treatment minus
  # outcome under control given mediator at each treatment (cols)
  #then take average across pages

  direct_control_mult[,,i] <- treatment_effects[,,i] - 
    matrix(rep(treatment_effects[icontrol,,i], nrow(newd)), nrow(newd), byrow = TRUE)
  direct_static_mult[,,i]  <- treatment_effects[,,i] - 
    matrix(rep(treatment_effects[istatic,,i], nrow(newd)), nrow(newd), byrow = TRUE)
  
}


#
# Calculate direct and indirect effects
#

#Rows are direct treatment effects on crusties,
# columns are indirect treatment effects via LSA,
# pages are simulations 

#Indirect effects (ACME - average causal mediation effect)
# summarize across columns (because we are averaging over different direct 
# treatment effects to get average mediation effect)
# The exp() turns differences on a log scale to multiples. 
#
y = indirect_control_mult %>%
  apply(., c(2,3), mean) %>%
  exp() %>%
  apply(., 1, 
      quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t() %>%
  data.frame() %>%
  setNames(paste0('indirect_', names(.))) %>%
  cbind(newd)
y
y[,c(6,2)]

#Indirect effects (ADE - average direct effect)
# summarize across rows (because we are averaging over different direct 
# treatment effects to get average mediation effect)
x = direct_control_mult %>%
  apply(., c(1,3), mean) %>% 
  exp() %>%
  apply(., 1, 
      quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t() %>%
  data.frame() %>%
  setNames(paste0('direct_', names(.))) %>%
  cbind(newd)
x
x[,c(6,2)]


