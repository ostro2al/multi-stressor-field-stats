#Estimate treatment effects using CBs custom function
# for simulating CIs from GAMs
# CJ Brown
# 2023-03-25

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

k_val <- 4
dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))
dat$ln_Day0_density <- log(dat$Day0_density)

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



#row IDs for treatments
icontrol <- which(newd$Treatment == "Control")
istatic <- which(newd$Treatment == "StaticStatic")


#Predictions with credible intervals
# A list of formulas to evaluate sample by sample
# the function then takes quantiles over these outcomes
# use names if you want the output dataframes and 
# columns for CIs to be named. 

forms <- c(effects = ~exp(x), #estimate median shoot density
           control_mult = ~ exp(x - x[icontrol]),
           #Estimate each treatment as multiple of control
          static_mult= ~ exp(x - x[istatic]),
          #Estimate each treatment as multiple of static
          prob_control = ~ x<x[icontrol],
          #Estimate probability treatment is less than control
           prob_static = ~ x<x[istatic])
        #Estimate probability treatment is less than static

#Specify sumamry functions, ie for effect sizes use quantiles to 
# get CIs,
# whereas for inequalities we just want to sum and divide by nsims
# to get the probability 
m3.1_functions <- c("quantile","quantile",
                    "quantile", "sum",
                    "sum")

gout <- simulate_gam_CIs(m3.1, 
                 newdata = newd,
                 forms = forms, 
                 random_var = "Block",
                 offset = mean(dat$ln_Day0_density),
                 probs = c(0.025, 0.5, 0.975),
                 nsims = 1000,
                 func_to_apply = m3.1_functions)

shoot_density_effects <- gout[[1]] %>% left_join(gout[[2]]) %>%
  left_join(gout[[3]]) %>%
  left_join(gout[[4]]) %>%
  left_join(gout[[5]])

shoot_density_effects


#
# LSA
#

m3.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = k_val) +
              Treatment + s(Week, by = Treatment, k = k_val) + 
              s(Block, bs = "re") + offset(Day0_Avg_LSA),
            data = dat)

#
# Create new dataframe for predictions 
#

newd2 <- with(dat, data.frame(Week = 6,
                             Block = 1,
                             Treatment = unique(Treatment),
                             Day0_Avg_LSA = mean(Day0_Avg_LSA))
)

#row IDs for treatments
icontrol2 <- which(newd2$Treatment == "Control")
istatic2 <- which(newd2$Treatment == "StaticStatic")


#Predictions with credible intervals
# A list of formulas to evaluate sample by sample
# the function then takes quantiles over these outcomes
# use names if you want the output dataframes and 
# columns for CIs to be named. 

forms <- c(effects = ~x, #estimate median LSA
           control_mult = ~ x / x[icontrol2],
           #Estimate each treatment as multiple of control
           static_mult= ~ x / x[istatic2],
           #Estimate each treatment as multiple of static
           prob_control = ~ x<x[icontrol2],
           #Estimate probability treatment is less than control
           prob_static = ~ x<x[istatic2])
#Estimate probability treatment is less than static

#Specify summary functions, ie for effect sizes use quantiles to 
# get CIs,
# whereas for inequalities we just want to sum and divide by nsims
# to get the probability 
m3.2_functions <- c("quantile","quantile",
                    "quantile", "sum",
                    "sum")

gout <- simulate_gam_CIs(m3.2, 
                         newdata = newd2,
                         forms = forms, 
                         random_var = "Block",
                         offset = mean(dat$Day0_Avg_LSA),
                         probs = c(0.025, 0.5, 0.975),
                         nsims = 1000,
                         func_to_apply = m3.2_functions)

LSA_effects <- gout[[1]] %>% left_join(gout[[2]]) %>%
  left_join(gout[[3]]) %>%
  left_join(gout[[4]]) %>%
  left_join(gout[[5]])

LSA_effects




# 
# Crustaceans
# 


#model 3
m3 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Block -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Block -> Avg_LSA
                Week -> Crustacean_abundance
                Treatment -> Crustacean_abundance
                Block -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )

#
# Total effect 
#
adjustmentSets(m3, "Treatment", "Crustacean_abundance", type="canonical")

#Condition on week and block to get total traetment effect (via any path)
# so fit model without shoots and LSA

dat$ln_Day0_crustacean <- log(dat$Day0_crustacean+ 0.01)
m3_crust_total <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
                        s(Week, by = Treatment, k = 4) + 
                        s(Block, bs = "re") + offset(ln_Day0_crustacean),
                      family = nb(), data = dat)

#now use simulate_gam_CIs() to get CIs for the total effect of Treatment

newd3 <- with(dat, data.frame(Week = 6,
                              Block = 1,
                              Treatment = unique(Treatment),
                              ln_Day0_crustacean = mean(ln_Day0_crustacean))
)

#row IDs for treatments (so we don't mix them up)
icontrol3 <- which(newd3$Treatment == "Control")
istatic3 <- which(newd3$Treatment == "StaticStatic")

#Predictions with credible intervals
forms <- c(effects = ~exp(x), #estimate median LSA
           control_mult = ~ exp(x - x[icontrol3]),
           #Estimate each treatment as multiple of control
           static_mult= ~ exp(x - x[istatic3]),
           #Estimate each treatment as multiple of static
           prob_control = ~ x<x[icontrol3],
           #Estimate probability treatment is less than control
           prob_static = ~ x<x[istatic3])
#Estimate probability treatment is less than static

#Specify sumamry functions, ie for effect sizes use quantiles to 
# get CIs,
# whereas for inequalities we just want to sum and divide by nsims
# to get the probability 
m3_crust_total_functions <- c("quantile","quantile",
                    "quantile", "sum",
                    "sum")

gout <- simulate_gam_CIs(m3_crust_total, 
                         newdata = newd3,
                         forms = forms, 
                         random_var = "Block",
                         offset = mean(dat$ln_Day0_crustacean),
                         probs = c(0.025, 0.5, 0.975),
                         nsims = 1000,
                         func_to_apply = m3_crust_total_functions)

crust_total_effects <- gout[[1]] %>% left_join(gout[[2]]) %>%
  left_join(gout[[3]]) %>%
  left_join(gout[[4]]) %>%
  left_join(gout[[5]])

crust_total_effects



#
# Direct effects of LSA and treatment
#

adjustmentSets(m3, "Avg_LSA", "Crustacean_abundance", type="canonical")


dat$ln_Day0_crustacean <- log(dat$Day0_crustacean+ 0.01)
m3_crust_directs <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
      s(Week, by = Treatment, k = 4) + s(Shoot_density, k = 4) + Avg_LSA +
      s(Block, bs = "re") + offset(ln_Day0_crustacean),
    family = nb(), data = dat)

#use this model to get direct effects of LSA on crusties (conditional on treatment)
# use same model to get direct effect of treatment on crusties (conditional on LSA)

newd4 <- with(dat, expand.grid(Week = 6,
                             Block = 1,
                             Treatment = unique(Treatment),
                             Shoot_density = mean(Shoot_density), 
                             Avg_LSA = c(mean(Avg_LSA) - sd(Avg_LSA),
                                         mean(Avg_LSA),
                                         mean(Avg_LSA) + sd(Avg_LSA)),
                             ln_Day0_crustacean = mean(ln_Day0_crustacean))
)


#set your condition, here all results will 
#be multiples of static treatment and mean LSA

istatic4 <- which(newd4$Treatment == "StaticStatic" & newd4$Avg_LSA == mean(dat$Avg_LSA))

icontrol4 <- which(newd4$Treatment == "Control" & newd4$Avg_LSA == mean(dat$Avg_LSA))

iinphase4 <- which(newd4$Treatment == "InPhase" & newd4$Avg_LSA == mean(dat$Avg_LSA))


#Now formula for results relative to static at mean LSA
forms2 <- c(static_LSA_mult = ~ exp(x - x[istatic4]),
            iinphase_LSA_mult = ~ exp(x - x[iinphase4]),
           #Estimate each treatment as multiple of static
           c(control_LSA_mult = ~ exp(x - x[icontrol4]),
           prob_LSA_static = ~ x<x[istatic4]))
#Estimate probability treatment is less than static



forms2_functions <- c("quantile", "quantile", "quantile", "sum")


gout <- simulate_gam_CIs(m3_crust_directs, 
                         newdata = newd4,
                         forms = forms2, 
                         random_var = "Block",
                         offset = mean(dat$ln_Day0_crustacean),
                         probs = c(0.025, 0.5, 0.975),
                         nsims = 1000,
                         func_to_apply = forms2_functions)

crust_direct_effects <- gout[[1]] %>% left_join(gout[[2]]) 


crust_direct_effects



#Here are conditional direct effects of treatments
# relative to the static treatment
# at the mean LSA
crust_direct_effects[crust_direct_effects$Avg_LSA == mean(dat$Avg_LSA),]

#Here are conditional direct effect of LSA for +- one SD around
#mean LSA and 
# at the static treatment
crust_direct_effects[crust_direct_effects$Treatment == "StaticStatic",]
crust_direct_effects[crust_direct_effects$Treatment == "InPhase",]


crust_direct_effects

