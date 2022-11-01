#Multi-stressor field experiment
#GAM analysis
#A. Ostrowski  CJ Brown 
# 2022-11-01

library(tidyverse)
library(mgcv)
library(visreg)
library(patchwork) #for plot layout

theme_set(theme_classic())

#load csv - data without control
dat <- read.csv("../Data/Data_NoControl.csv")
str(dat)
dat <- data.frame(unclass(dat), stringsAsFactors = TRUE)
str(dat)
dat$Plot.ID <- as.factor(dat$Plot.ID)
dat$Block <- as.factor(dat$Block)
str(dat)

sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")


#
#shoot density analysis - Gaussian/normal distribution
#

m1 <- gam(Shoot_density ~ s(Week, by = Block, k = 3) +
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + 
            offset(Day0_density),
          data = dat)
gam.check(m1)

#seems somewhat skewed

#
#shoot density analysis - poisson distribution
#

#noting shoot density is a count, so I used a coujnt model
# (poisson). SEems to model the increase of variance with the
# mean better than the normal distribution did
# this is effectively a model of the proportional change: 
# density/Day0_density, because the poisson uses a log link
#(also note we therefore need to log the offset)
# whereas the gaussian above is modelling the difference
# in the number of shoots:
# density - Day0_density

m1 <- gam(Shoot_density ~ s(Week, by = Block, k = 3) +
            Treatment + 
            # The 'by' allows random time trends by blocks
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + #random intercepts for blocks
            offset(log(Day0_density)),
          #note for poisson we need to log the offset
          family = "poisson",
          data = dat)

gam.check(m1)
plot(m1)
summary(m1)

#Signif treatmetn effects, but also signif
# block time trends
#Note the EDF values. Values >1 imply non-linearity in time trend
# e.g. inphase has slight non-linearity
# and by nonlinear I mean nonlinaer on log scale - all linear trends
# are non-linaer in the sense that we use a log linke, therefore they are 
#exponential curve fits

#make nice plots 
#Treatment effects
predout <- visreg(m1, xvar = "Week", 
       by = "Treatment",
       scale = "response",
       plot = FALSE,
       alpha = 0.05) # can plot directly with TRUE here
# or have more control if you save results and plot with ggplot
#returns results with 95% CIs (not SEs). 

g1 <- ggplot(predout$fit) +
  aes(x = Week, y = visregFit, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
              alpha = 0.2, color = NA)


#Block effects
predout_block <- visreg(m1, xvar = "Week", 
                  by = "Block",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05) # can plot directly with TRUE here


g2 <- ggplot(predout_block$fit) +
  aes(x = Week, y = visregFit, color = Block,
      group = Block, fill = Block)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
              alpha = 0.2, color = NA)

#magic of patchwork!

g1 + g2 + plot_annotation(tag_levels = 'A')

#for this I take there is: 
# all are losing seagrass 
# slower losses under static
# fasteer losses with outphase
# fastsest loses with inphase
# static signif different by ~week 3
#  inphase/outphase diverge by ~week 4
# Some signif block variation

# So you can try different things now and plot results
# e.g. see if block interaction matters?

m2 <- gam(Shoot_density ~ s(Week, k = 3) +
            Treatment + 
            # The 'by' allows random time trends by blocks
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + #random intercepts for blocks
            offset(log(Day0_density)),
          #note for poisson we need to log the offset
          family = "poisson",
          data = dat)
AIC(m1, m2) 
#block interaction seems to be very important!
#But you may want to plot results and see if m2 has
# a different interpretation for treatment
