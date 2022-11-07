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

#plots 
#Treatment effects
predout <- visreg(m2, xvar = "Week", 
                  by = "Treatment",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05) 

g1 <- ggplot(predout$fit) +
  aes(x = Week, y = visregFit, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
              alpha = 0.2, color = NA)

#Block effects
predout_block <- visreg(m2, xvar = "Week", 
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

g1 + g2 + plot_annotation(tag_levels = 'A')

#Similar model output from first model 
#same trends in treatment effects - all treatments reduced shoot density
#static slowest reduction
#out phase greater reduction, in phase greatest reduction
#effect of block still
#model 1 indicates there is variability of block effect over time,
#with greater reduction in shoots in blocks 2 and 4 starting at week 4



#
#AO additional response variables - Nov 11, 2022
#leaf surface area
m1 <- gam(Leaf_surface_area ~ s(Week, by = Block, k = 3) +
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + 
            offset(Day0_surface_area),
          data = dat)
gam.check(m1)

m2 <- gam(Leaf_surface_area ~ s(Week, by = Block, k = 3) +
            Treatment + 
            s(Block, bs = "re") +
            offset(Day0_surface_area),
          data = dat)

gam.check(m2)

m3 <- gam(Leaf_surface_area ~ s(Week, by = Block, k = 3) +
            Treatment + 
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + 
            offset(Day0_surface_area),
          data = dat)

gam.check(m3)
plot(m3)
summary(m3)

#plots
#Treatment effects
predout <- visreg(m3, xvar = "Week", 
                  by = "Treatment",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05)  

g1 <- ggplot(predout$fit) +
  aes(x = Week, y = visregFit, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
              alpha = 0.2, color = NA)
g1

#Block effects
predout_block <- visreg(m3, xvar = "Week", 
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
g2

g1 + g2 + plot_annotation(tag_levels = 'A')

#higher variability/larger CIs for LSA
#R-sq value is low for all three models
#significant block effect - greater reduction across weeks in block 2
#gradual reduction in static treatment
#greater reduction in in and out phase treatments, although overlap CIs with static treatment
#model summary indicates significant effect of both in and out phase treatments

#can we log transform LSA?
dat$lnarea <- log(dat$Leaf_surface_area)
dat$lnDay0area <- log(dat$Day0_surface_area)

m4 <- gam(lnarea ~ s(Week, by = Block, k = 3) + Treatment + 
            s(Week, by = Treatment, k = 3) + s(Block, bs = "re") +
            offset(lnDay0area), data = dat)

gam.check(m4)
plot(m4)
summary(m4)
#seems to be better fit..

#plots
#Treatment effects
predout <- visreg(m4, xvar = "Week", 
                  by = "Treatment",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05)  

g1 <- ggplot(predout$fit) +
  aes(x = Week, y = visregFit, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
              alpha = 0.2, color = NA)
g1

#Block effects
predout_block <- visreg(m4, xvar = "Week", 
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
g2

g1 + g2 + plot_annotation(tag_levels = 'A')
#significant effect of block 4 over time


#crustacean abundance
m1 <- gam(Total_crustacean ~ s(Week, by = Block, k = 3) +
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + 
            offset(Day0_crustacean),
          data = dat)
gam.check(m1)

m2 <- gam(Total_crustacean ~ s(Week, by = Block, k = 3) +
            Treatment + 
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + 
            offset(Day0_crustacean),
          data = dat)
gam.check(m2)

#poisson - this is count data..
m3 <- gam(Total_crustacean ~ s(Week, by = Block, k = 3) + Treatment + 
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
          family = "poisson",
          data = dat)
gam.check(m3)
plot(m3)
summary(m3)

m4 <- gam(Total_crustacean ~ s(Week, by = Block, k = 3) + Treatment + 
            s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
          family = "poisson",
          data = dat)

gam.check(m4)

#log transform data?
dat$lncrustacean <- log(dat$Total_crustacean + 0.01)
dat$lnDay0crustacean <- log(dat$Day0_crustacean + 0.01)

m5 <- gam(lncrustacean ~ s(Week, by = Block, k = 3) +
            Treatment + 
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + 
            offset(lnDay0crustacean),
          data = dat)
gam.check(m5)

#square root transform?
dat$SRcrustacean <- sqrt(dat$Total_crustacean)
dat$SRDay0crustacean <- sqrt(dat$Day0_crustacean)

m6 <- gam(SRcrustacean ~ s(Week, by = Block, k = 3) +
            Treatment + 
            s(Week, by = Treatment, k = 3) + 
            s(Block, bs = "re") + 
            offset(SRDay0crustacean),
          data = dat)
gam.check(m6)
plot(m6)
summary(m6)

#plots
#Treatment effects
predout <- visreg(m3, xvar = "Week", 
                  by = "Treatment",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05)  

g1 <- ggplot(predout$fit) +
  aes(x = Week, y = visregFit, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
              alpha = 0.2, color = NA)
g1

#Block effects
predout_block <- visreg(m3, xvar = "Week", 
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
g2

g1 + g2 + plot_annotation(tag_levels = 'A')

#model 3 seems best fit
#low R-sq values across all models
#significant effect of block over time - high variability across blocks and time
#effect across all treatments
#static treatment reduces abundance
#in and out phase treatments indicate greater reduction than static
