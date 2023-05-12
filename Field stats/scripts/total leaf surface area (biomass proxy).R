#Multi-stressor field experiment
#Supplementary material
#Figure S2. Change in total leaf surface area per plot


library(mgcv)
library(readr)
library(visreg)
library(ggplot2)
library(cowplot)
library(tidyverse)

#read csv
dat <- read.csv("data/Data_Control.csv")
str(dat)
dat <- data.frame(unclass(dat), stringsAsFactors = TRUE)
str(dat)
dat$Block <- as.factor(dat$Block)
str(dat)

sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

#identify individual and interactive effects of treatment and week, random block effect
m1 <- gam(Leaf_surface_area ~ s(Week, by = Block, k = 4) +
            Treatment + s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") + offset(Day0_surface_area),
          data = dat)

gam.check(m1)
plot(m1)
summary(m1)

dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))

#make model predictions
#treatment effects
predout <- visreg(m1, xvar = "Week", 
                  by = "Treatment",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05)

#plot predictions
ggplot(predout$fit) +
  aes(x = Week, y = visregFit - Day0_surface_area, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_ribbon(aes(ymin = visregLwr - Day0_surface_area, ymax = visregUpr - Day0_surface_area),
              alpha = 0.2, color = NA) +
  theme_cowplot()+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6"))+
  xlab("Week")+
  ylab("Model predictions \n change in total leaf surface area (cm^2 per 0.36 m^2 plot)")
