#Multi-stressor field experiment
#GAM analysis - updated models
#A. Ostrowski  CJ Brown 
#March 3, 2023

library(tidyverse)
library(mgcv)
library(visreg)
library(cowplot)

#read csv
dat <- read.csv("../Data/Data_Control.csv")
str(dat)
dat <- data.frame(unclass(dat), stringsAsFactors = TRUE)
str(dat)
dat$Plot.ID <- as.factor(dat$Plot.ID)
dat$Block <- as.factor(dat$Block)
str(dat)

sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

#shoot density
#identify individual and interactive effects of treatment and week, random block effect
m1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
            Treatment + s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") + offset(log(Day0_density)),
          family = "poisson", data = dat)

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
  aes(x = Week, y = visregFit - Day0_density, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr - Day0_density, ymax = visregUpr - Day0_density),
              alpha = 0.2, color = NA)+
  theme_cowplot()+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6"))+
  xlab("Week")+
  ylab("Model predictions \n change in shoot density per plot relative to day 0")


#leaf surface area
#identify individual and interactive effects of treatment and week, random block effect
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

m2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
            Treatment + s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") + offset(Day0_Avg_LSA),
          data = dat)

gam.check(m2)
plot(m2)
summary(m2)

dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))

#make model predictions
#treatment effects
predout <- visreg(m2, xvar = "Week", 
                  by = "Treatment",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05) 

#plot predictions
ggplot(predout$fit) +
  aes(x = Week, y = visregFit - Day0_Avg_LSA, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr - Day0_Avg_LSA, ymax = visregUpr - Day0_Avg_LSA),
              alpha = 0.2, color = NA) +
  theme_cowplot()+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6"))+
  xlab("Week")+
  ylab("Model predictions \n change in leaf surface area (cm^2) relative to day 0")


#crustacean abundance
#identify individual and interactive effects of treatment and week, random block effect
library(countreg)

rootogram <- function(yobs, ypred, brks = NA, hanging = FALSE, ...){
  
  if (is.na(brks)){
    brks <- seq(0-0.5, ceiling(max(c(ypred, yobs))) +0.5, by = 1)
  } 
  xwd <- diff(brks)[1]/2.1
  xy <- hist(ypred, breaks = brks, plot = F)
  xy2 <- hist(yobs, breaks = brks, plot = F)
  ylwr <-  sqrt(xy$counts) - sqrt(xy2$counts)
  ylim <- c(min(ylwr), sqrt(max(xy$counts)))
  xlim <- c(min(brks), max(brks))
  plot(xy$mids, sqrt(xy$counts), lty = 3, col = 'red', type = 'b', ylim = ylim,
       xlim = xlim, xlab ='Count', ylab = 'sqrt(Frequency)', ...)
  
  if (hanging){
    rect(xy$mids - xwd, sqrt(xy$counts), xy$mids + xwd,ylwr, col = 'grey')
  }else{
    rect(xy$mids - xwd, sqrt(xy2$counts), xy$mids + xwd,0, col = 'grey')
  }
  
  lines(xy$mids, sqrt(xy$counts), lty = 3, type = 'b', pch = 16, col = 'red')
  abline(h = 0)
}

dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

m3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
            s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
          family = nb(),
          data = dat)

dat$fit <- predict(m3, type = "response")
with(dat, rootogram(Crustacean_abundance, fit))

gam.check(m3)
plot(m3)
summary(m3)

dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))

#make model predictions
#treatment effects
predout <- visreg(m3, xvar = "Week", 
                  by = "Treatment",
                  scale = "response",
                  plot = FALSE,
                  alpha = 0.05)

#plot predictions
ggplot(predout$fit) +
  aes(x = Week, y = visregFit - Day0_crustacean, color = Treatment,
      group = Treatment, fill = Treatment)+ 
  geom_line() + 
  geom_ribbon(aes(ymin = visregLwr - Day0_crustacean, ymax = visregUpr - Day0_crustacean),
              alpha = 0.2, color = NA)+
  theme_cowplot()+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6"))+
  xlab("Week")+
  ylab("Model predictions \n change in crustacean abundance relative to day 0")
