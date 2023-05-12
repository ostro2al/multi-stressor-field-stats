#Multi-stressor field experiment
#Supplementary material
#Figure S1. Seagrass shoot survival rate (%)

library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)
library(tidyverse)

#read csv
dat <- read.csv("data/MortalityGrowth_Control.csv")
str(dat)
dat <- data.frame(unclass(dat), stringsAsFactors = TRUE)
str(dat)
dat$Block <- as.factor(dat$Block)
dat$Treatment <- as.factor(dat$Treatment)
str(dat)

dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

#effects of treatment on survival rate, random block effect
m1 <- lmer(Survival_rate ~ Treatment + (1|Block), dat = dat)

plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)

#formatted figure with CI
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-pred.se,
                conf.high=pred0+pred.se))
}

#Make df with all combos of observed treatments
dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))

newdata <- dat %>% 
  select(Treatment) %>%
  distinct()
newdata$mean <- predict(m1, re.form=NA,
                        newdata = newdata)

CIs <- easyPredCI(m1, newdata = newdata)
newdata$lwr <- CIs[,1]
newdata$upr <- CIs[,2]

#plot predictions
newdata %>%
  ggplot() +
  aes(x = Treatment, y = mean, colour = Treatment) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  theme_cowplot()+
  ylab("Model predictions \n seagrass shoot survival rate (%)") +
  xlab("Treatment") 
