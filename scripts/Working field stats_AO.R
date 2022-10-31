#Multi-stressor field experiment
#lmer analysis
#A. Ostrowski October 31, 2022


library(dplyr)
library(rstatix)
library(ggplot2)
library(emmeans)
library(tidyverse)
library(lme4)
library(lmerTest)

#load csv - data without control
dat <- read.csv("Data_NoControl.csv")
str(dat)
dat <- data.frame(unclass(dat), stringsAsFactors = TRUE)
str(dat)
dat$Plot.ID <- as.factor(dat$Plot.ID)
dat$Week <- as.factor(dat$Week)
dat$Block <- as.factor(dat$Block)
str(dat)

#shoot density analysis
#shoot density model and plot
sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

m1 <- lmer(Shoot_density ~ Week * Treatment + (1|Plot.ID) 
           + offset(Day0_density), dat = dat)

plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
anova(m1)

#transform data
#sqrt
dat$SRdensity <- sqrt(dat$Shoot_density)
hist(dat$SRdensity)
dat$SRDay0_density <- sqrt(dat$Day0_density)
hist(dat$SRDay0_density)

m2 <- lmer(SRdensity ~ Week * Treatment + (1|Plot.ID) 
           + offset(SRDay0_density), dat = dat)

shapiro.test(resid(m2))
hist(resid(m2))
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
summary(m2)
anova(m2)

#log transform - not a good fit
dat$lndensity <- log(dat$Shoot_density)
hist(dat$lndensity)
dat$lnDay0_density <- log(dat$Day0_density)
hist(dat$lnDay0_density)

m3 <- lmer(lndensity ~ Week * Treatment + (1|Plot.ID) 
           + offset(lnDay0_density), dat = dat)

shapiro.test(resid(m3))
hist(resid(m3))
plot(m3)
qqnorm(resid(m3))
qqline(resid(m3))
summary(m3)
anova(m3)

#investigating other lmers 
#add "block" as random effect nested within PlotID random effect
#shoot density model
sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

m4 <- lmer(SRdensity ~ Week * Treatment + (1|Block) + (1|Block:Plot.ID) 
           + offset(SRDay0_density), dat = dat)

shapiro.test(resid(m4))
hist(resid(m4))
plot(m4)
qqnorm(resid(m4))
qqline(resid(m4))
summary(m4)
anova(m4)

#plot
#formatted figure with SE
easyPredSE <- function(model,newdata=NULL,alpha=0.05) {
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

#Make df with all combos of observed treatments (m1)
dat$Treatment <- factor(dat$Treatment, levels = c("StaticStatic", "InPhase", "OutPhase"))

newdata <- dat %>% 
  select(Treatment, Week, Day0_density) %>%
  distinct()
newdata$Day0_density <- mean(dat$Day0_density)
newdata$mean <- predict(m1, re.form=NA,
                        newdata = newdata)

SEs <- easyPredSE(m1, newdata = newdata)
newdata$lwr <- SEs[,1]
newdata$upr <- SEs[,2]

newdata %>%
  ggplot() +
  aes(x = Week, y = mean - Day0_density, colour = Treatment, group = Treatment) +
  geom_errorbar(aes(ymin = lwr- Day0_density, ymax = upr- Day0_density), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  #geom_smooth(method = "lm", alpha = .15, aes(fill = Treatment))+
  theme_bw()+
  theme(text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 0, vjust=0.6, color = "black"),
        axis.text.y = element_text(color = "black"))+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  #geom_hline(yintercept = 0, colour = "grey42", size = 0.5, linetype = "longdash")+
  ylab("Model predictions \n change in shoot density relative to day 0") +
  xlab("Week") +
  scale_colour_discrete("Treatment") +
  theme(legend.position = c(0.88, 0.84))


#
##
#surface area model and plot
sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

m5 <- lmer(Leaf_surface_area ~ Week * Treatment + (1|Plot.ID) 
           + offset(Day0_surface_area), dat = dat)

shapiro.test(resid(m5))
hist(resid(m5))
plot(m5)
qqnorm(resid(m5))
qqline(resid(m5))
summary(m5)
anova(m5)

#sqrt transform
dat$SRarea <- sqrt(dat$Leaf_surface_area)
hist(dat$SRarea)
dat$SRDay0_area <- sqrt(dat$Day0_surface_area)
hist(dat$SRDay0_area)

m6 <- lmer(SRarea ~ Week * Treatment + (1|Plot.ID) 
           + offset(SRDay0_area), dat = dat)

shapiro.test(resid(m6))
hist(resid(m6))
plot(m6)
qqnorm(resid(m6))
qqline(resid(m6))
summary(m6)
anova(m6)

#block and plot ID random effects
m7 <- lmer(SRarea ~ Week * Treatment + (1|Block) + (1|Block:Plot.ID) 
           + offset(SRDay0_area), dat = dat)

plot(m7)
qqnorm(resid(m7))
qqline(resid(m7))
summary(m7)
anova(m7)

#plot
#formatted figure with SE
easyPredSE <- function(model,newdata=NULL,alpha=0.05) {
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

#Make df with all combos of observed treatments (m5)
dat$Treatment <- factor(dat$Treatment, levels = c("StaticStatic", "InPhase", "OutPhase"))

newdata <- dat %>% 
  select(Treatment, Week, Day0_surface_area) %>%
  distinct()
newdata$Day0_surface_area <- mean(dat$Day0_surface_area)
newdata$mean <- predict(m5, re.form=NA,
                        newdata = newdata)

SEs <- easyPredSE(m5, newdata = newdata)
newdata$lwr <- SEs[,1]
newdata$upr <- SEs[,2]

newdata %>%
  ggplot() +
  aes(x = Week, y = mean - Day0_surface_area, colour = Treatment, group = Treatment) +
  geom_errorbar(aes(ymin = lwr- Day0_surface_area, ymax = upr- Day0_surface_area), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  #geom_smooth(method = "lm", alpha = .15, aes(fill = Treatment))+
  theme_bw()+
  theme(text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 0, vjust=0.6, color = "black"),
        axis.text.y = element_text(color = "black"))+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  ylab("Model predictions \n change in leaf surface area relative to day 0") +
  xlab("Week") +
  scale_colour_discrete("Treatment") +
  theme(legend.position = c(0.88, 0.84))


#
##
#crustacean abundance model and plot
sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

hist(dat$Total_crustacean)

m8 <- lmer(Total_crustacean ~ Week * Treatment + (1|Plot.ID) 
           + offset(Day0_crustacean), dat = dat)

plot(m8)
qqnorm(resid(m8))
qqline(resid(m8))
summary(m8)
anova(m8)

#sqrt transform
dat$SRinvert <- sqrt(dat$Total_crustacean)
hist(dat$SRinvert)
dat$SRDay0_invert <- sqrt(dat$Day0_crustacean)
hist(dat$SRDay0_invert)

m9 <- lmer(SRinvert ~ Week * Treatment + (1|Plot.ID) 
           + offset(SRDay0_invert), dat = dat)

shapiro.test(resid(m9))
hist(resid(m9))
plot(m9)
qqnorm(resid(m9))
qqline(resid(m9))
summary(m9)
anova(m9)

#block and plot ID random effects
m10 <- lmer(SRinvert ~ Week * Treatment + (1|Block) + (1|Block:Plot.ID) 
           + offset(SRDay0_invert), dat = dat)

shapiro.test(resid(m10))
hist(resid(m10))
plot(m10)
qqnorm(resid(m10))
qqline(resid(m10))
summary(m10)
anova(m10)

#plot
#formatted figure with SE
easyPredSE <- function(model,newdata=NULL,alpha=0.05) {
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

#Make df with all combos of observed treatments (m8)
dat$Treatment <- factor(dat$Treatment, levels = c("StaticStatic", "InPhase", "OutPhase"))

newdata <- dat %>% 
  select(Treatment, Week, Day0_crustacean) %>%
  distinct()
newdata$Day0_crustacean <- mean(dat$Day0_crustacean)
newdata$mean <- predict(m8, re.form=NA,
                        newdata = newdata)

SEs <- easyPredSE(m8, newdata = newdata)
newdata$lwr <- SEs[,1]
newdata$upr <- SEs[,2]

newdata %>%
  ggplot() +
  aes(x = Week, y = mean - Day0_crustacean, colour = Treatment, group = Treatment) +
  geom_errorbar(aes(ymin = lwr- Day0_crustacean, ymax = upr- Day0_crustacean), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  #geom_smooth(method = "lm", alpha = .15, aes(fill = Treatment))+
  theme_bw()+
  theme(text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 0, vjust=0.6, color = "black"),
        axis.text.y = element_text(color = "black"))+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  ylab("Model predictions \n change in crustacean abundance relative to day 0") +
  xlab("Week") +
  scale_colour_discrete("Treatment") +
  theme(legend.position = c(0.88, 0.84))


#
##
#survival rate
#data without control
dat2 <- read.csv("MortalityGrowth_Data.csv")
str(dat2)
dat2 <- data.frame(unclass(dat2), stringsAsFactors = TRUE)
str(dat2)
dat2$Block <- as.factor(dat2$Block)
str(dat2)

dat2$Treatment <- relevel(factor(dat2$Treatment), ref = "StaticStatic")

m1 <- lmer(Survival_rate ~ Treatment + (1|Block), dat = dat2)

plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
anova(m1)

#formatted figure with SE
easyPredSE <- function(model,newdata=NULL,alpha=0.05) {
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

#Make df with all combos of observed treatments (m2)
dat2$Treatment <- factor(dat2$Treatment, levels = c("StaticStatic", "InPhase", "OutPhase"))

newdata <- dat2 %>% 
  select(Treatment) %>%
  distinct()
newdata$mean <- predict(m1, re.form=NA,
                        newdata = newdata)

SEs <- easyPredSE(m1, newdata = newdata)
newdata$lwr <- SEs[,1]
newdata$upr <- SEs[,2]

newdata %>%
  ggplot() +
  aes(x = Treatment, y = mean, colour = Treatment) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  theme_bw()+
  theme(text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 0, vjust=0.6, color = "black"),
        axis.text.y = element_text(color = "black"))+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  ylab("Model predictions \n seagrass survival rate (%)") +
  xlab("Treatment") 

#use binomial model instead of lmer
dat2$Treatment <- relevel(factor(dat2$Treatment), ref = "StaticStatic")

m2 <- glmer(cbind(Shoot_fin, Shoot_init) ~ Treatment + (1|Block), dat = dat2, family = "binomial")

summary(m2)
anova(m2)
plot(predict(m2))