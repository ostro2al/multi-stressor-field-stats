#Multi-stressor field experiment
#Mechanistic Analysis - using SEMss
#A. Ostrowski November 25, 2022


library(nlme) 
library(lavaan)
library(lme4)
library(piecewiseSEM)
library(dplyr)

dat <- read.csv("Data_Control.csv")
str(dat)
dat$Stressor_app <- as.numeric(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
str(dat)

dat <- dat %>% na.omit()
dat


### psem 
#with random block effect
m1_pSEM_randomList = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)
)

summary(m1_pSEM_randomList, .progressBar = F)


m2_pSEM_randomList = psem(
  lme(Shoot_density ~ Stressor_app + Week + Avg_LSA, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)
)

summary(m2_pSEM_randomList, .progressBar = F)


#without random block effect
m1_pSEM = psem(
  lm(Shoot_density ~ Stressor_app + Week, data = dat),
  lm(Avg_LSA ~ Stressor_app + Week, data = dat),
  lm(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, data = dat)
)

summary(m1_pSEM, .progressBar = F)


m2_pSEM = psem(
  lm(Shoot_density ~ Stressor_app + Week + Avg_LSA, data = dat),
  lm(Avg_LSA ~ Stressor_app  + Week, data = dat),
  lm(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, data = dat)
)

summary(m2_pSEM, .progressBar = F)

plot(m2_pSEM)

dSep(m2_pSEM)
coefs(m2_pSEM, standardize = "none")
residuals(m2_pSEM)
plot(residuals(m2_pSEM))


plot(m2_pSEM, return = FALSE,
     node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
     edge_attrs = data.frame(style = "solid", color = "black"),
     ns_dashed = T,
     alpha = 0.05,
     show = "std",
     digits = 3,
     add_edge_label_spaces = TRUE)

#model 2 has lowest AIC and BIC in psem - lowest without random block effect

#transform crustacean abundance
dat$lncrust <- log(dat$Crustacean_abundance + 0.001)
dat$SRcrust <- sqrt(dat$Crustacean_abundance)

m3_pSEM = psem(
  lm(Shoot_density ~ Stressor_app + Week + Avg_LSA, data = dat),
  lm(Avg_LSA ~ Stressor_app  + Week, data = dat),
  lm(lncrust ~ Stressor_app + Week + Shoot_density + Avg_LSA, data = dat)
)

summary(m3_pSEM, .progressBar = F)

plot(m3_pSEM)

dSep(m3_pSEM)
fisherC(m3_pSEM)
coefs(m3_pSEM, standardize = "none")
residuals(m3_pSEM)


#
#Latest code revisions (13.01.2023) - AO
#revised psem with random block effect and interaction between treatment and week

#load csv with data stacked to include day 0 data
dat <- read.csv("DataStacked_Control.csv")
str(dat)
dat$Stressor_app <- as.numeric(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
str(dat)

dat <- dat %>% na.omit()
dat

#new method to construct models
##break down individual regressions (of hypothesised paths of interest) to explore individual
#and interactive effects, then use AIC to select best regression to input into SEM

#shoot density
m1 <- lmer(Shoot_density ~ Stressor_app + Week + Avg_LSA + (1 | Block), data = dat)
m2 <- lmer(Shoot_density ~ Stressor_app * Week + Avg_LSA + (1 | Block), data = dat)
m3 <- lmer(Shoot_density ~ Stressor_app * Week * Avg_LSA + (1 | Block), data = dat)
m4 <- lmer(Shoot_density ~ Stressor_app + Week * Avg_LSA + (1 | Block), data = dat)

#m2 lowest AIC, m3 highest df
AIC(m1, m2, m3, m4)
#   df      AIC
#m1  6 1654.328
#m2  7 1646.950
#m3 10 1648.833
#m4  7 1658.982

plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))


#average leaf surface area
m1 <- lmer(Avg_LSA ~ Stressor_app * Week + (1 | Block), data = dat)
m2 <- lmer(Avg_LSA ~ Stressor_app + Week + (1 | Block), data = dat)

#m2 lower AIC, m1 higher df
AIC(m1, m2)
#   df      AIC
#m1  6 1262.265
#m2  5 1261.839

plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))


#crustacean abundance
m1 <- lmer(Crustacean_abundance ~ Week + Shoot_density + Avg_LSA + (1 | Block), data = dat)
m2 <- lmer(Crustacean_abundance ~ Week * Shoot_density + Avg_LSA + (1 | Block), data = dat)

#m1 lowest AIC, m2 higher df
AIC(m1, m2)
#df      AIC
#m1  6 1098.990
#m2  7 1108.985

plot(m2)
qqnorm(resid(m2))
qqline(resid(m2)) #not good fit

#log transform crustacean abundance
dat$lncrust <- log(dat$Crustacean_abundance + 0.001)
#square root transform crustacean abundance - seems best fit
dat$SRcrust <- sqrt(dat$Crustacean_abundance)

m1 <- lmer(SRcrust ~ Week + Shoot_density + Avg_LSA + (1 | Block), data = dat)
m2 <- lmer(SRcrust ~ Week * Shoot_density + Avg_LSA + (1 | Block), data = dat)

#m1 has lowest AIC, m2 highest df
AIC(m1, m2)
#   df      AIC
#m1  6 541.2747
#m2  7 554.6181

plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))

m1 <- lmer(lncrust ~ Week + Shoot_density + Avg_LSA + (1 | Block), data = dat)
m2 <- lmer(lncrust ~ Week * Shoot_density + Avg_LSA + (1 | Block), data = dat)

#m1 lowest AIC, m2 highest df
AIC(m1, m2)
#   df      AIC
#m1  6 890.2207
#m2  7 902.2558


##
#now try SEM with regressions selected based on lowest AIC - second is best (using SQRT transform crustacean variable)
mnew_pSEM_randomList = psem(
  lme(Shoot_density ~ Stressor_app * Week + Avg_LSA, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)
)

summary(mnew_pSEM_randomList, .progressBar = F)

plot(mnew_pSEM_randomList)


mnew2_pSEM_randomList = psem(
  lme(Shoot_density ~ Stressor_app * Week + Avg_LSA, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(SRcrust ~ Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)
)

summary(mnew2_pSEM_randomList, .progressBar = F)

plot(mnew2_pSEM_randomList)


mnew3_pSEM_randomList = psem(
  lme(Shoot_density ~ Stressor_app * Week + Avg_LSA, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(lncrust ~ Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)
)

summary(mnew3_pSEM_randomList, .progressBar = F)

plot(mnew3_pSEM_randomList)