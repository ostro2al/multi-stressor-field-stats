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

#####--------------------------- ##########
### CB Psem models
# - mediation of crustaceans (no direct treatment effect)

m1_pSEM_randomList_CB = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat),
  Avg_LSA %~~% Shoot_density
)

summary(m1_pSEM_randomList_CB, .progressBar = F)
#suggests week effect needed for crusts

m2_pSEM_randomList_CB = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat),
  Avg_LSA %~~% Shoot_density
)

summary(m2_pSEM_randomList_CB, .progressBar = F)
#fits ok. 
plot(m2_pSEM_randomList_CB)

#compare to no week effect on crusts, but with stressor app effect

m3_pSEM_randomList_CB = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Stressor_app + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat),
  Avg_LSA %~~% Shoot_density
)

summary(m3_pSEM_randomList_CB, .progressBar = F)
#doesn't fit

#end cb models
#####--------------------------- ##########

### psem 
#Latest code revisions (25.01.2023) - AO
#revised psem with random block effect

#load csv with data stacked to include day 0 data
dat <- read.csv("DataStacked_Control.csv")
str(dat)
dat$Stressor_app <- as.numeric(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
str(dat)

#try different pathway hypotheses
#should we also look at effects of LSA on density instead of correlaton between the variables??
#week * treatment interaction, mediation of crustaceans
m1 = psem(
  lme(Shoot_density ~ Stressor_app * Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app * Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat),
  Avg_LSA %~~% Shoot_density
)

summary(m1, .progressBar = F)

#no interaction, mediation of crustaceans (no direct treatment effect)
m2 = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat),
  Avg_LSA %~~% Shoot_density
)

summary(m2, .progressBar = F)
#doesn't fit, suggests week effect needed for crusts

#add effect of week on crustaceans
m3 = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat),
  Avg_LSA %~~% Shoot_density
)

summary(m3, .progressBar = F)
anova(m3)
plot(m3)
#check residuals 
plot(resid(lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat)))
plot(resid(lme(Avg_LSA ~ Stressor_app + Week, random = ~ 1 | Block, data = dat)))
plot(resid(lme(Crustacean_abundance ~ Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)))

#try week and stressor app effect on crustaceans
m4 = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat),
  Avg_LSA %~~% Shoot_density
)

summary(m4, .progressBar = F)
#same results as m3, no sig effect of treatment on crustaceans

#
#compare all models that fit: m1, m3, m4 
anova(m1, m3, m4)
#Chi-square Difference Test
#AIC    BIC Fisher.C  Fisher.C.Diff  DF.diff      P.value    
#1    39.106 98.461    1.106                               
#vs 2 33.106 83.089    1.106             0       0       1 
#vs 3 34.000 87.107    0.000         1.106       2  0.5752 

AIC(m1, m3)
#df    AIC
#x 18 37.106
#y 16 33.106

#m1 summary
#AIC      BIC
#37.106   93.337
#Fisher's C = 1.106 with P-value = 0.575 and on 2 degrees of freedom

#m3
#AIC      BIC
#33.106   83.089
#Fisher's C = 1.106 with P-value = 0.575 and on 2 degrees of freedom

#m4
#AIC      BIC
#34.000   87.107
#Fisher's C = 0 with P-value = 1 and on 0 degrees of freedom


#
##try running pSEM with GAMs 
str(dat)
dat$Stressor_app <- as.numeric(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
dat$Block <- as.numeric(dat$Block)
str(dat)
dat

sort(dat$Week)

m1_gam = psem(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Stressor_app + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Stressor_app + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
        Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), family = nb(), data = dat)
)
#issues getting the psem to run

nrow(dat)
length(dat[,1])

length(unique(dat[,1]))==length(dat[,1])

dups <- which(duplicated(dat))
length(dups)
unique(dat[dups])

#
#for reference
#individual GAMs already tested - for input into pSEM
m1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
            Treatment + 
            s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") +
            offset(log(Day0_density)),
          family = "poisson",
          data = dat)

m2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
            Treatment + 
            s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") + 
            offset(Day0_Avg_LSA),
          data = dat)

m3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
            s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
          family = nb(),
          data = dat)