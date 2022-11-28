#Multi-stressor field experiment
#Mechanistic Analysis - using SEMss
#A. Ostrowski November 25, 2022


library(ape)
library(caper)
library(nlme) 
library(lavaan)
library(lme4)
library(piecewiseSEM)

dat <- read.csv("Data_Control.csv")
str(dat)
dat$Stressor_app <- as.numeric(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
str(dat)

dat <- dat %>% na.omit()
dat

# Create interaction term
dat$treat_week_int <- dat$Stressor_app * dat$Week

View(dat)

# List structured equations for lavaan
m1 = '
  Shoot_density ~ Stressor_app + Week + Avg_LSA 
  Avg_LSA ~ Stressor_app  + Week + Shoot_density 
  Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA
'

# Fit vcov SEM
m1.sem = sem(m1, dat, estimator = "MLM")
varTable(m1.sem)

# Summary output with standardized coefficients
summary(m1.sem, standardize = TRUE)

inspect(m1.sem, "rsquare")
#Shoot_density         Avg_LSA        Crustacean_abundance 
#  0.343                0.258                0.083 
#does not give SEs, model not identified


m2 = '
  Shoot_density ~ Stressor_app + Week 
  Avg_LSA ~ Stressor_app  + Week  
  Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA
'

# Fit vcov SEM
m2.sem = sem(m2, dat, estimator = "MLM")
varTable(m2.sem)

# Summary output with standardized coefficients
summary(m2.sem, standardize = TRUE)

inspect(m2.sem, "rsquare")
#Shoot_density         Avg_LSA        Crustacean_abundance 
#  0.297                0.094                0.079 


m3 = '
  Shoot_density ~ Stressor_app + Week + Crustacean_abundance
  Avg_LSA ~ Stressor_app  + Week + Crustacean_abundance
  Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA
'

# Fit vcov SEM
m3.sem = sem(m3, dat, estimator = "MLM")
varTable(m3.sem)

# Summary output with standardized coefficients
summary(m3.sem, standardize = TRUE)

inspect(m3.sem, "rsquare")
#Shoot_density         Avg_LSA        Crustacean_abundance 
#  0.189               -0.877               -1.332 
#does not give SEs, model not identified


m4 = '
  Shoot_density ~ Stressor_app + Week + Avg_LSA
  Avg_LSA ~ Stressor_app  + Week
  Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA
'

# Fit vcov SEM
m4.sem = sem(m4, dat, estimator = "MLM")
varTable(m4.sem)

# Summary output with standardized coefficients
summary(m4.sem, standardize = TRUE)

inspect(m4.sem, "rsquare")
#Shoot_density         Avg_LSA        Crustacean_abundance 
#  0.428                0.094                0.083 

parameterestimates(m4.sem, standardize = TRUE) #CIs for parameters
fitted(m4.sem) #look at cov table
residuals(m4.sem) #look at residuals
fitmeasures(m4.sem) #fit indices


##
#try models 2, 4 with piecewiseSEM

### psem 

#with random block effect
m2_pSEM_randomList = psem(
  lme(Shoot_density ~ Stressor_app + Week, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)
)

summary(m2_pSEM_randomList, .progressBar = F)


m4_pSEM_randomList = psem(
  lme(Shoot_density ~ Stressor_app + Week + Avg_LSA, random = ~ 1 | Block, data = dat),
  lme(Avg_LSA ~ Stressor_app  + Week, random = ~ 1 | Block, data = dat),
  lme(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, random = ~ 1 | Block, data = dat)
)

summary(m4_pSEM_randomList, .progressBar = F)


#without random block effect
m2_pSEM = psem(
  lm(Shoot_density ~ Stressor_app + Week, data = dat),
  lm(Avg_LSA ~ Stressor_app + Week, data = dat),
  lm(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, data = dat)
)

summary(m2_pSEM, .progressBar = F)


m4_pSEM = psem(
  lm(Shoot_density ~ Stressor_app + Week + Avg_LSA, data = dat),
  lm(Avg_LSA ~ Stressor_app  + Week, data = dat),
  lm(Crustacean_abundance ~ Stressor_app + Week + Shoot_density + Avg_LSA, data = dat)
)

summary(m4_pSEM, .progressBar = F)

plot(m4_pSEM)

dSep(m4_pSEM)
fisherC(m4_pSEM)
coefs(m4_pSEM, standardize = "none")
residuals(m4_pSEM)
plot(residuals(m4_pSEM))
partialResid(Shoot_density ~ Avg_LSA, m4_pSEM)
partialResid(Crustacean_abundance ~ Avg_LSA, m4_pSEM)

plot(m4_pSEM, return = FALSE,
     node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
     edge_attrs = data.frame(style = "solid", color = "black"),
     ns_dashed = T,
     alpha = 0.05,
     show = "std",
     digits = 3,
     add_edge_label_spaces = TRUE)


#no differences in overall conclusions between lavaan and psem models
#no differences in overall conclusions between lme and lm models 
#model 4 has lowest AIC and BIC in psem - lowest without random block effect


#crustaceans are not normally distributed, non-linear responses
#can I transform variable?
#is this why we are not seeing significant correlation between seagrass variables and crustaceans?

dat$lncrust <- log(dat$Crustacean_abundance + 0.001)
dat$SRcrust <- sqrt(dat$Crustacean_abundance)

m5_pSEM = psem(
  lm(Shoot_density ~ Stressor_app + Week + Avg_LSA, data = dat),
  lm(Avg_LSA ~ Stressor_app  + Week, data = dat),
  lm(lncrust ~ Stressor_app + Week + Shoot_density + Avg_LSA, data = dat)
)

summary(m5_pSEM, .progressBar = F)

plot(m5_pSEM)

dSep(m5_pSEM)
fisherC(m5_pSEM)
coefs(m5_pSEM, standardize = "none")
residuals(m5_pSEM)


#model now indicates indirect effect of treatment on crustacean abundance, mediated by LSA
#can I do this? only transform one variable??
