#Multi-stressor field experiment
#Mechanistic Analysis - using GAM SEMs
#A. Ostrowski January 30, 2023


library(dplyr)
library(nlme) 
library(lme4)
library(mgcv)
library(piecewiseSEM)
library(MASS)
library(dagitty)

#load csv with data stacked to include day 0 data to see difference in model outputs
dat <- read.csv("DataStacked_Control.csv")
str(dat)
dat$Stressor_app <- as.numeric(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
dat$Block <- as.numeric(dat$Block)
str(dat)

sort(dat$Week)

#try different pathway hypotheses - pSEM using GAMs
#models of interest
DAG1 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Week -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
plot(DAG1)

DAG2 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
plot(DAG2)

DAG3 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Week -> Crustacean_abundance
                Treatment -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
plot(DAG3)

#
#d-sep tests for each model of interest
m1 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Week -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )

#test implied independencies of model against correlation matrix
corr <- lavCor(dat)
localTests(m1, sample.cov=corr, sample.nobs=nrow(dat)) 
#                                     estimate   p.value        2.5%     97.5%
#  Crs_ _||_ Trtm | A_LS, Sht_, Week 0.1039625 0.1841698 -0.04960906 0.2528021
#Trtm _||_ Week                    0.0000000 1.0000000 -0.15141514 0.1514151
#data do not provide evidence against implied conditional independence being tested

plotLocalTestResults(localTests(m1, sample.cov=corr, sample.nobs=nrow(dat)))

#d-sep tests
coef(summary(lm(Crustacean_abundance ~ Treatment + Avg_LSA + Shoot_density, dat)))
coef(summary(lm(Week ~ Treatment, dat)))

#determine independence statements in model - basis set of regressions needed to ensure model is fully tested
latents(m1) <- c() # undo out setting of X as latent from (f)
impliedConditionalIndependencies(m1, type="basis.set")
#A_LS _||_ Sht_ | Trtm, Week
#Crs_ _||_ Trtm | A_LS, Sht_, Week
#Sht_ _||_ A_LS | Trtm, Week
#Trtm _||_ Week
#Week _||_ Trtm   
#note last two are the same so there are 4 independence statements

#individual regressions to test for d-sep
coef(summary(lm(Avg_LSA ~ Treatment + Week + Shoot_density, data=dat)))
coef(summary(lm(Crustacean_abundance ~ Avg_LSA + Shoot_density + Week + Treatment, data=dat)))
coef(summary(lm(Shoot_density ~ Treatment + Week + Avg_LSA, data=dat)))
coef(summary(lm(Week ~ Treatment, data=dat)))

#repeat process for DAGs 2 and 3
m2 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
#test implied independencies of model against correlation matrix
corr <- lavCor(dat)
localTests(m2, sample.cov=corr, sample.nobs=nrow(dat)) 
#                              estimate    p.value        2.5%       97.5%
#Crs_ _||_ Trtm | A_LS, Sht_  0.1154297 0.13879155 -0.03755204  0.26321103
#Crs_ _||_ Week | A_LS, Sht_ -0.1892739 0.01444585 -0.33230314 -0.03804943
#Trtm _||_ Week               0.0000000 1.00000000 -0.15141514  0.15141514
#data provide evidence against implied conditional independence being tested for week effect 
#on crustaceans mediated by LSA and density - should include in model

plotLocalTestResults(localTests(m2, sample.cov=corr, sample.nobs=nrow(dat)))

#d-sep tests
coef(summary(lm(Crustacean_abundance ~ Treatment + Avg_LSA + Shoot_density, dat)))
coef(summary(lm(Crustacean_abundance ~ Week + Avg_LSA + Shoot_density, dat)))
coef(summary(lm(Week ~ Treatment, dat)))

#determine independence statements in model - basis set of regressions needed to ensure model is fully tested
latents(m2) <- c() # undo out setting of X as latent from (f)
impliedConditionalIndependencies(m2, type="basis.set")
#A_LS _||_ Sht_ | Trtm, Week
#Crs_ _||_ Trtm, Week | A_LS, Sht_
#Sht_ _||_ A_LS | Trtm, Week
#Trtm _||_ Week
#Week _||_ Trtm
#note last two are the same so there are 4 independence statements

#individual regressions to test for d-sep
coef(summary(lm(Avg_LSA ~ Treatment + Week + Shoot_density, data = dat)))
coef(summary(lm(Crustacean_abundance ~ Avg_LSA + Shoot_density + Treatment + Week, data = dat)))
coef(summary(lm(Shoot_density ~ Treatment + Week + Avg_LSA, data = dat)))
coef(summary(lm(Week ~ Treatment, data = dat)))

m3 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Week -> Crustacean_abundance
                Treatment -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
#test implied independencies of model against correlation matrix
corr <- lavCor(dat)
localTests(m3, sample.cov=corr, sample.nobs=nrow(dat)) 
#             estimate p.value       2.5%     97.5%
#Trtm _||_ Week        0       1 -0.1514151 0.1514151
#data do not provide evidence against implied conditional independence being tested

plotLocalTestResults(localTests(m3, sample.cov=corr, sample.nobs=nrow(dat)))

#d-sep tests
coef(summary(lm(Week ~ Treatment, dat)))

#determine independence statements in model - basis set of regressions needed to ensure model is fully tested
latents(m3) <- c() # undo out setting of X as latent from (f)
impliedConditionalIndependencies(m3, type="basis.set")
#A_LS _||_ Sht_ | Trtm, Week
#Sht_ _||_ A_LS | Trtm, Week
#Trtm _||_ Week
#Week _||_ Trtm
#note last two are the same so there are 3 independence statements

#individual regressions to test for d-sep
coef(summary(lm(Avg_LSA ~ Treatment + Week + Shoot_density, data = dat)))
coef(summary(lm(Shoot_density ~ Treatment + Week + Avg_LSA, data = dat)))
coef(summary(lm(Week ~ Treatment, data = dat)))

#model should include tests for:
#effect of week, treatment & density on LSA
#effect of week, treatment & LSA on density
#effect of week and LSA on crustacean abundance


#
#trying GAMs in pSEM - Shipley & Douma 2020 method
#copy and paste r code from paper
#R function to calculate the AIC, loglikelihood and chi-
#square statistics.
MLX2<-function(submodels,saturated.submodels,data){
  #submodels is a list containing the submodels in your model
  # in the form of linear, generalized linear, generalized additive
  # mixed models, or any other model object that has AIC and
  #log-likelihood attributes.
  #saturated.submodels is a list containing the submodels of your
  # model that defines the saturated submodels (or otherwise)
  #into which your model is properly nested.
  #data is the data set
  #
  #number of submodels in full model
  nobs<-dim(data)[1]
  nsub<-length(submodels)
  error.mes1<-error.mes2<-rep(FALSE,nsub)
  #if there is an error in estimating a model, then error.mes==TRUE
  #and don't calculate statistics
  for(i in 1:nsub){
    #check if the submodels and the saturated submodels are
    #of the same class, and stop if not.
    stopifnot(class(submodels[[i]])==
                class(saturated.submodels[[i]]))
    if(all(class(submodels[[i]])=="lm")){
      error.mes1[i]<-error.mes2[i]<-FALSE
    }
    else{
      error.mes1[i]<-!submodels[[i]]$converged & submodels[[i]]$boundary
      error.mes2[i]<-!saturated.submodels[[i]]$converged & saturated.submodels[[i]]$boundary
    }}
  out<-data.frame(submodel=1:nsub,logLikelihoods=rep(NA,nsub),
                  k=rep(NA,nsub),AICs=rep(NA,nsub),n.free.parameters=
                    rep(NA,nsub))
  out.saturated<-data.frame(submodel=1:nsub,logLikelihoods=rep(NA,nsub),
                            k=rep(NA,nsub),AICs=rep(NA,nsub),
                            n.free.parameters=rep(NA,nsub))
  
  #get likelihoods, AIC & k and store in "out"
  for(i in 1:nsub){
    out$logLikelihoods[i]<-logLik(submodels[[i]])
    out$AICs[i]<-AIC(submodels[[i]])
    out$k[i]<-nobs-df.residual(submodels[[i]])
    out$n.free.parameters[i]<-attributes(logLik(submodels[[i]]))$df
    out.saturated$n.free.parameters[i]<-attributes(logLik(saturated.submodels[[i]]))$df
    out.saturated$logLikelihoods[i]<-logLik(saturated.submodels[[i]])
    out.saturated$AICs[i]<-AIC(saturated.submodels[[i]])
    out.saturated$k[i]<-nobs-df.residual(saturated.submodels[[i]])
  }
  #Overall k, AIC and likelihood for models:
  model.AIC<-sum(out$AIC)
  model.LL<-sum(out$logLikelihoods)
  #model df is the difference in the # of free parameters in the
  #less constrained model relative to the # of free parameters in the
  # more constrained (nested) model
  model.df<-sum(out.saturated$n.free.parameters)-
    sum(out$n.free.parameters)
  n.free.parameters<-sum(out$n.free.parameters)
  n.saturated.free.parameters<-sum(out.saturated$n.free.parameters)
  saturated.model.AIC<-sum(out.saturated$AIC)
  saturated.model.LL<-sum(out.saturated$logLikelihoods)
  saturated.model.df<-sum(out.saturated$k)
  # the MLX2 statistic is the difference in likelihoods between the
  #more constrained (nested) model and the less constrained model;
  #usually the saturated model
  X2<--2*(model.LL-saturated.model.LL)
  if(X2<0)X2<-NA
  #  df<-saturated.model.df-model.df
  null.prob<-NA
  #Only calculate null prob if the X2 is valid with valid df
  if(!is.na(X2)&(model.df>0))null.prob<-1-pchisq(X2,model.df)
  #check if any models had errors in estimation
  error.flag<-sum(error.mes1)+sum(error.mes2)
  list(model.X2=X2,model.df=model.df,null.prob=null.prob,
       model.loglikelihood=model.LL,
       n.free.parameters=n.free.parameters,
       n.free.parameters.saturated=n.saturated.free.parameters,
       saturated.model.loglikelihood=saturated.model.LL,
       model.AIC=model.AIC,submodels=submodels,error.flag=error.flag)
}

#now run pSEM on DAGs from above
#using the gam smoother for both the substantive model and its saturated version
#model 1
m1a <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Stressor_app + Avg_LSA + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Stressor_app + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "gaussian", data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Avg_LSA +
        Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
    family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week, by = Block, k = 4) +
          Stressor_app + Avg_LSA + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
          Stressor_app + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), family = "gaussian", data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
          Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
        family = nb(), data = dat)),
  data = dat)

#3 submodels - model 1a
m1a.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
               Stressor_app + Avg_LSA + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), family = "poisson", data = dat)
m1a.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
               Stressor_app + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), data = dat)
m1a.3 <-  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Avg_LSA + 
                Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
              family = nb(), data = dat)

submodels<-list(m1=m1a.1,m2=m1a.2,m3=m1a.3)
#likelihood values and AIC statistics
mod1a<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                  k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod1a$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod1a[,3:4]<-AIC(m1a.1,m1a.2,m1a.3)
mod1a
#Overall AIC and likelihood for model 1:
mod1a.AIC<-sum(mod1a$AIC)
mod1a.AIC
mod1a.LL<-LL.fit<-sum(mod1a$logLikelihoods)
mod1a.LL

#switching the direction of the arrow between density and LSA since have the same causal order
m1b <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Stressor_app + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Stressor_app + Shoot_density + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Avg_LSA + 
        Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
      family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week, by = Block, k = 4) +
          Stressor_app + Avg_LSA + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
          Stressor_app + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
          Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
        family = nb(), data = dat)),
  data = dat)

#3 submodels - model 1b
m1b.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
               Stressor_app + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), family = "poisson", data = dat)
m1b.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
               Stressor_app + Shoot_density + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), data = dat)
m1b.3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Avg_LSA + 
               Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
             family = nb(), data = dat)

submodels<-list(m1=m1b.1,m2=m1b.2,m3=m1b.3)
#likelihood values and AIC statistics
mod1b<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                  k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod1b$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod1b[,3:4]<-AIC(m1b.1,m1b.2,m1b.3)
mod1b
#Overall AIC and likelihood for model 1:
mod1b.AIC<-sum(mod1b$AIC)
mod1b.AIC
mod1b.LL<-LL.fit<-sum(mod1b$logLikelihoods)
mod1b.LL

  
#model 2  
m2a <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Stressor_app + Avg_LSA + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Stressor_app + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), data = dat),
  gam(Crustacean_abundance ~ Avg_LSA + 
        Shoot_density + + s(Block, bs = "re"), 
      family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week, by = Block, k = 4) +
          Stressor_app + Avg_LSA + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
          Stressor_app + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
          Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
        family = nb(), data = dat)),
  data=dat)


#3 submodels - model 2a
m2a.1 <-  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
                Stressor_app + Avg_LSA + s(Week, by = Stressor_app, k = 4) + 
                s(Block, bs = "re"), family = "poisson", data = dat)
m2a.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
               Stressor_app + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), data = dat)
m2a.3 <-  gam(Crustacean_abundance ~ Avg_LSA + 
                Shoot_density + + s(Block, bs = "re"), 
              family = nb(), data = dat)

submodels<-list(m1=m2a.1,m2=m2a.2,m3=m2a.3)
#likelihood values and AIC statistics
mod2a<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                  k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod2a$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod2a[,3:4]<-AIC(m2a.1,m2a.2,m2a.3)
mod2a
#Overall AIC and likelihood for model 1:
mod2a.AIC<-sum(mod2a$AIC)
mod2a.AIC
mod2a.LL<-LL.fit<-sum(mod2a$logLikelihoods)
mod2a.LL

#switching the direction of the arrow between density and LSA since have the same causal order
m2b <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Stressor_app + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Stressor_app + Shoot_density + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), data = dat),
  gam(Crustacean_abundance ~ Avg_LSA + 
        Shoot_density + + s(Block, bs = "re"), 
      family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week, by = Block, k = 4) +
          Stressor_app + Avg_LSA + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
          Stressor_app + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
          Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
        family = nb(), data = dat)),
  data=dat)

#3 submodels - model 2b
m2b.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
               Stressor_app + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), family = "poisson", data = dat)
m2b.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
               Stressor_app + Shoot_density + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), data = dat)
m2b.3 <- gam(Crustacean_abundance ~ Avg_LSA + 
               Shoot_density + + s(Block, bs = "re"), 
             family = nb(), data = dat)

submodels<-list(m1=m2b.1,m2=m2b.2,m3=m2b.3)
#likelihood values and AIC statistics
mod2b<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                  k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod2b$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod2b[,3:4]<-AIC(m2b.1,m2b.2,m2b.3)
mod2b
#Overall AIC and likelihood for model 1:
mod2b.AIC<-sum(mod2b$AIC)
mod2b.AIC
mod2b.LL<-LL.fit<-sum(mod2b$logLikelihoods)
mod2b.LL


#model 3
m3a <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
      Stressor_app + Avg_LSA + s(Week, by = Stressor_app, k = 4) + 
      s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
      Stressor_app + s(Week, by = Stressor_app, k = 4) + 
      s(Block, bs = "re"), data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
      Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
      family = nb(), data = dat)),
  saturated.submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Stressor_app + Avg_LSA + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Stressor_app + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
        Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
      family = nb(), data = dat)),
  data=dat)

#3 submodels - model 3a
m3a.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
               Stressor_app + Avg_LSA + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), family = "poisson", data = dat)
m3a.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
               Stressor_app + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), data = dat)
m3a.3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
               Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
             family = nb(), data = dat)

submodels<-list(m1=m3a.1,m2=m3a.2,m3=m3a.3)
#likelihood values and AIC statistics
mod3a<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                  k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod3a$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod3a[,3:4]<-AIC(m3a.1,m3a.2,m3a.3)
mod3a
#Overall AIC and likelihood for model 1:
mod3a.AIC<-sum(mod3a$AIC)
mod3a.AIC
mod3a.LL<-LL.fit<-sum(mod3a$logLikelihoods)
mod3a.LL

#switching the direction of the arrow between density and LSA since have the same causal order
m3b <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Stressor_app + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Stressor_app + Shoot_density + s(Week, by = Stressor_app, k = 4) + 
        s(Block, bs = "re"), data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
        Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
      family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week, by = Block, k = 4) +
          Stressor_app + Avg_LSA + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
          Stressor_app + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
          Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
        family = nb(), data = dat)),
  data=dat)

#3 submodels - model 3b
m3b.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
              Stressor_app + s(Week, by = Stressor_app, k = 4) + 
              s(Block, bs = "re"), family = "poisson", data = dat)
m3b.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
               Stressor_app + Shoot_density + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), data = dat)
m3b.3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
               Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
             family = nb(), data = dat)

submodels<-list(m1=m3b.1,m2=m3b.2,m3=m3b.3)
#likelihood values and AIC statistics
mod3b<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                 k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod3b$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod3b[,3:4]<-AIC(m3b.1,m3b.2,m3b.3)
mod3b
#Overall AIC and likelihood for model 1:
mod3b.AIC<-sum(mod3b$AIC)
mod3b.AIC
mod3b.LL<-LL.fit<-sum(mod3b$logLikelihoods)
mod3b.LL

#compare overall AIC across all models - mod1b and mod3b have lowest AIC
#mod1a
#[1] 5349.437
#mod1b
#[1] 5905.095
#mod2a
#[1] 5360.266
#mod2b
#[1] 5915.924
#mod3a
#[1] 5349.437
#mod3b
#[1] 5905.095