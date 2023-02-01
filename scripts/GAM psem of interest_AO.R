#Multi-stressor field experiment
#Mechanistic Analysis - using GAM SEMs
#A. Ostrowski Feb 1, 2023


library(dplyr)
library(mgcv)
library(MASS)
library(dagitty)

#load csv
dat <- read.csv("DataStacked_Control.csv")
str(dat)
dat$Stressor_app <- as.numeric(dat$Stressor_app)
dat$Week <- as.numeric(dat$Week)
dat$Block <- as.numeric(dat$Block)
str(dat)

sort(dat$Week)

#hypothesised causal path diagram
m1 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Week -> Crustacean_abundance
                Treatment -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
plot(m1)

#determine independence statements in model - basis set of regressions needed to ensure model is fully tested
latents(m1) <- c() # undo out setting of X as latent from (f)
impliedConditionalIndependencies(m1, type="basis.set")

#individual regressions to test for d-sep
coef(summary(lm(Avg_LSA ~ Treatment + Week + Shoot_density, data = dat)))
coef(summary(lm(Shoot_density ~ Treatment + Week + Avg_LSA, data = dat)))
coef(summary(lm(Week ~ Treatment, data = dat)))

#
#trying GAMs in pSEM - Shipley & Douma 2020 method
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
m1a <- MLX2(submodels=list(
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
m1a.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
               Stressor_app + Avg_LSA + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), family = "poisson", data = dat)
m1a.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
               Stressor_app + s(Week, by = Stressor_app, k = 4) + 
               s(Block, bs = "re"), data = dat)
m1a.3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Stressor_app + Avg_LSA + 
               Shoot_density + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
             family = nb(), data = dat)

#determine direct significant paths
summary(m1a.1)
summary(m1a.2)
summary(m1a.3)

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

#NOTE: I also looked at reverse direction for relationship between LSA and density, but 
#model above had lowest AIC
