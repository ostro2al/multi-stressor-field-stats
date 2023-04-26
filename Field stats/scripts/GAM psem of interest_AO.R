#Multi-stressor field experiment
#GAM SEMs
#A. Ostrowski  CJ Brown 
#March 3, 2023


library(mgcv)
library(dagitty)

#load csv
#dat <- read.csv("../Data/DataStacked_Control.csv")
dat <- read.csv("data/Data_Control.csv")
str(dat)
dat <- data.frame(unclass(dat), stringsAsFactors = TRUE)
str(dat)
dat$Plot.ID <- as.factor(dat$Plot.ID)
dat$Block <- as.factor(dat$Block)
dat$Week <- as.numeric(dat$Week)
str(dat)

sort(dat$Week)

#model of interest
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
plot(m3)

#using GAMs in SEMs
#R function to calculate the AIC, loglikelihood (Shipley & Douma 2020)
MLX2<-function(submodels,saturated.submodels,data){
  nobs<-dim(data)[1]
  nsub<-length(submodels)
  error.mes1<-error.mes2<-rep(FALSE,nsub)
  for(i in 1:nsub){
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
  model.AIC<-sum(out$AIC)
  model.LL<-sum(out$logLikelihoods)
  model.df<-sum(out.saturated$n.free.parameters)-
    sum(out$n.free.parameters)
  n.free.parameters<-sum(out$n.free.parameters)
  n.saturated.free.parameters<-sum(out.saturated$n.free.parameters)
  saturated.model.AIC<-sum(out.saturated$AIC)
  saturated.model.LL<-sum(out.saturated$logLikelihoods)
  saturated.model.df<-sum(out.saturated$k)
  X2<--2*(model.LL-saturated.model.LL)
  if(X2<0)X2<-NA
  null.prob<-NA
  if(!is.na(X2)&(model.df>0))null.prob<-1-pchisq(X2,model.df)
  error.flag<-sum(error.mes1)+sum(error.mes2)
  list(model.X2=X2,model.df=model.df,null.prob=null.prob,
       model.loglikelihood=model.LL,
       n.free.parameters=n.free.parameters,
       n.free.parameters.saturated=n.saturated.free.parameters,
       saturated.model.loglikelihood=saturated.model.LL,
       model.AIC=model.AIC,submodels=submodels,error.flag=error.flag)
}

#GAM SEM of interest (lowest AIC)
#model 3
m3 <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Treatment + s(Week, by = Treatment, k = 4) + 
        s(Block, bs = "re") + offset(log(Day0_density)),
      family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Treatment + s(Week, by = Treatment, k = 4) + 
        s(Block, bs = "re") + offset(Day0_Avg_LSA),
      data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
        s(Week, by = Treatment, k = 4) + s(Shoot_density, k = 4) + Avg_LSA +
        s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
      family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week, by = Block, k = 4) + Treatment + Avg_LSA + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) + Treatment + Crustacean_abundance + s(Week, by = Stressor_app, k = 4) + 
          s(Block, bs = "re"), data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + Avg_LSA + s(Shoot_density, k = 4) + 
          s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), family = nb(), data = dat)),
     data = dat)
#3 submodels 
m3.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
              Treatment + s(Week, by = Treatment, k = 4) + 
              s(Block, bs = "re") + offset(log(Day0_density)),
            family = "poisson", data = dat)
m3.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
              Treatment + s(Week, by = Treatment, k = 4) + 
              s(Block, bs = "re") + offset(Day0_Avg_LSA),
            data = dat)
m3.3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
              s(Week, by = Treatment, k = 4) + s(Shoot_density, k = 4) + Avg_LSA +
              s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
            family = nb(), data = dat)

#determine direct significant causal paths
summary(m3.1)
summary(m3.2)
summary(m3.3)

submodels<-list(m1=m3.1,m2=m3.2,m3=m3.3)
#likelihood values and AIC statistics
mod3<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                 k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod3$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod3[,3:4]<-AIC(m3.1,m3.2,m3.3)
mod3
#Overall AIC and likelihood for model 1:
mod3.AIC<-sum(mod3$AIC)
mod3.AIC
mod3.LL<-LL.fit<-sum(mod3$logLikelihoods)
mod3.LL


#other GAM SEMs tested
#model 1
m1 <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Treatment + s(Week, by = Treatment, k = 4) + 
        s(Block, bs = "re") + offset(log(Day0_density)),
      family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Treatment + s(Week, by = Treatment, k = 4) + 
        s(Block, bs = "re") + offset(Day0_Avg_LSA),
      data = dat),
  gam(Crustacean_abundance ~ s(Shoot_density, k = 4) + Avg_LSA +
        s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
      family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week,by = Block, k = 4) + Treatment + Avg_LSA + Crustacean_abundance + s(Week, by = Treatment, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) + Treatment + Crustacean_abundance + s(Week, by = Treatment, k = 4) + 
          s(Block, bs = "re"), data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + Avg_LSA + 
          s(Shoot_density, k = 4) + s(Week, by = Treatment, k = 4) + s(Block, bs = "re"), 
        family = nb(), data = dat)),
  data=dat)
#3 submodels
m1.1 <-  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
               Treatment + s(Week, by = Treatment, k = 4) + 
               s(Block, bs = "re") + offset(log(Day0_density)),
             family = "poisson", data = dat)
m1.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
              Treatment + s(Week, by = Treatment, k = 4) + 
              s(Block, bs = "re") + offset(Day0_Avg_LSA),
            data = dat)
m1.3 <-  gam(Crustacean_abundance ~ s(Shoot_density, k = 4) +
             Avg_LSA + s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
             family = nb(), data = dat)

#determine direct significant paths
summary(m1.1)
summary(m1.2)
summary(m1.3)

submodels<-list(m1=m1.1,m2=m1.2,m3=m1.3)
#likelihood values and AIC statistics
mod1<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                 k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod1$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod1[,3:4]<-AIC(m1.1,m1.2,m1.3)
mod1
#Overall AIC and likelihood for model 1:
mod1.AIC<-sum(mod1$AIC)
mod1.AIC
mod1.LL<-LL.fit<-sum(mod1$logLikelihoods)
mod1.LL


#model 2
m2 <- MLX2(submodels=list(
  gam(Shoot_density ~ s(Week, by = Block, k = 4) +
        Treatment + s(Week, by = Treatment, k = 4) + 
        s(Block, bs = "re") + offset(log(Day0_density)),
      family = "poisson", data = dat),
  gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
        Treatment + s(Week, by = Treatment, k = 4) + 
        s(Block, bs = "re") + offset(Day0_Avg_LSA),
      data = dat),
  gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + s(Shoot_density, k = 4) + Avg_LSA +
        s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
      family = nb(), data = dat)),
  saturated.submodels=list(
    gam(Shoot_density ~ s(Week, by = Block, k = 4) + Treatment + Avg_LSA + Crustacean_abundance + s(Week, by = Treatment, k = 4) + 
          s(Block, bs = "re"), family = "poisson", data = dat),
    gam(Avg_LSA ~ s(Week, by = Block, k = 4) + Treatment + Crustacean_abundance + s(Week, by = Treatment, k = 4) + 
          s(Block, bs = "re"), data = dat),
    gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + Avg_LSA + 
          s(Shoot_density, k = 4) + s(Week, by = Treatment, k = 4) + s(Block, bs = "re"), 
        family = nb(), data = dat)),
  data = dat)
#3 submodels 
m2.1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
              Treatment + s(Week, by = Treatment, k = 4) + 
              s(Block, bs = "re") + offset(log(Day0_density)),
            family = "poisson", data = dat)
m2.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
              Treatment + s(Week, by = Treatment, k = 4) + 
              s(Block, bs = "re") + offset(Day0_Avg_LSA),
            data = dat)
m2.3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + s(Shoot_density, k = 4) + Avg_LSA +
              s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
            family = nb(), data = dat)

#determine direct significant paths
summary(m2.1)
summary(m2.2)
summary(m2.3)

submodels<-list(m1=m2.1,m2=m2.2,m3=m2.3)
#likelihood values and AIC statistics
mod2<-data.frame(submodel=1:3,logLikelihoods=rep(NA,3),
                 k=rep(NA,3),AICs=rep(NA,3))
#get likelihoods
for(i in 1:3){
  mod2$logLikelihoods[i]<-logLik(submodels[[i]])
}
#get AIC and k
mod2[,3:4]<-AIC(m2.1,m2.2,m2.3)
mod2
#Overall AIC and likelihood for model 1:
mod2.AIC<-sum(mod2$AIC)
mod2.AIC
mod2.LL<-LL.fit<-sum(mod2$logLikelihoods)
mod2.LL
