
library(tidyverse)
library(mgcv)
library(visreg)
library(cowplot)

#read csv
dat <- read.csv("data/Data_Control.csv")
str(dat)

#
# Parameters
#

nsims <- 1000
k_val <- 4
dat$Treatment <- factor(dat$Treatment, levels = c("Control", "StaticStatic", "InPhase", "OutPhase"))

#
# Shoot density effect sizes 
#

m3.1 <- gam(Shoot_density ~ s(Week, by = Block, k = k_val) +
              Treatment + s(Week, by = Treatment, k = k_val) + 
              s(Block, bs = "re") + offset(log(Day0_density)),
            family = "poisson", data = dat)

#
# Create new dataframe for predictions 
#

newd <- with(dat, data.frame(Week = 6,
                   Block = 1,
                   Treatment = unique(Treatment),
                   Day0_density = mean(Day0_density))
)

newd$pred <- predict(m3.1, newd, type = "response")
newd

#get code from ?predict.gam
Xp <- predict(m3.1, newd, type="lpmatrix") 

rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

br <- rmvn(nsims, coef(m3.1),m3.1$Vp) ## 1000 replicate param. vectors
res <- matrix(0,ncol = nsims, nrow = nrow(newd),nsims)
res3 <- res2 <- rep(0, nsims)
#static relative to in phase
for (i in 1:1000){ 
  pr <- Xp %*% br[i,] ## replicate predictions
  res[,i] <- exp(pr)
  res2[i] <- pr[2] < pr[3]
  res3[i] <- exp(pr[2]) - exp(pr[3])
}

quantile(res3, probs = c(0.025, 0.5, 0.975))

sum(res2)/nsims
cbind(newd,t(apply(res, 1, quantile, probs = c(0.025, 0.5, 0.975))))

mean(res);var(res)

#static relative to out phase
for (i in 1:1000){ 
  pr <- Xp %*% br[i,] ## replicate predictions
  res[,i] <- exp(pr)
  res2[i] <- pr[2] < pr[4]
  res3[i] <- exp(pr[2]) - exp(pr[4])
}

quantile(res3, probs = c(0.025, 0.5, 0.975))

sum(res2)/nsims
cbind(newd,t(apply(res, 1, quantile, probs = c(0.025, 0.5, 0.975))))

mean(res);var(res)


#
# Surface area effect sizes 
#
m3.2 <- gam(Avg_LSA ~ s(Week, by = Block, k = 4) +
              Treatment + s(Week, by = Treatment, k = 4) + 
              s(Block, bs = "re") + offset(Day0_Avg_LSA),
            data = dat)

#
# Create new dataframe for predictions 
#

newd2 <- with(dat, data.frame(Week = 6,
                             Block = 1,
                             Treatment = unique(Treatment),
                             Day0_Avg_LSA = mean(Day0_Avg_LSA))
)

newd2$pred <- predict(m3.2, newd2, type = "response")
newd2

#get code from ?predict.gam
Xp <- predict(m3.2, newd2, type="lpmatrix") 

rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

br <- rmvn(nsims, coef(m3.2),m3.2$Vp) ## 1000 replicate param. vectors
res <- matrix(0,ncol = nsims, nrow = nrow(newd2),nsims)
res3 <- res2 <- rep(0, nsims)
#static relative to in phase
for (i in 1:1000){ 
  pr <- Xp %*% br[i,] ## replicate predictions
  res[,i] <- exp(pr)
  res2[i] <- pr[2] < pr[3]
  res3[i] <- exp(pr[2]) - exp(pr[3])
}

quantile(res3, probs = c(0.025, 0.5, 0.975))

sum(res2)/nsims
cbind(newd2,t(apply(res, 1, quantile, probs = c(0.025, 0.5, 0.975))))

mean(res);var(res)


#
# Crustacean abundance effect sizes 
#NOTE: this is the total treatment effect on crustaceans-removed LSA and density effect in pSEM
m3.3 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
              s(Week, by = Treatment, k = 4) +
              s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
            family = nb(), data = dat)

#
# Create new dataframe for predictions 
#

newd3 <- with(dat, data.frame(Week = 6,
                              Block = 1,
                              Treatment = unique(Treatment),
                              Day0_crustacean = mean(Day0_crustacean))
)

newd3$pred <- predict(m3.3, newd3, type = "response")
newd3

#get code from ?predict.gam
Xp <- predict(m3.3, newd3, type="lpmatrix") 

rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

br <- rmvn(nsims, coef(m3.3),m3.3$Vp) ## 1000 replicate param. vectors
res <- matrix(0,ncol = nsims, nrow = nrow(newd3),nsims)
res3 <- res2 <- rep(0, nsims)
#static relative to in phase
for (i in 1:1000){ 
  pr <- Xp %*% br[i,] ## replicate predictions
  res[,i] <- exp(pr)
  res2[i] <- pr[2] < pr[3]
  res3[i] <- exp(pr[2]) - exp(pr[3])
}

quantile(res3, probs = c(0.025, 0.5, 0.975))

sum(res2)/nsims
cbind(newd3,t(apply(res, 1, quantile, probs = c(0.025, 0.5, 0.975))))

mean(res);var(res)


#
# Crustacean abundance effect sizes 
#LSA and density effects on crustaceans
m3.4 <- gam(Crustacean_abundance ~ s(Week, by = Block, k = 4) + Treatment + 
              s(Week, by = Treatment, k = 4) + s(Shoot_density, k = 4) + Avg_LSA +
              s(Block, bs = "re") + offset(log(Day0_crustacean + 0.01)),
            family = nb(), data = dat)

#
# Create new dataframe for predictions 
#

newd4 <- with(dat, data.frame(Week = 6,
                              Block = 1,
                              Shoot_density = mean(Shoot_density, by = Treatment),
                              Avg_LSA = mean(Avg_LSA, by = Treatment),
                              Treatment = unique(Treatment),
                              Day0_crustacean = mean(Day0_crustacean))

