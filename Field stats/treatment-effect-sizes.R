
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


