#Multi-stressor field study - Calculate effect sizes
#A. Ostrowski and C. Brown
#March 17, 2023

library(mgcv)
library(dagitty)
library(visreg)
library(tidyverse)


#read csv
dat <- read.csv("Data_Control.csv")
str(dat)
dat <- data.frame(unclass(dat), stringsAsFactors = TRUE)
str(dat)
dat$Block <- as.factor(dat$Block)
str(dat)

sort(dat$Week)
dat$Treatment <- relevel(factor(dat$Treatment), ref = "StaticStatic")

#treatment effect on shoot density
m1 <- gam(Shoot_density ~ s(Week, by = Block, k = 4) +
            Treatment + s(Week, by = Treatment, k = 4) + 
            s(Block, bs = "re") + offset(log(Day0_density)),
          family = "poisson", data = dat)

gam.check(m1)
plot(m1)
summary(m1)


#calculate effect sizes
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*%  matrix(rnorm(m * n), m, n))
}
#1000 draws of the parameters from covariance matrix 
br <- rmvn(1000, coef(m1), m1$Vp)

#predictions
xfit <-  visreg(m1, xvar = "Shoot_density", partial = FALSE, 
                plot = TRUE, by = "Treatment")
#newdat1 <- xfit$fit %>% filter((Shoot_density %% 1) == 0)
newdat1 <- xfit$fit
lp1 <- predict(m1, newdata = newdat1, type = "lpmatrix")


treat <- match(subset(newdat1, Treatment == 'Control')$Shoot_density, 
               subset(newdat1, Treatment == 'StaticStatic')$Shoot_density,
               subset(newdat1, Treatment == 'InPhase')$Shoot_density,
               subset(newdat1, Treatment == 'OutPhase')$Shoot_density)
control <-  which(newdat1$Treatment == 'Control')


#Simulate posterior
B20mult <-  matrix(NA, nrow = nrow(newdat1),
                   ncol = 1000)
res <- res2 <- matrix(NA, nrow = nrow(newdat1)/2,
                      ncol = 1000)
res3 <- res4 <- rep(NA, nrow(newdat1))


for (i in 1:1000){
  val1 <- lp1 %*% br[i,]  
  B20mult[, i] <- 10^val1
  #Treatment
  res[,i] <- val1[treat] < val1[control] #Is it higher under control?
  res2[,i] <- 10^(val1[control] - val1[treat]) #how much higher?
  res3[i] <- val1[1] < val1[5] #prob its higher 
  res4[i] <- 10^(val1[5])/10^(val1[1]) #how much higher? 
}

CItrend <- data.frame(t(apply(B20mult, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat1)
CItrend$Treatment <- factor(CItrend$Treatment, levels = c("Control","StaticStatic", "InPhase", "OutPhase"), labels = c("Control", "StaticStatic", "InPhase", "OutPhase"))
treat_probCIs <- data.frame(prob_treat = apply(res, 1, sum)/1000) %>%
  cbind(newdat1)
treat_diff_magCIs <- data.frame(t(apply(res2, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat1)

#effect multiples
quantile(res4, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

#density effect at treatment control and static
quantile(B20mult[6,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(B20mult[10,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

#plot results
g1 <- ggplot(CItrend) + 
  aes(x = Treatment, y = X50.) +
  #geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.5) + 
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5) + 
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 2, 
                 position = position_dodge(0.25)) +
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1, 
                 position = position_dodge(0.25)) +
  geom_point(size = 3, position = position_dodge(0.25)) +
  geom_hline(yintercept = 1)+ 
  theme_classic() + labs(colour="", fill="") +
  ylab("Change in shoot density") + 
  #xlab("Week") + scale_y_continuous(breaks = seq(0,40, 5), 
   #                                  limits = c(0, 10)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), 
        legend.position = c(.90, .90)) 
g1






#3 submodels 
m3.1 <- gam(Shoot_density ~ Stressor_app + s(Week, by = Stressor_app, k = 4) + 
              s(Block, bs = "re"), family = "poisson", data = dat)
m3.2 <- gam(Avg_LSA ~ Stressor_app + s(Week, by = Stressor_app, k = 4) + 
              s(Block, bs = "re"), data = dat)
m3.3 <- gam(Crustacean_abundance ~ Stressor_app + Avg_LSA +
              s(Shoot_density, k = 4) + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
            family = nb(), data = dat)
#m3.4 <- gam(Crustacean_abundance ~ Stressor_app + s(Week, by = Stressor_app, k = 4) + s(Block, bs = "re"), 
#family = nb(), data = dat)


#determine direct significant causal paths
summary(m3.1)
summary(m3.2)
summary(m3.3)
#summary(m3.4)


