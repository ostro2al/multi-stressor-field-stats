#check raw data for LSA effect on crustaceans
#plot LSA by crustacean abundance at week 6, facet by treatment

library(rstatix)
library(dplyr)
library(ggplot2)
library(cowplot)

#read csv
dat <- read.csv("data/Data_Control.csv")

#plot crustacean abundance vs LSA to verify weak mediation effect
dat %>%
  ggplot()+
  (aes(x = Avg_LSA, y = Crustacean_abundance, colour = Treatment))+
  facet_grid(Treatment~Week)+
  geom_point(size = 2.5, position = position_dodge(0.8)) +
  stat_smooth(se = FALSE)+
  theme_cowplot()+
  xlab("Leaf surface area (cm^2)")+
  ylab("Crustacean abundance")
