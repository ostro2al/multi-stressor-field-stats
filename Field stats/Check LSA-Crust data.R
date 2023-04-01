#check raw data for LSA effect on crustaceans
#plot LSA by crustacean abundance at week 6, facet by treatment

library(rstatix)
library(dplyr)
library(ggplot2)

#read csv
dat <- read.csv("data/Data_Control.csv")

#plot raw data for week 6
dat %>%
  filter(Week == 6) %>%
  ggplot()+
  (aes(x = Avg_LSA, y = Crustacean_abundance))+
  facet_grid(~Treatment)+
  geom_point() +
  theme_cowplot()

#get mean data for LSA and crust
options(tibble.print_min = Inf)
dat2 <- dat %>% 
  group_by(Treatment, Week) %>%
  get_summary_stats(Avg_LSA, type = "mean_se")
names(dat2)[5] <- "LSA"
dat2 


dat3 <- dat %>% 
  group_by(Treatment, Week) %>%
  get_summary_stats(Crustacean_abundance, type = "mean_se")
names(dat3)[5] <- "Crust"
dat3 <- dat3[,!names(dat3) %in% c("Week", "Treatment", "variable", "n", "se")]
dat3  

dat4 <- cbind(dat2, dat3)
dat4

#mean data
dat4 %>%
  filter(Week == 6) %>%
  ggplot()+
  (aes(x = LSA, y = Crust, group = Treatment, colour = Treatment))+
  #facet_grid(~Treatment)+
  geom_point() +
  theme_cowplot()