#Multi-stressor field experiment
#DAGs - hypothesised causal pathways, data consistency, adjestment sets
#A. Ostrowski  CJ Brown 
#March 3, 2023


library(dagitty)
library(ggdag)

#read csv
dat <- read.csv("../Data/DataStacked_Control.csv")
str(dat)
dat$Stressor_app <- as.factor(dat$Stressor_app)
dat$Treatment <- as.factor(dat$Treatment)
dat$Week <- as.factor(dat$Week)
dat$Block <- as.factor(dat$Block)
str(dat)

sort(dat$Week)

#try different pathway hypotheses - models of interest
#model 1
m1 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Block -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Block -> Avg_LSA
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Block -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
plot(m1)

#d-sep tests
localTests(m1, dat, type="cis.chisq")

#inspect paths between variables
paths(m1, "Treatment", "Shoot_density", directed = FALSE)
paths(m1, "Week", "Shoot_density", directed = FALSE)
paths(m1, "Treatment", "Avg_LSA", directed = FALSE)
paths(m1, "Week", "Avg_LSA", directed = FALSE)
paths(m1, "Avg_LSA", "Crustacean_abundance", directed = FALSE)
paths(m1, "Shoot_density", "Crustacean_abundance", directed = FALSE)

#identify adjustment sets
adjustmentSets(m1, "Treatment", "Shoot_density", type="canonical")
adjustmentSets(m1, "Treatment", "Avg_LSA", type="canonical")
adjustmentSets(m1, "Week", "Shoot_density", type="canonical")
adjustmentSets(m1, "Week", "Avg_LSA", type="canonical")
adjustmentSets(m1, "Block", "Shoot_density", type="canonical")
adjustmentSets(m1, "Block", "Avg_LSA", type="canonical")
adjustmentSets(m1, "Block", "Crustacean_abundance", type="canonical")
adjustmentSets(m1, "Avg_LSA", "Crustacean_abundance", type="canonical")
adjustmentSets(m1, "Shoot_density", "Crustacean_abundance", type="canonical")


#model 2
m2 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Block -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Block -> Avg_LSA
                Week -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Block -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
plot(m2)

#d-sep tests
localTests(m2, dat, type="cis.chisq")

#inspect paths between variables
paths(m2, "Treatment", "Shoot_density", directed = FALSE)
paths(m2, "Week", "Shoot_density", directed = FALSE)
paths(m2, "Treatment", "Avg_LSA", directed = FALSE)
paths(m2, "Week", "Avg_LSA", directed = FALSE)
paths(m2, "Week", "Crustacean_abundance", directed = FALSE)
paths(m2, "Avg_LSA", "Crustacean_abundance", directed = FALSE)
paths(m2, "Shoot_density", "Crustacean_abundance", directed = FALSE)

#identify adjustment sets
adjustmentSets(m2, "Treatment", "Shoot_density", type="canonical")
adjustmentSets(m2, "Treatment", "Avg_LSA", type="canonical")
adjustmentSets(m2, "Week", "Shoot_density", type="canonical")
adjustmentSets(m2, "Week", "Avg_LSA", type="canonical")
adjustmentSets(m2, "Week", "Crustacean_abundance", type="canonical")
adjustmentSets(m2, "Block", "Shoot_density", type="canonical")
adjustmentSets(m2, "Block", "Avg_LSA", type="canonical")
adjustmentSets(m2, "Block", "Crustacean_abundance", type="canonical")
adjustmentSets(m2, "Avg_LSA", "Crustacean_abundance", type="canonical")
adjustmentSets(m2, "Shoot_density", "Crustacean_abundance", type="canonical")


#model 3
m3 <- dagitty( "dag {Treatment -> Shoot_density
                Week -> Shoot_density
                Block -> Shoot_density
                Treatment -> Avg_LSA
                Week -> Avg_LSA
                Block -> Avg_LSA
                Week -> Crustacean_abundance
                Treatment -> Crustacean_abundance
                Block -> Crustacean_abundance
                Avg_LSA -> Crustacean_abundance
                Shoot_density -> Crustacean_abundance
                Avg_LSA <-> Shoot_density
                }" )
plot(m3)

#d-sep tests
localTests(m3, dat, type="cis.chisq")

#inspect paths between variables
paths(m3, "Treatment", "Shoot_density", directed = FALSE)
paths(m3, "Week", "Shoot_density", directed = FALSE)
paths(m3, "Treatment", "Avg_LSA", directed = FALSE)
paths(m3, "Week", "Avg_LSA", directed = FALSE)
paths(m3, "Week", "Crustacean_abundance", directed = FALSE)
paths(m3, "Treatment", "Crustacean_abundance", directed = FALSE)
paths(m3, "Avg_LSA", "Crustacean_abundance", directed = FALSE)
paths(m3, "Shoot_density", "Crustacean_abundance", directed = FALSE)

#identify adjustment sets
adjustmentSets(m3, "Treatment", "Shoot_density", type="canonical")
adjustmentSets(m3, "Treatment", "Avg_LSA", type="canonical")
adjustmentSets(m3, "Week", "Shoot_density", type="canonical")
adjustmentSets(m3, "Week", "Avg_LSA", type="canonical")
adjustmentSets(m3, "Week", "Crustacean_abundance", type="canonical")
adjustmentSets(m3, "Treatment", "Crustacean_abundance", type="canonical")
adjustmentSets(m3, "Block", "Shoot_density", type="canonical")
adjustmentSets(m3, "Block", "Avg_LSA", type="canonical")
adjustmentSets(m3, "Block", "Crustacean_abundance", type="canonical")
adjustmentSets(m3, "Avg_LSA", "Crustacean_abundance", type="canonical")
adjustmentSets(m3, "Shoot_density", "Crustacean_abundance", type="canonical")
