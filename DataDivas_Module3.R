
#1) establishing the data tibble

library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)
anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")

# merge & filter
anole2 <- anole%>%
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%
  na.omit()

# log transformation 
anole.log <- anole2%>%
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)

#residuals column for linear model 
anole.log.lm  <- lm(HTotal~SVL,anole.log)
anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm))

# loading + plotting phylogenetic tree
anole.tree <- read.tree("anole.tre")
plot(anole.tree,cex=0.4)

# PGLS under Brownian Motion (BM) 
pgls.BM2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

# phylogenetically corrected residuals 
anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM2))

# anole.log tibble with residuals
print(anole.log)


#2) lm of ArbPD/PH on HTotal

# lm assessing ArbPD on HTotal
model_pd <- lm(HTotal~SVL+ArbPD, data=anole.log)
summary(model_pd)
# lm assessing PH on HTotal
model_ph <-  lm(HTotal~SVL+PH, data=anole.log)
summary(model_ph)


#3) plot the residuals of both lm against Ecomorph2

# mutating anole.log to include residuals from both models
anole.log <- anole.log %>%
  mutate(res_pd = residuals(model_pd),res_ph = residuals(model_ph))

# plotting the residuals from model_pd against ecomorph
pd_plot <-  anole.log %>% 
  ggplot(aes(x=Ecomorph2, y=res_pd))+
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", size=3, color="blue")+
  geom_jitter(width =0.3, alpha = 0.5) +
  theme_minimal()
print(pd_plot) 

# plotting residuals from model_ph against ecomorph
ph_plot <- anole.log %>% 
  ggplot(aes(x=Ecomorph2, y= res_ph))+
  geom_boxplot()+
  stat_summary(fun=mean, geom="point",size=3,color="red")+
  geom_jitter(width =0.3, alpha = 0.5) +
  theme_minimal()
print(ph_plot)


#4) PGLS HTotal-SVL relationship covarieties
anole.tree <- read.tree("anole.tre")

#hind limb-SVL + perch height
pgls_height <- gls(HTotal ~ SVL + PH, 
                   correlation = corBrownian(1, phy = anole.tree, form = ~Species), 
                   data = anole.log, 
                   method = "ML")
summary(pgls_height)

#hindlimb-SVL+ perch diameter
pgls_diameter <- gls(HTotal ~ SVL + ArbPD, 
                     correlation = corBrownian(1, phy = anole.tree, form = ~Species), 
                     data = anole.log, 
                     method = "ML")
summary(pgls_diameter)

# hindlimb-SVL + both perch height & perch diameter
pgls_combined <- gls(HTotal ~ SVL + PH + ArbPD, 
                     correlation = corBrownian(1, phy = anole.tree, form = ~Species), 
                     data = anole.log, 
                     method = "ML")
summary(pgls_combined)

#5) testing model fitness with AICc and AICw
# AICc
anole.pgls.aic <- AICc(pgls_height, pgls_diameter, pgls_combined)
print(anole.pgls.aic)
# height: -64.77956, diameter: -73.81081, combined:-75.52571

# AICw
anole.pgls.aicw <- aicw(anole.pgls.aic$AICc)
print(anole.pgls.aicw)
# height: .003247185, diameter: .296905077, both: .699847738

# Based on AIC scores, the pgls_combined model best fit our data since it has the lowest AICc score and the highest AICw score, therefore
# both covarities contribute to hindlimb length variation

summary(pgls_height)
summary(pgls_diameter)
summary(pgls_combined)
# Intercept:.3889, SVL: .000, PH: .0505, ArbPD: .0005
# based on the p values, SVL and perch diameter are significant predictors
# perch height's effect is less certain, but has marginal significance 
# the combined model is the best fit to our data

#7) produce a plot visualizing the effect of your covariate(s) and factors on the hindlimb residuals of the best fitting PGLS model.
anole.log <- anole.log %>%
  mutate(pgls_combined_res = residuals(pgls_combined))

# combined plot
combo_plot <- ggplot(anole.log, aes(x = ArbPD, y = pgls_combined_res, color = Ecomorph2))+
  geom_point()+  
  geom_smooth(method = "lm", se=FALSE) +
  facet_wrap(~Ecomorph2, scales = "free_x")+
  theme_minimal()

print(combo_plot)

