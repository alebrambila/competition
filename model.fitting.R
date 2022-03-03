library(rstan)
library(reshape2)
library(brms)
library(tidybayes)
library(modelr)
library(nlstools)
library(tidyverse)
library(minpack.lm)
library(grid)
library(gridExtra)
library(ggpubr)


# Jeff notes:
# Alejandro: for seedlings, gaussian not working, but also should try binomial one - maybe more appropriate

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('data.cleaning.R')

# This script is to estimate the parameters necessary to be able to perform invasion analysis and additional 
# simulation of an annual (Lolium multiflorum) and perennial (Festuca roemeri) in competition at various densities.
#  Here I estimate the inherent growth rates and competition coefficients in both warmed and ambient plots.

###### 1. SIMPLE ANNUALS: SEEDS-IN:SEEDS-OUT ----
## annual modeling uses the tibble, 'annuals' from 'data.cleaning.r'. 

# seeded_a refers to the actual number of annual seeds added, regardless of gopher damage
#percap is the number seeds out/ seeds in
#seeded_am2 refers to the number of seeds added/m2
# pm2 refers to the m2 density of ACTUAL adult perrennials in the plot
#starting_pm2 refers to the density/m2 of adult perennials that should have been in the plot ignoring gophers - the only source of adult mortality over the experiment
# dpm2 refers to the number of damaged adult perennials per meter square
# nam2 refers to the number of 'new adults' per meter square, seedlings that were left as phytometers to see in the next year
#now the model runs off of pm2, but it could be adapted for dpm2, nam2, or starting pm2 

### visualize ----
ggplot(annuals, aes(x=seeded_am2, y=starting_pm2, color=warmtrt)) +
  geom_jitter(width=.1)+facet_wrap(~time)#+
 # geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
 # scale_colour_manual(values = c("dodgerblue", "darkred"))
 
# JD: look at data
head(annuals); dim(annuals)
# how many plots / block?  answer: 12
annuals %>% group_by(block) %>%
  summarize(n.plots = n_distinct(plotid))
table(annuals$warmtrt, annuals$block)

# per capita
ggplot(annuals, aes(x=starting_pm2, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=as.factor(time)), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x , se=F) +xlab("perennial adults")+ylab("annual percapita seeds out")

ggplot(subset(annuals), aes(x=seeded_am2, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=as.factor(time)), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x, se=F)+xlab("annual seeds in")#+ylab("annual percapita seeds out")

#added in seedling perennials (not sure it's going to be helpful)
#ggplot(subset(annuals), aes(x=seeded_s, y=percap, color=warmtrt)) +
#  geom_jitter(aes(shape=as.factor(time)), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
#  geom_smooth(method = 'lm',formula = y ~ x , se=F)+xlab("perennial seeds in")#+ylab("annual percapita seeds out")

### BRM fits ----
#Simple annual gaussian model (no block/alphaAS)
#annuals.sc <- mutate(annuals, seeded_a=seeded_a/1000)
annual.simple.gaussian <- brm(bf(percap ~ lambdaA*100 / (1 + alphaAA*seeded_am2 + alphaAP*starting_pm2), # + alphaAS*seeded_s #switched seeded_a to seeded_am2 to get density, not sure if to switch pm2 to starting_pm2
                                 lambdaA + alphaAA + alphaAP ~ warmtrt + (1|time),
                                 nl=TRUE), 
                              data = annuals,
                              family = gaussian, 
                              prior = c(prior(normal(0, 1), nlpar = "lambdaA"), 
                                        prior(normal(0, .1), nlpar = "alphaAA"),
                                        # prior(normal(0, .1), nlpar = "alphaAS"),  #added term for seedling competitive effect
                                        prior(normal(0, .1), nlpar = "alphaAP")),
                              inits="0",
                              cores=4, 
                              chains=4,
                              iter=5000, 
                              thin=1,
                              refresh=100,
                              control = list(adapt_delta = 0.9, max_treedepth = 16))

plot(annual.simple.gaussian)

#look at the plot of relationship between no competition lambda and intraspecific competition
# pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
# pairs(annual.simple.gaussian, pars = c(“lambdaA”, “alphaAA”, “alphaAP”)

saveRDS(annual.simple.gaussian, file="a030122.rds")
#annual.simple.gaussian<-readRDS(file="a020122.rds")
## Get parameters ----
get_variables(annual.simple.gaussian)

## predict and plot COUNTERFACTUALS----

### predict over pm2 (hold seeded_a steady)
dat.new.annual.max <- expand.grid(
  seeded_am2 = max(annuals$seeded_am2, na.rm=TRUE)
  # seeded_a = 0
  # seeded_a = mean(annuals$seeded_a, na.rm=TRUE) #0
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) # 0
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm")
)

dat.new.annual.mean <- expand.grid(
  #seeded_a = max(annuals$seeded_a, na.rm=TRUE)
  # seeded_a = 0
  seeded_am2 = mean(annuals$seeded_am2, na.rm=TRUE) #0
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) # 0
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm")
)

dat.new.annual.0 <- expand.grid(
  seeded_am2 = 0
  # seeded_a = 0
  # seeded_a = mean(annuals$seeded_a, na.rm=TRUE) #0
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) # 0
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm")
)

pred.annual.gaussian.0 <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual.0, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual.0) 
pred.annual.gaussian.mean <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual.mean, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual.mean) 
pred.annual.gaussian.max <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual.max, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual.max) 

a<-ggplot(data=filter(pred.annual.gaussian.0), aes(x = starting_pm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with seeded_a=0")+
  theme_classic()
a1<-ggplot(data=filter(pred.annual.gaussian.0, starting_pm2>0), aes(x = starting_pm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with seeded_a=0")+
  theme_classic()
b<-ggplot(data=filter(pred.annual.gaussian.mean), aes(x = starting_pm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with mean(seeded_a)")+
  theme_classic()
c<-ggplot(data=filter(pred.annual.gaussian.max), aes(x = starting_pm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with max(seeded_a)")+
  theme_classic()


# predict over seeded_a (hold pm2 steady)

dat.new.annual.max2 <- expand.grid(
  seeded_am2 = seq(0,max(annuals$seeded_am2, na.rm=T), length.out=100)
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) #0
  ,starting_pm2 = max(annuals$starting_pm2, na.rm=TRUE) # 0
  ,warmtrt = c("amb","warm")
)
dat.new.annual.mean2 <- expand.grid(
  seeded_am2 = seq(0,max(annuals$seeded_am2, na.rm=T), length.out=100)
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) #0
  ,starting_pm2 = mean(annuals$starting_pm2, na.rm=TRUE) # 0
  ,warmtrt = c("amb","warm")
)
dat.new.annual.02 <- expand.grid(
  seeded_am2 = seq(0,max(annuals$seeded_am2, na.rm=T), length.out=100)
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) #0
  ,starting_pm2 = 0 # 0
  ,warmtrt = c("amb","warm")
)


pred.annual.gaussian.02 <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual.02, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual.02) 
pred.annual.gaussian.mean2 <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual.mean2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual.mean2) 
pred.annual.gaussian.max2 <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual.max2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual.max2) 

f<-ggplot(data=filter(pred.annual.gaussian.max2), aes(x = seeded_am2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with max(pm2)")+
  theme_classic()
e<-ggplot(data=filter(pred.annual.gaussian.mean2), aes(x = seeded_am2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with mean(pm2)")+
  theme_classic()
d<-ggplot(data=filter(pred.annual.gaussian.02), aes(x = seeded_am2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with pm2=0")+
  theme_classic()
d1<-ggplot(data=filter(pred.annual.gaussian.02, seeded_am2>50), aes(x = seeded_am2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=annuals, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Per capita fecundity with pm2=0")+
  theme_classic()

#visualize conditional plots first with all data
grid.arrange(a, b, c, d, e, f, nrow=2, ncol=3)
grid.arrange(a1, b, c, d1, e, f, nrow=2, ncol=3)# then dropping low values in zero conditionals

#rest of annuals - jeff looking at fits. can skip to 279 for perennial adults ----
allfits.annuals <- annual.simple.gaussian %>%
  spread_draws(`b_.*`, regex = TRUE) %>% 
  mutate(
    lam_amb=b_lambdaA_Intercept,
    lam_warm=b_lambdaA_Intercept + b_lambdaA_warmtrtwarm,
    
    alphaAA_amb=b_alphaAA_Intercept,
    alphaAA_warm=b_alphaAA_Intercept + b_alphaAA_warmtrtwarm,
    
    alphaAP_amb=b_alphaAP_Intercept,
    alphaAP_warm=b_alphaAP_Intercept + b_alphaAP_warmtrtwarm,
    
    # alphaAS_amb=b_alphaAS_Intercept,
    # alphaAS_warm=b_alphaAS_Intercept + b_alphaAS_warmtrtwarm
  ) %>%
  dplyr::select(-contains("b_")) %>% 
  pivot_longer(-c(.chain, .iteration, .draw), 
               names_to = c("param", "treatment"),
               names_sep = "_")%>%
  mutate(treatment=ifelse(treatment=="amb", "ambient", "warmed"))

fitsum.annuals<-allfits.annuals%>%
  group_by(param, treatment) %>% 
  # median_qi(.width = c(.95)) 
  mean_qi(.width = c(.66)) 

#look at parameter posterior distributions
ggplot(fitsum.annuals, aes(x=value, y=interaction(treatment, param)))+geom_point(aes(color=treatment))+geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=treatment), width=.2)

ggplot(fitsum.annuals, aes(x=value, y=treatment))+geom_point(aes(color=treatment))+geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=treatment), width=.1) +facet_wrap(~param, scales="free")

#automatically .66 and .95 
ggplot(subset(allfits.annuals, .iteration<10000), aes(x=value, y=treatment))+
  stat_dotsinterval()+facet_wrap(~param, scales="free")+ 
  theme(text=element_text(size=20))


#plot fit lines
bev.annuals <- function(x, lam, AA,AP,N_P) {lam  / (1+AA*x + AP*N_P)} # Beverton-Holt
bev.annuals2 <- function(x, lam, AA,AP,N_A) {lam  / (1+AA*N_A + AP*x)} # Beverton-Holt




################################
### 2. ADULT PERENNIAL FECUNDITY ----
#dat_p is the tibble from data.cleaning here

#visualizations----
p1<-ggplot(dat_p, aes(x=starting_pm2, y=fecundity, color=time)) +
  geom_jitter(aes(shape=time))+
  geom_smooth(aes(linetype=time), method = 'lm',formula = y ~ x , se=F)+
  # scale_colour_manual(values = c("dodgerblue", "darkred"))+
  labs(x="adult perennial density", y="adult perennial fecundity")

p2<-ggplot(dat_p, aes(x=seeded_am2, y=fecundity, color=time)) + # warmtrt
  geom_jitter(aes(shape=time))+
  geom_smooth(aes(linetype=time), method = 'lm',formula = y ~ x, se=F)+
  # scale_colour_manual(values = c("dodgerblue", "darkred"))+
  labs(x="annual seeds in", y="adult perennial fecundity")

ggarrange(p1, p2, common.legend = T)



#Fit BRM----

head(dat_p)
hist(dat_p$fecundity)
hist(dat_p$pm2)
hist(dat_p$seeded_a)

# make scaled version - model runs great if you use these; just complicates interpretation
dat_p$fecundity_scaled <- as.vector(scale(dat_p$fecundity))
dat_p$seeded_a_scaled <- as.vector(scale(dat_p$seeded_a))
dat_p$pm2_scaled <- as.vector(scale(dat_p$pm2))


# find reasonable possible param values by fitting simple nls model
eq1 <- function(x1,x2, lam, pa, pp) {
  lam/(1+ pa*x1 + pp*x2)
}
eq1(2,2, 10, 3,3)  # testing the equation works: 3^2 + 1 = 10

m.nls <- nls(fecundity ~ eq1(seeded_a, pm2, lam,pa,pp), data = dat_p, start = list(lam=20000, pa = 1, pp = 1), trace = F)
summary(m.nls)

# lam=15570
# pa=.0007
# pp=0.09

# graph nls results 
# not working - need to troubleshoot
# seeded.sim <- seq(from = 0, to = 1500, length = 100)
# pm.sim <- seq(from = 0, to = 10, length = 100)
# plot(dat_p$seeded_a, dat_p$fecundity)
# pred1 <- predict(m.nls, list(x1 = seeded.sim, x2=pm.sim))
# lines(seeded.sim, pred1, col = "green")


# alt way to graph
# check out reasonable possible param values
ggplot(dat_p, aes(x= seeded_am2, y=fecundity)) +
  geom_point() +
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x), # BH
              method.args = list(start = list(a = 20000, b = .1)),
              se = FALSE, 
              color='orange')

# check out reasonable possible param values
ggplot(dat_p, aes(x= pm2, y=fecundity)) +
  geom_point() +
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x), # BH
              method.args = list(start = list(a = 20000, b = .1)),
              se = FALSE, 
              color='orange')

### 2.1 Perennial fecundity scaled
# new SIMPLE

#dat_p tibble is the one to use for modeling
#each row is for an adult phytometer
#uses mostly the same column names as above
#type is an important column - previously we have only run this with type="adult".  
# the other types present are adult_g (obviously damaged by gophers) and newadult (a seedling that has become an adult)


adult.simple.gaussian <- brm(bf(fecundity/20000 ~ lambdaP / (1+alphaPA*seeded_am2 + alphaPP*starting_pm2), # lambdaP*20000
  lambdaP + alphaPA + alphaPP ~ warmtrt + (1|time), 
    nl=TRUE), 
    data = dat_p,
    family = gaussian, 
    prior = c(prior(normal(0, 1),nlpar = "lambdaP"), 
    prior(normal(.001, .01), nlpar = "alphaPA"),
    prior(normal(.1, .1), nlpar = "alphaPP")),
    inits = "0",  
    cores=3, 
    chains=3,
    iter=5000, 
    thin=2,
     control = list(adapt_delta = 0.9, max_treedepth = 16),
    refresh=100,
)

conditional_effects(adult.simple.gaussian)
adult.simple.gaussian
plot(adult.simple.gaussian)
#saveRDS(adult.simple.gaussian, file="adult_simple_JD.rds")
saveRDS(adult.simple.gaussian, file="p030122.rds")
#readRDS(adult.simple.gaussian, file="p020122.rds")



# Plot model results - JD - need to remember how to do this, and how to back-convert the parameters and predictions given the transformations in the model





## predict and plot COUNTERFACTUALS----

### predict over pm2 (hold seeded_a steady)
dat.new.adult.max <- expand.grid(
  seeded_am2 = max(dat_p$seeded_am2, na.rm=TRUE)
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm"))
dat.new.adult.mean <- expand.grid(
  seeded_am2 = mean(dat_p$seeded_am2, na.rm=TRUE) #0
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm"))
dat.new.adult.0 <- expand.grid(
  seeded_am2 = 0
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm"))

pred.adult.gaussian.0 <- 
  as.data.frame(predict(adult.simple.gaussian, newdata = dat.new.adult.0, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.adult.0) 
pred.adult.gaussian.mean <- 
  as.data.frame(predict(adult.simple.gaussian, newdata = dat.new.adult.mean, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.adult.mean) 
pred.adult.gaussian.max <- 
  as.data.frame(predict(adult.simple.gaussian, newdata = dat.new.adult.max, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.adult.max) 

#add in estimate*1000 to get scaling corrected.
a<-ggplot(data=filter(pred.adult.gaussian.0), aes(x = starting_pm2, y = Estimate*20000,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5*20000, ymax=Q95*20000, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=dat_p, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult perennial fecundity with seeded_a=0")+
  theme_classic()
b<-ggplot(data=filter(pred.adult.gaussian.mean), aes(x = starting_pm2, y = Estimate*20000,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5*20000, ymax=Q95*20000, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=dat_p, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult perennial fecundity with mean(seeded_a)")+
  theme_classic()
c<-ggplot(data=filter(pred.adult.gaussian.max), aes(x = starting_pm2, y = Estimate*20000,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5*20000, ymax=Q95*20000, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=dat_p, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult perennial fecundity with max(seeded_a)")+
  theme_classic()


# predict over seeded_a (hold pm2 steady)

dat.new.adult.max2 <- expand.grid(
  seeded_am2 = seq(0,max(dat_p$seeded_am2, na.rm=TRUE), length.out=100)
  ,starting_pm2 = max(dat_p$starting_pm2, na.rm=TRUE) # 0
  ,warmtrt = c("amb","warm")
)
dat.new.adult.mean2 <- expand.grid(
  seeded_am2 = seq(0,max(dat_p$seeded_am2, na.rm=TRUE), length.out=100)
  ,starting_pm2 = mean(dat_p$starting_pm2, na.rm=TRUE) # 0
  ,warmtrt = c("amb","warm")
)
dat.new.adult.02 <- expand.grid(
  seeded_am2 = seq(0,max(dat_p$seeded_am2, na.rm=TRUE), length.out=100)
  ,starting_pm2 = 0
  ,warmtrt = c("amb","warm")
)


pred.adult.gaussian.02 <- 
  as.data.frame(predict(adult.simple.gaussian, newdata = dat.new.adult.02, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.adult.02) 
pred.adult.gaussian.mean2 <- 
  as.data.frame(predict(adult.simple.gaussian, newdata = dat.new.adult.mean2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.adult.mean2) 
pred.adult.gaussian.max2 <- 
  as.data.frame(predict(adult.simple.gaussian, newdata = dat.new.adult.max2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.adult.max2) 

f<-ggplot(data=filter(pred.adult.gaussian.max2), aes(x = seeded_am2, y = Estimate*20000,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5*20000, ymax=Q95*20000, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=dat_p, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult perennial fecundity with max(pm2)")+
  theme_classic()
e<-ggplot(data=filter(pred.adult.gaussian.mean2), aes(x = seeded_am2, y = Estimate*20000,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5*20000, ymax=Q95*20000, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=dat_p, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult perennial fecundity with mean(pm2)")+
  theme_classic()
d<-ggplot(data=filter(pred.adult.gaussian.02), aes(x = seeded_am2, y = Estimate*20000,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5*20000, ymax=Q95*20000, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=dat_p, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult perennial fecundity with pm2=0")+
  theme_classic()
grid.arrange(a, b, c, d, e, f, nrow=2, ncol=3)

## Get parameters: jeff's work. can skip to 610 for perennial seedlings ----
get_variables(adult.simple.gaussian)

allfits.adult <- adult.simple.gaussian %>%
  spread_draws(`b_.*`, regex = TRUE) %>% 
  mutate(
    lam_amb=b_lambdaP_Intercept,
    lam_warm=b_lambdaP_Intercept + b_lambdaP_warmtrtwarm,
    
    alphaPA_amb=b_alphaPA_Intercept,
    alphaPA_warm=b_alphaPA_Intercept + b_alphaPA_warmtrtwarm,
    
    alphaPP_amb=b_alphaPP_Intercept,
    alphaPP_warm=b_alphaPP_Intercept + b_alphaPP_warmtrtwarm,
    
    # alphaAS_amb=b_alphaAS_Intercept,
    # alphaAS_warm=b_alphaAS_Intercept + b_alphaAS_warmtrtwarm
  ) %>%
  dplyr::select(-contains("b_")) %>% 
  pivot_longer(-c(.chain, .iteration, .draw), 
               names_to = c("param", "treatment"),
               names_sep = "_")%>%
  mutate(treatment=ifelse(treatment=="amb", "ambient", "warmed"))

fitsum.adult<-allfits.adult%>%
  group_by(param, treatment) %>% 
  # median_qi(.width = c(.95)) 
  mean_qi(.width = c(.66)) 

#look at parameter posterior distributions
ggplot(fitsum.adult, aes(x=value, y=interaction(treatment, param)))+geom_point(aes(color=treatment))+geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=treatment), width=.2)

ggplot(fitsum.adult, aes(x=value, y=treatment))+geom_point(aes(color=treatment))+geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=treatment), width=.1) +facet_wrap(~param, scales="free")

#automatically .66 and .95 
ggplot(subset(allfits.adult, .iteration>5000), aes(x=value, y=treatment))+
  stat_dotsinterval()+facet_wrap(~param, scales="free")+ 
  theme(text=element_text(size=20))


#plot fit lines
bev.adult.self <- function(x, lam, PA,PP,N_A) {lam  / (1+PP*x + PA*N_A)} # Beverton-Holt
bev.adult.inter <- function(x, lam, PA,PP,N_P) {lam  / (1+PP*N_P + PA*x)} # Beverton-Holt


#plot annual fecundity when other competitor is at mean value
ggplot(subset(datp, seeded_a<500), aes(x = density_p, y = out_p, color=warmtrt)) + 
  geom_jitter()+
#  geom_function(fun=bev.adult.self, color='dodgerblue', size=1.5,
#                args=c("lam"=5000*fitsum.adult$value[5], "PP"=fitsum.adult$value[3], 
#                       # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
#                       "PA"=fitsum.annuals$value[1], "N_A"=mean(datp$seeded_a)
#                ))+
#  geom_function(fun=bev.adult.self, color='darkred', size=1.5,
#                args=c("lam"=5000*fitsum.adult$value[6], "PP"=fitsum.adult$value[4], 
#                       # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
#                       "PA"=fitsum.adult$value[2], "N_A"=mean(datp$seeded_a)
#                ))+
  geom_smooth(method=lm, se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+
  ylab("Per capita adult perennial fecundity")+
  xlab("Adult perennials/m2")+ 
  theme(text=element_text(size=16))

ggplot(subset(datp, density_p<3), aes(x = seeded_a, y = out_p, color=warmtrt)) + 
  geom_jitter()+
#  geom_function(fun=bev.adult.inter, color='dodgerblue', size=1.5,
##                args=c("lam"=5000*fitsum.adult$value[5], "PP"=fitsum.adult$value[3], 
#                       # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
#                       "PA"=fitsum.annuals$value[1], "N_P"=mean(datp$density_p)
#                ))+
#  geom_function(fun=bev.adult.inter, color='darkred', size=1.5,
#                args=c("lam"=5000*fitsum.adult$value[6], "PP"=fitsum.adult$value[4], 
#                       # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
#                       "PA"=fitsum.adult$value[2], "N_P"=mean(datp$density_p)
#                ))+
  geom_smooth(se=F, method=lm)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+
  ylab("Per capita adult perennial fecundity")+
  xlab("Seeded annuals/m2")+ 
  theme(text=element_text(size=16))


ggplot(subset(datp, comptrt=="none"), aes(x=warmtrt, y=out_p)) +
  geom_boxplot() +
  geom_jitter(width=.2)+xlab("")+
  ylab("Perennial percapita fecundity no competition")+ 
  theme(text=element_text(size=16))




#old ambient/warm ----
perennial.lambda <- brm(bf(out_p+1 ~ lambdaP*5000 / (1+alphaPA*seeded_a + alphaPP*density_p)+1, # increased scale to 5000
                           lambdaP ~ 1 + (1|block),
                           alphaPA ~  1 + (1|block),# didnt add seedling competition on adults, dont think it would help here
                           alphaPP ~  1 + (1|block),
                           nl=TRUE),
                        family=poisson, # switched to poisson
                        data = subset(datp, warmtrt=="amb"),
                        prior = c(prior(normal(1, .5), lb=0, nlpar = "lambdaP"), 
                                  prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaPA"), # added bounds to the priors
                                  prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaPP")),
                        # inits = "0",  
                        cores=4, 
                        chains=4,
                        refresh=100,
                        iter=10000, 
                        thin=3
                        ,control = list(adapt_delta = 0.99, max_treedepth = 16)
)

perennial.lambda <- brm(bf(out_p+1 ~ lambdaP*5000 / (1+alphaPA*seeded_a + alphaPP*density_p)+1, # increased scale to 5000
                           lambdaP ~ 1 + (1|block),
                           alphaPA ~  1 + (1|block),# didnt add seedling competition on adults, dont think it would help here
                           alphaPP ~  1 + (1|block),
                           nl=TRUE),
                        family=poisson, # switched to poisson
                        data = subset(datp, warmtrt=="warm"),
                        prior = c(prior(normal(1, .5), lb=0, nlpar = "lambdaP"), 
                                  prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaPA"), # added bounds to the priors
                                  prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaPP")),
                        # inits = "0",  
                        cores=4, 
                        chains=4,
                        refresh=100,
                        iter=10000, 
                        thin=3
                        ,control = list(adapt_delta = 0.99, max_treedepth = 16)
)
#plw1<-perennial.lambda
#stancode(perennial.lambda)
#perennial.lambda
#plot(perennial.lambda)
#fixef(perennial.lambda)
#conditional_effects(readPL)%>% plot(points=T)
conditional_effects(readPL, effects = "seeded_a:warmtrt")%>% plot(points=T)
conditional_effects(readPL, effects = "density_p:warmtrt")%>% plot(points=T)

#the closes to a decent fit i've gotten so far, can't replicate it when I do the annual and warmed seperately. 
readPL<-readRDS("PL.rds") 


saveRDS(perennial.lambda.w, file="PL070121w.rds")
saveRDS(perennial.lambda.a, file="PL070121a.rds")

adult.warm.model<-readRDS("PL070121w.rds")
adult.amb.model<-readRDS("PL070121a.rds") #this is ambient fitted alone, way lower ESS?


### 2.2 Perennial fecundity logscaled ----
perennial.lambda.logscale <- brm(bf(log(1+out_p) ~ log(exp(log(lambdaP) - log((1+alphaPA*(seeded_a+1) + alphaPP*(density_p+1))))),
                           lambdaP ~ warmtrt + (1|block),
                           alphaPA ~  warmtrt + (1|block),# didnt add seedling competition on adults, dont think it would help
                           alphaPP ~  warmtrt + (1|block),
                           nl=TRUE),
                        data = subset(datp),
                        prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaP"), 
                                  prior(normal(0, .2), nlpar = "alphaPA"),
                                  prior(normal(0, .2), nlpar = "alphaPP")),
                        inits = "0",  
                        cores=4, 
                        chains=4,
                        refresh=100,
                        iter=20000, 
                        # thin=2
                        #  ,control = list(adapt_delta = 0.99, max_treedepth = 17)
)

perennial.lambda.logscale
plot(perennial.lambda.logscale)

conditional_effects(perennial.lambda.logscale)%>% plot(points=T)
conditional_effects(perennial.lambda.logscale, effects = "seeded_a:warmtrt")%>% plot(points=T)
conditional_effects(perennial.lambda.logscale, effects = "density_p:warmtrt")%>% plot(points=T)

savedPLL<-perennial.lambda.logscale
#saveRDS(savedPLL, file="PLL.rds")
#readPLL<-readRDS("PLL.rds")

#############################################
### 3. PERENNIAL SEEDLINGS  ----

#vizualisation ----

#seeds:adults vs. seeds
s1<-ggplot(seedlings, aes(x=seeded_sm2, y=fall.g/seeded_sm2, color=warmtrt)) + #switch the gopehr effect to spring not fall!!!
  geom_point(aes(shape=time))+ geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("seeded perennials")+ylab("seed:adult survival rate")

#seeds:adults vs. perennials
s2<-ggplot(seedlings, aes(x=starting_pm2, y=fall.g/seeded_sm2, color=warmtrt)) +
  geom_point(aes(shape=time))+ geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("adult perennial density")+ylab("seed:adult survival rate")

#seeds:adults vs. annuals
s3<-ggplot(seedlings, aes(x=seeded_am2, y=fall.g/seeded_sm2, color=warmtrt)) +
  geom_point(aes(shape=time))+ geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("seeded annuals")+ylab("seed:adult survival rate")

ggarrange(s1, s2, s3, common.legend = T, nrow=1, ncol=3)

# model: fall.g/seeded_s ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*pm2), 

### 3.1 Simple (seeds in:adults out) perennial seedlings binomial ----
#seedling.binomial<- brm(bf(as.integer(fall.g)|trials(as.integer(seeded_s)) ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*pm2), 
#                           lambdaS +alphaSA +alphaSP+alphaSS~ warmtrt + (1|time), 
#                          nl=TRUE),
#                        family=binomial,
#                        data = seedlings,  
#                        prior = c(prior(normal(0.05, .05),nlpar = "lambdaS"), 
##                                  prior(normal(0, .1), nlpar = "alphaSA"),
#                                prior(normal(0, .1), nlpar = "alphaSS"),
#                                  prior(normal(0, .1), nlpar = "alphaSP")),
#                     #   inits = "0",  
#                        cores=4, 
#                        chains=4,
#                        iter=5000, 
#                        thin=1,
#                        control = list(adapt_delta = 0.95, max_treedepth = 16))


#Gaussian alternative 
fall.g/seeded_s

# seedling.simple.gaussian<- brm(bf(fall.g/seeded_s ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*pm2), 
                                  
# eq2 <- function(x1,x2,x3, lam, sa, ss, sp) {
#   lam/(1+ sa*x1 + ss*x2 + sp*x3)
# }
# m2.nls <- nls(fall.g/seeded_s ~ eq2(seeded_a, seeded_s, pm2, lam,sa, ss, sp), data = seedlings, start = list(lam=.03, sa = 1, ss = 1, sp=1), trace = F)
# summary(m2.nls)

#seedlings is the one to use for modeling
# seeded_a refers to the actual number of annual seeds added, regardless of gopher damage
#percap is the number seeds out/ seeds in
#seeded_am2 refers to the number of seeds added/m2# pm2 refers to the m2 density of ACTUAL adult perrennials in the plot
#starting_pm2 refers to the density/m2 of adult perennials that should be in the plot ignoring gophers
# dpm2 refers to the number of damaged adult perennials per meter square
# nam2 refers to the number of 'new adults' per meter square, seedlings that were left as phytometers to see in the next year
#now the model runs off of pm2, but it could be adapted for dpm2, nam2, or starting pm2 - see how their fecundities vary below:


seedling.simple.gaussian<- brm(bf(fall.g/seeded_sm2 ~ lambdaS / (1+alphaSA*seeded_am2 + alphaSS*seeded_sm2 + alphaSP*starting_pm2), #update data terms with new from above
                                  lambdaS +alphaSA +alphaSP+alphaSS~ warmtrt + (1|time), 
                                  nl=TRUE),
                      family=gaussian,
                      data = seedlings,   #running this with limited dataset as in the figures above (only in seedling comptrts)
                      prior = c(prior(normal(.05, .1), nlpar = "lambdaS"), 
                                prior(normal(0.1, .1),    nlpar = "alphaSA"),
                                prior(normal(0.1, .1),    nlpar = "alphaSS"),
                                prior(normal(0.1, .1),    nlpar = "alphaSP")),
                         inits = "0",  
                      cores=3, 
                      chains=3,
                      iter=5000, 
                      thin=2,
                      control = list(adapt_delta = 0.9, max_treedepth = 16))

plot(seedling.simple.gaussian)
seedling.simple.gaussian
saveRDS(seedling.simple.gaussian, file="s030122.rds")
#readRDS(seedling.simple.gaussian, file="s020122.rds")
#try counterfactual plots anyway
## predict and plot COUNTERFACTUALS----

### predict over pm2 (hold seeded_a and seeded_s steady)
dat.new.seedling.max <- expand.grid(
  seeded_am2 = max(seedlings$seeded_am2, na.rm=TRUE),
  seeded_sm2 = as.integer(max(seedlings$seeded_sm2, na.rm=T))
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm"))
dat.new.seedling.mean <- expand.grid(
  seeded_am2 = mean(seedlings$seeded_am2, na.rm=TRUE),
  seeded_sm2 = as.integer(mean(seedlings$seeded_sm2, na.rm=T))
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm"))
dat.new.seedling.0 <- expand.grid(
  seeded_am2 = 0,
  seeded_sm2 = as.integer(0)
  ,starting_pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm"))

pred.seedling.gaussian.0 <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.0, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.0) 
pred.seedling.gaussian.mean <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.mean, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.mean) 
pred.seedling.gaussian.max <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.max, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.max) 

a<-ggplot(data=filter(pred.seedling.gaussian.0), aes(x = starting_pm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=starting_pm2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with seeded_a&seeded_s=0")+
  theme_classic()
b<-ggplot(data=filter(pred.seedling.gaussian.mean), aes(x = starting_pm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=starting_pm2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with mean(seeded_a&seeded_s)")+
  theme_classic()
c<-ggplot(data=filter(pred.seedling.gaussian.max), aes(x = starting_pm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=starting_pm2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with max(seeded_a&seeded_s)")+
  theme_classic()


# predict over seeded_a (hold density_p and seeded_s steady)

dat.new.seedling.max2 <- expand.grid(
  seeded_am2 = seq(0,max(seedlings$seeded_am2), length.out=100)
  ,starting_pm2 = max(seedlings$starting_pm2, na.rm=TRUE) # 0
  ,seeded_sm2 = max(seedlings$seeded_sm2, na.rm=TRUE)
  ,warmtrt = c("amb","warm")
)
dat.new.seedling.mean2 <- expand.grid(
  seeded_am2 = seq(0,max(seedlings$seeded_am2), length.out=100)
  ,starting_pm2 = mean(seedlings$starting_pm2, na.rm=TRUE) # 0
  ,seeded_sm2 = mean(seedlings$seeded_sm2, na.rm=TRUE)
  ,warmtrt = c("amb","warm")
)
dat.new.seedling.02 <- expand.grid(
  seeded_am2 = seq(0,max(seedlings$seeded_am2), length.out=100)
  ,starting_pm2 = 0
  ,seeded_sm2 = 0
  ,warmtrt = c("amb","warm")
)

pred.seedling.gaussian.02 <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.02, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.02) 
pred.seedling.gaussian.mean2 <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.mean2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.mean2) 
pred.seedling.gaussian.max2 <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.max2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.max2) 

f<-ggplot(data=filter(pred.seedling.gaussian.max2), aes(x = seeded_am2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=seeded_am2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with max(density_p&seeded_s)")+
  theme_classic()
e<-ggplot(data=filter(pred.seedling.gaussian.mean2), aes(x = seeded_am2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=seeded_am2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with mean(density_p&seeded_s)")+
  theme_classic()
d<-ggplot(data=filter(pred.seedling.gaussian.02), aes(x = seeded_am2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=seeded_am2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with density_p&seeded_s=0")+
  theme_classic()

# predict over seeded_s (hold density_p and seeded_a steady)

dat.new.seedling.max2 <- expand.grid(
  seeded_sm2 = seq(0,max(seedlings$seeded_sm2), length.out=100)
  ,starting_pm2 = max(seedlings$starting_pm2, na.rm=TRUE) # 0
  ,seeded_am2 = max(seedlings$seeded_am2, na.rm=TRUE)
  ,warmtrt = c("amb","warm")
)
dat.new.seedling.mean2 <- expand.grid(
  seeded_sm2 = seq(0,max(seedlings$seeded_sm2), length.out=100)
  ,starting_pm2 = mean(seedlings$starting_pm2, na.rm=TRUE) # 0
  ,seeded_am2 = mean(seedlings$seeded_am2, na.rm=TRUE)
  ,warmtrt = c("amb","warm")
)
dat.new.seedling.02 <- expand.grid(
  seeded_sm2 = seq(0,max(seedlings$seeded_sm2), length.out=100)
  ,starting_pm2 = 0
  ,seeded_am2 = 0
  ,warmtrt = c("amb","warm")
)

pred.seedling.gaussian.02 <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.02, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.02) 
pred.seedling.gaussian.mean2 <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.mean2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.mean2) 
pred.seedling.gaussian.max2 <- 
  as.data.frame(predict(seedling.simple.gaussian, newdata = dat.new.seedling.max2, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling.max2) 

i<-ggplot(data=filter(pred.seedling.gaussian.max2), aes(x = seeded_sm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=seeded_sm2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with max(density_p&seeded_a)")+
  theme_classic()
h<-ggplot(data=filter(pred.seedling.gaussian.mean2), aes(x = seeded_sm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=seeded_sm2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with mean(density_p&seeded_a)")+
  theme_classic()
g<-ggplot(data=filter(pred.seedling.gaussian.02), aes(x = seeded_sm2, y = Estimate,  color=warmtrt)) + 
  geom_smooth(aes(ymin=Q5, ymax=Q95, fill=warmtrt), stat="identity", alpha = 1/5, size = 1/4) +
  geom_point()+
  geom_jitter(data=seedlings, aes(x=seeded_sm2, y=fall.g/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival with density_p&seeded_a=0")+
  theme_classic()

ggarrange(a, b, c, d, e, f, g, h, i, nrow=3, ncol=3)

### STOP HERE ###

##alternates----
#try splitting warmtrt and add back bounds
seedling.simple.gaussian.amb.b<- brm(bf(fall20_s/seeded_s.g ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                                  lambdaS+alphaSA+alphaSP+alphaSS ~ 1, #1+ (1|block), 
                                  nl=TRUE),
                               family=gaussian,
                               data = subset(dat.sumsurv, warmtrt=="amb"),   #running this with limited dataset as in teh figures above (only in seedling comptrts)
                               prior = c(prior(normal(.03, .03), nlpar = "lambdaS"), 
                                         prior(normal(0, .01),    nlpar = "alphaSA"),
                                         prior(normal(0, .01),   nlpar = "alphaSS"),
                                         prior(normal(0, .01),     nlpar = "alphaSP")),
                               #   inits = "0",  
                               cores=4, 
                               chains=4,
                               iter=15000, 
                               thin=5,
                               control = list(adapt_delta = 0.99, max_treedepth = 18))

seedling.simple.gaussian.warm.b<- brm(bf(fall20_s/seeded_s.g ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                                      lambdaS+alphaSA+alphaSP+alphaSS ~ 1, #1+ (1|block), 
                                      nl=TRUE),
                                   family=gaussian,
                                   data = subset(dat.sumsurv, warmtrt=="warm"),   #running this with limited dataset as in teh figures above (only in seedling comptrts)
                                   prior = c(prior(normal(.03, .03),  nlpar = "lambdaS"), 
                                             prior(normal(0, .01),    nlpar = "alphaSA"),
                                             prior(normal(0, .01),   nlpar = "alphaSS"),
                                             prior(normal(0, .01),    nlpar = "alphaSP")),
                                   #   inits = "0",  
                                   cores=4, 
                                   chains=4,
                                   iter=15000, 
                                   thin=5,
                                   control = list(adapt_delta = 0.99, max_treedepth = 18))


## Get parameters ----
get_variables(seedling.simple.gaussian.amb)
get_variables(seedling.simple.gaussian.amb)


allfits.seedling.amb <- seedling.simple %>%
  spread_draws(`b_.*`, regex = TRUE) 
allfits.seedling.warm <- seedling.simple %>%
  spread_draws(`b_.*`, regex = TRUE) 
allfits.seedling<-allfits.seedling.amb%>% 
  mutate(
    lam_amb=b_lambdaS_Intercept,
    lam_warm=b_lambdaS_Intercept + allfits.seedling.warm$b_lambdaS_Intercept,
    
    alphaSA_amb=b_alphaSA_Intercept,
    alphaSA_warm=b_alphaSA_Intercept + allfits.seedling.warm$b_alphaSA_Intercept,
    
    alphaSP_amb=b_alphaSP_Intercept,
    alphaSP_warm=b_alphaSP_Intercept + allfits.seedling.warm$b_alphaSP_Intercept,
    
     alphaSS_amb=b_alphaSS_Intercept,
     alphaSS_warm=b_alphaSS_Intercept + allfits.seedling.warm$b_alphaSS_Intercept
  ) %>%
  dplyr::select(-contains("b_")) %>% 
  pivot_longer(-c(.chain, .iteration, .draw), 
               names_to = c("param", "treatment"),
               names_sep = "_")%>%
  mutate(treatment=ifelse(treatment=="amb", "ambient", "warmed"))

fitsum.seedling<-allfits.seedling%>%
  group_by(param, treatment) %>% 
  # median_qi(.width = c(.95)) 
  mean_qi(.width = c(.66)) 

#look at parameter posterior distributions
ggplot(fitsum.seedling, aes(x=value, y=interaction(treatment, param)))+geom_point(aes(color=treatment))+geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=treatment), width=.2)

ggplot(fitsum.seedling, aes(x=value, y=treatment))+geom_point(aes(color=treatment))+geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=treatment), width=.1) +facet_wrap(~param, scales="free")

#automatically .66 and .95 
ggplot(subset(allfits.seedling), aes(x=value, y=treatment))+
  stat_dotsinterval()+facet_wrap(~param, scales="free")+ 
  theme(text=element_text(size=20))



#plot fit lines
bev.seedling.self <-   function(x, lam, SA,SP,SS, N_A, N_P) {lam  / (1+SS*x + SP*N_P + SA*N_A)} # Beverton-Holt
bev.seedling.interp <- function(x, lam, SA,SP,SS, N_S, N_A) {lam  / (1+SS*N_S + SP*x + SA*N_A)} # Beverton-Holt
bev.seedling.intera <- function(x, lam, SA,SP,SS, N_S, N_P) {lam  / (1+SP*N_S + SP*N_P+ SA*x)} # Beverton-Holt


#plot seedling survival vary (why is the prediction so small??)
ggplot(subset(dat.sumsurv, seeded_s>100&seeded_a<5000&density_p<3), aes(x=seeded_s.g, y=fall20_s/seeded_s.g, color=warmtrt))+
  geom_jitter()+
#  geom_function(fun=bev.seedling.self, color='dodgerblue', size=1.5,
#                args=c("lam"=fitsum.seedling$value[7], "SS"=fitsum.seedling$value[5], 
#                        "SA"=fitsum.seedling$value[1], "N_A"=mean(dat.sumsurv$seeded_a),
#                       "SP"=fitsum.seedling$value[3], "N_P"=mean(dat.sumsurv$density_p)
#                ))+
#  geom_function(fun=bev.seedling.self, color='darkred', size=1.5,
##                args=c("lam"=fitsum.seedling$value[8], "SS"=fitsum.seedling$value[6], 
#                       "SA"=fitsum.seedling$value[2], "N_A"=mean(dat.sumsurv$seeded_a),
#                       "SP"=fitsum.seedling$value[4], "N_P"=mean(dat.sumsurv$density_p)
#                ))+
  geom_smooth(method=lm, se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+
  ylab("Proportion of seeds surviving to adulthood")+
  xlab("Seeded perennials/m2")+ 
  theme(text=element_text(size=16))

#vary adults
ggplot(subset(dat.sumsurv, seeded_s>100&seeded_a<5000), aes(x=density_p, y=fall20_s/seeded_s.g, color=warmtrt))+
  geom_jitter()+
#  geom_function(fun=bev.seedling.interp, color='dodgerblue', size=1.5,
#                args=c("lam"=fitsum.seedling$value[7], "SS"=fitsum.seedling$value[5], 
#                       "SA"=fitsum.seedling$value[1], "N_A"=mean(dat.sumsurv$seeded_a),
#                       "SP"=fitsum.seedling$value[3], "N_S"=mean(dat.sumsurv$seeded_s)
#                ))+
#  geom_function(fun=bev.seedling.interp, color='darkred', size=1.5,
#                args=c("lam"=fitsum.seedling$value[8], "SS"=fitsum.seedling$value[6], 
#                       "SA"=fitsum.seedling$value[2], "N_A"=mean(dat.sumsurv$seeded_a),
#                       "SP"=fitsum.seedling$value[4], "N_S"=mean(dat.sumsurv$seeded_s)
#                ))+
  geom_smooth(method=lm, se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+
  ylab("")+
  xlab("Adult perennial density/m2")+ 
  theme(text=element_text(size=16))

ggplot(subset(dat.sumsurv, comptrt=="none"), aes(x=warmtrt, y=fall20_s/seeded_s.g)) +
  geom_boxplot()+geom_jitter(width=.2)+
  ylab("perennial seedling survival no competition")+ 
  theme(text=element_text(size=16))


#vary annuals
ggplot(subset(dat.sumsurv, seeded_s>100&density_p<3), aes(x=seeded_a, y=fall20_s/seeded_s.g, color=warmtrt))+
  geom_jitter()+
#  geom_function(fun=bev.seedling.intera, color='dodgerblue', size=1.5,
#                args=c("lam"=fitsum.seedling$value[7], "SS"=fitsum.seedling$value[5], 
#                       "SA"=fitsum.seedling$value[1], "N_S"=mean(dat.sumsurv$seeded_s),
#                       "SP"=fitsum.seedling$value[3], "N_P"=mean(dat.sumsurv$density_p)
#                ))+
#  geom_function(fun=bev.seedling.intera, color='darkred', size=1.5,
#                args=c("lam"=fitsum.seedling$value[8], "SS"=fitsum.seedling$value[6], 
#                       "SA"=fitsum.seedling$value[2], "N_S"=mean(dat.sumsurv$seeded_s),
#                       "SP"=fitsum.seedling$value[4], "N_P"=mean(dat.sumsurv$density_p)
#                ))+
geom_smooth(method=lm, se=F)+ 
 scale_colour_manual(values = c("dodgerblue", "darkred"))+
  ylab("")+
  xlab("Seeded annuals/m2")+ 
  theme(text=element_text(size=16))


# try to predict seedling as a factor of perennial density using predict() rather than bev holt in case of transformation issues
  dat.new.seedling <- expand.grid(
    #seeded_a = max(annuals$seeded_a, na.rm=TRUE)
    # seeded_a = 0
    seeded_a = mean(dat.sumsurv$seeded_a, na.rm=TRUE) #0
    ,seeded_s.g= as.integer(mean(dat.sumsurv$seeded_s.g, na.rm=TRUE)) # 0
    ,density_p = seq(0,10, length.out=20)
    ,warmtrt = c("amb","warm")
  )

# try predicting seedling lines
pred.seedling <- 
  as.data.frame(predict(seedling.simple, newdata = dat.new.seedling, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.seedling) 

ggplot(pred.seedling, aes(x = density_p, y = Estimate, color=warmtrt)) +  #why is the estimate so big? 
  geom_line() + 
  geom_jitter(data=dat.sumsurv, aes(x=density_p, y=fall20_s/(seeded_s.g)))+
  ylab("seedling survival")



#separate warmed and amb ----
seedling.binomial.w<- brm(bf(fall20_s|trials(seeded_s.g) ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                             lambdaS +alphaSA +alphaSP+alphaSS~ 1+ (1|block), nl=TRUE),
                          family=binomial,
                          data = subset(dat.sumsurv, seeded_s>150&warmtrt=="warm"),   #running this with limited dataset as in teh figures above (only in seedling comptrts)
                          prior = c(prior(normal(.015, .015), lb=0, nlpar = "lambdaS"), 
                                    prior(normal(0, .1), lb=0, ub=1, nlpar = "alphaSA"),
                                    prior(normal(0, .1), lb=0, ub=1, nlpar = "alphaSS"),
                                    prior(normal(0, .1), lb=0, ub=1, nlpar = "alphaSP")),
                          inits = "0",  
                          cores=4, 
                          chains=4,
                          iter=15000, 
                          thin=5,
                          control = list(adapt_delta = 0.9, max_treedepth = 18))
)

savedps2<-readRDS("PS.rds") #previous best, model w/ warmed and ambient

saveRDS(seedling.binomial.a, file="PS070321a")
saveRDS(seedling.binomial.w, file="PS070321w")

seedling.a<-readRDS("PS070321a")
seedling.w<-readRDS("PS070321w")



 # conditional_effects(savedps2)%>% plot(points=T) # I am getting errors when I try to see conditional effects of the binomial model.  It says that it is exporting 
#  conditional_effects(savedPS, conditions = data.frame(trials = 100), effects = "seeded_a:warmtrt")%>% plot(points=T)
#  conditional_effects(savedPS, effects = "density_p:warmtrt")%>% plot(points=T)
#  plot(conditional_effects(savedPS, conditions = data.frame(trials = 10)), points = T)
  #saveRDS(savedPS, file="PS.rds")
 # savedps2<-readRDS("PS.rds")
  
### 3.2 Spring (seeds in:stems out) perennial seedlings binomial (first half of 3.1) (SKIP - DIDNT SEEM LIKE 3.2 or 3.3 added much)----
  sprsur.binomial<- brm(bf(spring_20s|trials(seeded_s) ~ sprsurS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                             sprsurS ~ warmtrt+ (1|block), 
                             alphaSA ~  warmtrt+ (1|block), 
                             alphaSP ~  warmtrt+ (1|block), 
                             alphaSS ~  warmtrt+ (1|block), nl=TRUE),
                          family=binomial,
                          data = subset(dat.sprsurv, seeded_s>150),   #running this with limited dataset as in teh figures above (only in seedling comptrts)
                          prior = c(prior(normal(0, .5), lb=0, nlpar = "sprsurS"), 
                                    prior(normal(0, .1), nlpar = "alphaSA"),
                                    prior(normal(0, .1), nlpar = "alphaSS"),
                                    prior(normal(0, .1), nlpar = "alphaSP")),
                          inits = "0",  
                          cores=4, 
                          chains=4,
                          iter=10000, 
                          thin=5,
                          control = list(adapt_delta = 0.9, max_treedepth = 15))
  )
  
  sprsur
  summary(sprsur)
  plot(sprsur)
  fixef(sprsur)
  conditional_effects(sprsur)
  
 # saveRDS(sprsur.binomial, file="SPR.rds")
  #savedspr<-readRDS("SPR.rds")

### 3.3 Summer (stems in:adults out) perennial seedlings binomial (second half of 3.1) ----
sumsur <- brm(bf(fall20_s|trials(spring20_s.g) ~ sumsurS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                            sumsurS ~ warmtrt + (1|block), 
                            alphaSA ~  warmtrt + (1|block), 
                            alphaSP ~  warmtrt + (1|block), 
                            alphaSS ~  warmtrt + (1|block), nl=TRUE),
                         family=binomial,
                         data = subset(dat.sumsurv, seeded_s>150),
                         prior = c(prior(normal(1, 1), lb=0, nlpar = "sumsurS"), 
                                   prior(normal(0, .1), nlpar = "alphaSA"),
                                   prior(normal(0, .1), nlpar = "alphaSS"),
                                   prior(normal(0, .1), nlpar = "alphaSP")),
                         inits = "0",  
                         cores=4, 
                         chains=4,
                         iter=5000, 
                         thin=5,
                         control = list(adapt_delta = 0.99, max_treedepth = 18))
  
  plot(seedling_sumsur)
  summary(seedling_sumsur)
  
  saveRDS(sumsur, file="SUM.rds")
  savedsum<-readRDS("SUM.rds")

  
  
############################################
##OLD/PREVIOUS MODELS (just for reference)##
############################################
#########################################
##### ANNUAL STEMS:SEEDS-OUT LAMBDA ----

###visualizations ----
ggplot(dat, aes(x=density_p, y=out_a, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=1)+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("Na(t+1)")+xlab("Np(t))")

ggplot(dat, aes(x=density_a, y=out_a, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=1)+
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("Na(t+1)")+xlab("Na(t))")

# check replication within blocks
table(dat$warmtrt, dat$block)

### BRM fit AB: SCALED ----
annual_lambda <- brm(bf(out_a ~ lambdaA*3000 / (1+alphaAA*density_a + alphaAP*density_p), 
                        lambdaA ~ warmtrt + (1|block), 
                        alphaAA ~  warmtrt + (1|block), 
                        alphaAP ~  warmtrt + (1|block), nl=TRUE),
                     data = dat,
                     prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaAA"),
                               prior(normal(0, .1), nlpar = "alphaAP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=5000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 16))
annual_lambda
plot(annual_lambda)
fixef(annual_lambda)
conditional_effects(annual_lambda)

# just ambient scaled
annual_lambda.amb <- brm(bf(out_a ~ lambdaA*3000 / (1+alphaAA*density_a + alphaAP*density_p), 
                            lambdaA ~ 1 + (1|block), 
                            alphaAA ~  1 + (1|block), 
                            alphaAP ~  1 + (1|block), nl=TRUE),
                         data = subset(dat, warmtrt="amb"),
                         prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaA"), 
                                   prior(normal(0, .1), nlpar = "alphaAA"),
                                   prior(normal(0, .1), nlpar = "alphaAP")),
                         inits = "0",  cores=4, chains=4, iter=5000, thin=5,
                         control = list(adapt_delta = 0.99, max_treedepth = 16))

# just warmed scaled
annual_lambda.amb <- brm(bf(out_a ~ lambdaA*3000 / (1+alphaAA*density_a + alphaAP*density_p), 
                            lambdaA ~ 1 + (1|block), 
                            alphaAA ~  1 + (1|block), 
                            alphaAP ~  1 + (1|block), nl=TRUE),
                         data = subset(dat, warmtrt="warm"),
                         prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaA"), 
                                   prior(normal(0, .1), nlpar = "alphaAA"),
                                   prior(normal(0, .1), nlpar = "alphaAP")),
                         inits = "0", cores=4, chains=4,iter=5000,thin=5,
                         control = list(adapt_delta = 0.99, max_treedepth = 16))
# just ambient log-scaled
annual.ambient <- brm(bf(log(1+out_a) ~ exp(log(lambdaA) - log((1+alphaAA*(density_a+1) + alphaAP*(density_p+1)))), 
                         lambdaA ~ 1+(1|block),
                         alphaAA ~ 1+ (1|block), 
                         alphaAP ~ 1+  (1|block), nl=TRUE),
                      data = subset(annuals, warmtrt=="amb"),
                      prior = c(prior(normal(8, 2), lb=0, nlpar = "lambdaA"), 
                                prior(normal(0, .1), nlpar = "alphaAA"),
                                prior(normal(0, .1), nlpar = "alphaAP")),
                      inits = "0", cores=4, chains=4, iter=5000, 
                      control = list(adapt_delta = 0.99, max_treedepth = 18))

# just warmed log-scaled
annual.ambient <- brm(bf(log(1+out_a) ~ exp(log(lambdaA) - log((1+alphaAA*(density_a+1) + alphaAP*(density_p+1)))), 
                         lambdaA ~ 1+(1|block),
                         alphaAA ~ 1+(1|block), 
                         alphaAP ~ 1+(1|block), nl=TRUE),
                      data = subset(annuals, warmtrt=="warm"),
                      prior = c(prior(normal(8, 2), lb=0, nlpar = "lambdaA"), 
                                prior(normal(0, .1), nlpar = "alphaAA"),
                                prior(normal(0, .1), nlpar = "alphaAP")),
                      inits = "0",  cores=4, chains=4, iter=5000, 
                      control = list(adapt_delta = 0.99, max_treedepth = 18))


##### ANNUAL SEEDS-IN:STEMS ----
### visualize ----
ggplot(dat.surv.a, aes(x=seeded_a, y=spring20_a)) +
  geom_jitter(aes(color=warmtrt), width=100)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("Na(t)")+xlab("na(t))")

ggplot(dat.surv.a, aes(x=warmtrt, y=spring20_a/seeded_a)) +
  geom_boxplot()+
  # geom_point(color='blue') +
  geom_jitter(width=.1) # 

# check replication within blocks
table(dat.surv$warmtrt, dat.surv$block) # 3 per block


### BRM fit ----
annual_sprsur <- brm(bf(spring20_a|trials(seeded_a) ~ warmtrt + (1|block)), 
                     data = dat.surv.a,
                     family=binomial,
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=5000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 18))

stancode(annual_sprsur)
annual_sprsur
summary(annual_sprsur)
plot(annual_sprsur)
fixef(annual_sprsur)
conditional_effects(annual_sprsur)

# Graph probabilities ~ treatment ----
nd <- tibble(warmtrt = c("amb","warm"), block=NA, seeded_a=1)
n_iter<-100

fitted(annual_sprsur,newdata = nd, summary = TRUE ) %>%
  as_tibble() 

f <-
  fitted(annual_sprsur,
         newdata  = nd,
         summary = F,
         nsamples = n_iter) %>% 
  as_tibble() %>%
  rename(amb = V1, warm=V2) %>%
  pivot_longer(cols=amb:warm, names_to="treatment", values_to="probability")
f

f %>% ggplot(aes(x=treatment, y=probability))+
  stat_pointinterval(.width = c(.66, .95))


## 1.1 Simple (seeds in:out) annuals percapita and scaled----
annual.simple.a <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a + alphaAP*pm2 + alphaAS*seeded_s), #changed from old version: annual.simple <- brm(bf(plotseeds ~ (lambdaA*seeded_a) / (1+alphaAA*seeded_a + alphaAP*pm2), 
                          lambdaA ~ 1 + (1|block), # only lambda vary by warmtrt
                          alphaAA + alphaAP + alphaAS ~ 1 + (1|block), #just one overall alpha, simplify!
                          nl=TRUE), 
                       data = subset(annuals, warmtrt=="amb"),
                       family = poisson, #poisson distribution based on lina's suggestion.
                       prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaA"), 
                                 prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaAA"),
                                 prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaAS"),  #added term for seedling competitive effect
                                 prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaAP")),
                       inits = "0",  
                       cores=4, 
                       chains=4,
                       iter=20000, 
                       thin=1,
                       refresh=100,
                       control = list(adapt_delta = 0.99, max_treedepth = 18))

annual.simple.w <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a + alphaAP*pm2 + alphaAS*seeded_s), #changed from old version: annual.simple <- brm(bf(plotseeds ~ (lambdaA*seeded_a) / (1+alphaAA*seeded_a + alphaAP*pm2), 
                          lambdaA ~ 1 + (1|block), # only lambda vary by warmtrt
                          alphaAA + alphaAP + alphaAS ~ 1 + (1|block), 
                          nl=TRUE), 
                       data = subset(annuals, warmtrt=="warm"),
                       family = poisson, 
                       prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaA"), 
                                 prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaAA"),
                                 prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaAS"),  #added term for seedling competitive effect
                                 prior(normal(0, 1), lb=0, ub=1, nlpar = "alphaAP")),
                       inits = "0",  
                       cores=4, 
                       chains=4,
                       iter=20000, 
                       thin=1,
                       refresh=100,
                       control = list(adapt_delta = 0.99, max_treedepth = 18))

#save both ambient and warmed
saveRDS(annual.simple, file="A063021a.rds") #use this when I have a good run to save it as a file, can add date to filename
saveRDS(annual.simple, file="A063021w.rds") #use this when I have a good run to save it as a file, can add date to filename

#read back in ambient and warmed from most recent saved file
annual.warm.model<-readRDS("A063021w.rds")
annual.amb.model<-readRDS("A063021a.rds")

conditional_effects(annual.warm.model)%>% plot(points=T)



### 1.2 Simple annuals logscaled (SKIP)----
# example from hallett: m1A <- as.formula(log(AVseedout +1) ~ log(ag*(AVseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))
annual.logscale <- brm(bf(log(1+plotseeds) ~ log((seeded_a+1)*exp(log(lambdaA) - log((1+alphaAA*(seeded_a+1) + alphaAP*(pm2+1)+ alphaAS*(sm2+1))))), 
                          lambdaA ~ warmtrt+(1|block),
                          alphaAA ~ warmtrt+ (1|block), 
                          alphaAS ~ warmtrt+ (1|block), 
                          alphaAP ~ warmtrt+  (1|block), nl=TRUE),
                       data = annuals,
                       prior = c(prior(normal(8, 3), lb=0, nlpar = "lambdaA"), 
                                 prior(normal(0, .1), nlpar = "alphaAA"),
                                 prior(normal(0, .1), nlpar = "alphaAS"),
                                 prior(normal(0, .1), nlpar = "alphaAP")),
                       inits = "0",  
                       cores=4, 
                       chains=4,
                       iter=10000, 
                       control = list(adapt_delta = 0.95, max_treedepth = 16))
savedAL<-annual.logscale
#saveRDS(savedAL, file="AL.rds")
#readAL<-readRDS("AL.rds")


### 1.3 Simple annuals percapita & logscaled
annual.logscale.percap <- brm(bf(log(1+percap) ~ log(exp(log(lambdaA) - log((1+alphaAA*(seeded_a+1) + alphaAP*(pm2+1)+ alphaAS*(sm2+1))))), 
                                 lambdaA ~ warmtrt+(1|block),
                                 alphaAA ~ warmtrt+ (1|block), 
                                 alphaAS ~ warmtrt+ (1|block), 
                                 alphaAP ~ warmtrt+  (1|block), nl=TRUE),
                              data = annuals,
                              prior = c(prior(normal(8, 3), lb=0, nlpar = "lambdaA"), 
                                        prior(normal(0, .1), nlpar = "alphaAA"),
                                        prior(normal(0, .1), nlpar = "alphaAS"),
                                        prior(normal(0, .1), nlpar = "alphaAP")),
                              inits = "0",  
                              cores=4, 
                              chains=4,
                              iter=5000, 
                              control = list(adapt_delta = 0.95, max_treedepth = 16))

savedALP<-annual.logscale.percap
#saveRDS(savedALP, file="ALP.rds")
#readAL<-readRDS("ALP.rds")


### Final Params: seedling competitive effects, not sure how these work----
# value from mordecai is .5, oddly much higher than other competitive effects in the model...?
#alpha_as=theta*alpha_aa
#alpha_ps=theta*alpha_pa

#theta~Uniform(0, 1)
#normal(old) seedling model----
#seedling.simple <- brm(bf(fall20_s.g ~ (lambdaS*seeded_s) / (1+alphaSA*seeded_a + alphaSP*density_p+ alphaSS*seeded_s), 
#                        lambdaS ~ warmtrt + (1|block), 
#                        alphaSA ~  warmtrt + (1|block), 
#                        alphaSP ~  warmtrt + (1|block),
#                        alphaSS ~  warmtrt + (1|block),
#                        nl=TRUE),
#                     data = seedlings,
#                     prior = c(prior(normal(.001, .0002), lb=0, nlpar = "lambdaS"), 
#                               prior(normal(0, .1), nlpar = "alphaSA"),
#                               prior(normal(0, .1), nlpar = "alphaSS"),
#                               prior(normal(0, .1), nlpar = "alphaSP")),
#                     inits = "0",  
#                     cores=4, 
#                     chains=4,
#                     iter=100000, 
#                     thin=5,
#                     control = list(adapt_delta = 0.99, max_treedepth = 20))

plot(seedling.simple)

