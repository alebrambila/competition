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
#annuals.sc <- mutate(annuals, seeded_a=seeded_a/1000)
testannuals<-

annual.simple.poisson <- brm(bf(as.integer(percap) ~  lambdaA / (1 + alphaAA*seeded_am2 + alphaAP*starting_pm2), # + alphaAS*seeded_s #switched seeded_a to seeded_am2 to get density, not sure if to switch pm2 to starting_pm2
                                 lambdaA + alphaAA + alphaAP ~ 0 +warmtrt + time,
                                 # + time,
                                 nl=TRUE), 
                              data = annuals,
                              family = poisson, 
                              prior = c(prior(normal(0, 1), nlpar = "lambdaA"),
                                        prior(normal(.1, .1), nlpar = "alphaAA"),
                                        prior(normal(.1, .1), nlpar = "alphaAP")),
                            #  inits="0",
                              cores=4, 
                              chains=4,
                              iter=10000, 
                              thin=1,
                              refresh=100,
                              control = list(adapt_delta = 0.9, max_treedepth = 16))

annual.simple.poisson <- brm(bf(round(percap) ~  lambdaA / (1 + alphaAA*seeded_am2 + alphaAP*starting_pm2), # *100
                                lambdaA + alphaAA + alphaAP  ~ 0 + warmtrt + time,
                                nl=TRUE),
                             data = annuals,
                             family = poisson,
                             prior = c(prior(normal(4, 4), nlpar = 'lambdaA'),
                                       prior(normal(0.1, .1), nlpar = 'alphaAP'),
                                       prior(normal(0.1, .1), nlpar = 'alphaAA')
                             ),
                             inits='0',
                              cores=4,
                              chains=4,
                              iter=5000,
                              thin=1,
                              refresh=100,
                              control = list(adapt_delta = 0.9, max_treedepth = 16))


plot(annual.simple.poisson)

pairs(annual.simple.poisson, pars = c('b_lambdaA_warmtrtamb', 'b_lambdaA_warmtrtwarm', 'b_alphaAA_warmtrtamb','b_alphaAA_warmtrtwarm', 'b_alphaAP_warmtrtamb','b_alphaAP_warmtrtwarm'))
pairs(x, pars = NA, variable = NULL, regex = FALSE, fixed = FALSE, ...)
#look at the plot of relationship between no competition lambda and intraspecific competition
# pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
# pairs(annual.simple.gaussian, pars = c(“lambdaA”, “alphaAA”, “alphaAP”)
saveRDS(annual.simple.poisson, file="a032222b.rds")
annual.simple.poisson<-readRDS(file="a032222b.rds")

## Get parameters ----
get_variables(annual.simple.gaussian)

plot(conditional_effects(annual.simple.gaussian), ask = FALSE)



## predict and plot COUNTERFACTUALS----

### predict over pm2 (hold seeded_a steady)
dat.new.annual.max <- expand.grid(
  seeded_am2 = max(annuals$seeded_am2, na.rm=TRUE)
  # seeded_a = 0
  # seeded_a = mean(annuals$seeded_a, na.rm=TRUE) #0
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) # 0
  ,starting_pm2 = seq(0,10, length.out=200)
  ,warmtrt = c("amb","warm")
)

dat.new.annual.mean <- expand.grid(
  #seeded_a = max(annuals$seeded_a, na.rm=TRUE)
  # seeded_a = 0
  seeded_am2 = mean(annuals$seeded_am2, na.rm=TRUE) #0
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) # 0
  ,starting_pm2 = seq(0,10, length.out=200)
  ,warmtrt = c("amb","warm")
)

dat.new.annual.0 <- expand.grid(
  seeded_am2 = 50
  # seeded_a = 0
  # seeded_a = mean(annuals$seeded_a, na.rm=TRUE) #0
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) # 0
  ,starting_pm2 = seq(0,10, length.out=200)
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

# extract parameters
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
dat_p<-filter(dat_p, !id%in% c("s1", "s2", "s3", "s4"))
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
hist(dat_p$starting_pm2)
hist(dat_p$seeded_am2)

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
ggplot(dat_p, aes(x= starting_pm2, y=fecundity)) +
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

adult.simple.poisson <- brm(bf(fecundity ~ lambdaP / (1+alphaPA*seeded_am2 + alphaPP*starting_pm2), # lambdaP*20000
                                lambdaP + alphaPA + alphaPP ~ 0 + warmtrt + time, 
                                nl=TRUE), 
                             data = dat_p,
                             family = poisson, 
                             prior = c(prior(normal(10, 10),nlpar = "lambdaP"), 
                                       prior(normal(.1, .1), nlpar = "alphaPA"),
                                       prior(normal(.1, .1), nlpar = "alphaPP")),
                         #    inits = "0",  
                             cores=3, 
                             chains=3,
                             iter=3000, 
                             thin=1,
                             control = list(adapt_delta = 0.9, max_treedepth = 16),
                             refresh=100,
)

conditional_effects(adult.simple.poisson)
adult.simple.poisson
plot(adult.simple.poisson)
conditional_effects(adult.simple.poisson)

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

## Get parameters 
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

#find germination fraction
seedlings$fall.g<-as.integer(seedlings$fall.g)
seedling.simple.poisson<- brm(bf(fall.g  ~ survivalS * seeded_sm2 / (1+alphaSA*seeded_am2 + alphaSS*seeded_sm2 + alphaSP*starting_pm2), #update data terms with new from above
                                  survivalS +alphaSA +alphaSP+alphaSS~ 0 + warmtrt + time, 
                                  nl=TRUE),
                               family=poisson,
                               data = seedlings,   #running this with limited dataset as in the figures above (only in seedling comptrts)
                               prior = c(prior(normal(.0, 1), nlpar = "survivalS"), 
                                         prior(normal(0.1, .1),    nlpar = "alphaSA"),
                                         prior(normal(0.1, .1),    nlpar = "alphaSS"),
                                         prior(normal(0.1, .1),    nlpar = "alphaSP")),
                               inits = "0",  
                               cores=3, 
                               chains=3,
                               iter=5000, 
                               thin=2,
                               control = list(adapt_delta = 0.9, max_treedepth = 16))

plot(seedling.simple.poisson)
seedling.simple.poisson
saveRDS(seedling.simple.gaussian, file="s032222.rds")
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

## Get parameters ----
get_variables(seedling.simple.gaussian.amb)
get_variables(seedling.simple.gaussian.amb)


allfits.seedling <- seedling.simple.gaussian %>%
  spread_draws(`b_.*`, regex = TRUE) %>% 
  mutate(
    lam_amb=b_lambdaS_Intercept,
    lam_warm=b_lambdaS_Intercept + b_lambdaS_warmtrtwarm,  #lam_warm=b_lambdaS_Intercept + allfits.seedling.gaussian$b_lambdaS_Intercept,
    
    alphaSA_amb=b_alphaSA_Intercept,   
    alphaSA_warm=b_alphaSA_Intercept + b_alphaSA_warmtrtwarm, #alphaSA_warm=b_alphaSA_Intercept + allfits.seedling.gaussiana$b_alphaSA_Intercept,
    
    alphaSP_amb=b_alphaSP_Intercept,
    alphaSP_warm=b_alphaSP_Intercept + b_alphaSP_warmtrtwarm, #alphaSP_warm=b_alphaSP_Intercept + allfits.seedling.warm$b_alphaSP_Intercept,
    
    alphaSS_amb=b_alphaSS_Intercept,
    alphaSS_warm=b_alphaSS_Intercept + b_alphaSS_warmtrtwarm #alphaSS_warm=b_alphaSS_Intercept + allfits.seedling.warm$b_alphaSS_Intercept
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


savedps2<-readRDS("PS.rds") #previous best, model w/ warmed and ambient

saveRDS(seedling.binomial.a, file="PS070321a")
saveRDS(seedling.binomial.w, file="PS070321w")

seedling.a<-readRDS("PS070321a")
seedling.w<-readRDS("PS070321w")


plot(seedling.simple)