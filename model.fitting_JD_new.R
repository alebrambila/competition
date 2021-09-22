library(rstan)
library(reshape2)
library(brms)
library(tidybayes)
library(modelr)
library(minpack.lm)
library(nlstools)
library(tidyverse)
library(minpack.lm)
library(grid)
library(gridExtra)

library(patchwork)
library(GGally)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('data.cleaning.R')

# This script is to estimate the parameters necessary to be able to perform invasion analysis and additional 
# simulation of an annual (Lolium multiflorum) and perennial (Festuca roemeri) in competition at various densities.
#  Here I estimate the inherent growth rates and competition coefficients in both warmed and ambient plots.

#Models 
# 1.1 Simple (seeds in:out) annuals percapita and scaled
# 1.2 (skip) Simple annuals logscaled 
# 1.3 (skip) Simple annuals percapita & logscaled
# 2.1 Perennial fecundity scaled
# 2.2 (skip) Perennial fecundity logscaled
# 3.1 Simple (seeds in:adults out) perennial seedlings binomial
# 3.2 (skip) pring (seeds in:stems out) perennial seedlings binomial (first half of 3.1)
# 3.3 (skip) Summer (stems in:adults out) perennial seedlings binomial (second half of 3.1)


###### 1. SIMPLE ANNUALS: SEEDS-IN:SEEDS-OUT ----
### data----
dat<-left_join(fecundity_phyt_a, density_spring20)%>%
  mutate(out_a=seeds, density_a=am2, density_p=pm2, density_s=sm2)%>%
  select(1, 6, 7, 8, out_a, density_a, density_p, density_s)

dat.surv.a <- sprsur2020a<-filter(spr_sur2020, seeded_a!=0)%>%  # just naming the data something different 
  dplyr::select(-spring20_s)%>%
  mutate(sprsur_a=spring20_a/seeded_a)

annuals<-left_join(fecundity_plot_a, dplyr::select(dat.surv.a, 1,9,10))
annuals<-left_join(annuals, plotkey)%>%
  mutate(seeded_a=ifelse(comptrt=="none"|comptrt=='seedling perennials', 4/0.06019467, ifelse(comptrt=="adult perennials", 8/0.06019467, seeded_a)))%>% # germination correction factor for phytometers.  how many stems I had in spring divided by germination factor to give me how many seeds were added.   
  mutate(seeded_s=ifelse(comptrt=="none"|comptrt=='annuals', 4/0.06019467, ifelse(comptrt=="adult perennials", 8/0.06019467, ifelse(comptrt=="seedling perennials", 4500, seeded_s))))%>% # germination correction factor for phytometers.  how many stems I had in spring divided by germination factor to give me how many seeds were added.   
  mutate(plotseeds=seeds*count_a)
annuals<-left_join(annuals, density_spring20)

annuals$block <- as.factor(annuals$block)
annuals$percap <- annuals$plotseeds/annuals$seeded_a


### visualize ----
ggplot(annuals, aes(x=seeded_a, y=plotseeds, color=warmtrt)) +
  geom_point(aes(shape=comptrt), width=1)+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("total annual seed production")+xlab("seeded annuals")

head(annuals); dim(annuals)
# how many plots / block?  answer: 12
annuals %>% group_by(block) %>%
  summarize(n.plots = n_distinct(plotid))
table(annuals$warmtrt, annuals$block)

# per capita
aa <- ggplot(subset(annuals), aes(x=pm2, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F) +xlab("perennial adults")+ylab("annual percapita seeds out")

bb<- ggplot(subset(annuals), aes(x=seeded_a, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+xlab("annual seeds in")+ylab("annual percapita seeds out")

cc<-ggplot(subset(annuals), aes(x=seeded_s, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+xlab("seeded_s")+ylab("annual percapita seeds out")

aa+bb+cc
ggsave(filename="annuals.explore1.pdf", width=12, height=4)


a1<-ggplot(annuals, aes(x=percap)) + geom_histogram() 
a2<-ggplot(annuals, aes(x=seeded_a)) + geom_histogram() 
a3<-ggplot(annuals, aes(x=seeded_s)) + geom_histogram() 
a4<-ggplot(annuals, aes(x=pm2)) + geom_histogram() 

(a1 | a2) / (a3 | a4)
ggsave(filename="annuals.explore2.pdf", width=9, height=9)

ggplot2::theme_set(ggplot2::classic()) 
ggpairs(dplyr::select(annuals, seeded_a, seeded_s, pm2,percap))
ggsave(filename="annuals.pairs.pdf", width=6, height=6)

a<-annuals %>% ggplot(aes(x=seeded_a, y=seeded_s,color=warmtrt))+
  geom_jitter()
b<-annuals %>% ggplot(aes(x=pm2, y=seeded_s,color=warmtrt))+
  geom_jitter()

a+b
ggsave(filename="annuals.explore3.pdf", width=9, height=4)

### gaussian fit ----
## 1.1 Simple (seeds in:out) annuals percapita and scaled

# scale seeded_a
# annuals.backup <- annuals
annuals <- mutate(annuals, seeded_a=seeded_a/1000)
hist(annuals$seeded_a)

annual.simple.gaussian <- brm(bf(percap ~ lambdaA*100 / (1 + alphaAA*seeded_a + alphaAP*pm2 + alphaAS*seeded_s),
                        lambdaA ~ warmtrt , #+ (1|block),
                        alphaAA + alphaAP + alphaAS ~ warmtrt,# + (1|block), #
                        nl=TRUE), 
                     data = annuals,
                     family = gaussian, #poisson, 
                     prior = c(prior(normal(1, 1), nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaAA"),
                               prior(normal(0, .1), nlpar = "alphaAS"), 
                               prior(normal(0, .1), nlpar = "alphaAP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=5000, 
                     thin=1,
                     refresh=100,
                     control = list(adapt_delta = 0.99, max_treedepth = 18))

annual.simple.gaussian
saveRDS(annual.simple.gaussian, file="annual.simple.gaussian.rds")

### poisson fit ----
## 1.1 Simple (seeds in:out) annuals percapita and scaled
annual.simple.poisson <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a + alphaAP*pm2 ), # + alphaAS*seeded_s
                          lambdaA ~ warmtrt + (1|block),
                          alphaAA + alphaAP  ~ warmtrt + (1|block), # + alphaAS
                          nl=TRUE), 
                       data = annuals,
                       family = poisson, #poisson distribution based on lina's suggestion.
                       prior = c(prior(normal(1, 1), nlpar = "lambdaA"), 
                                 prior(normal(0, 1), nlpar = "alphaAA"),
                                 # prior(normal(0, 1), nlpar = "alphaAS"),  #added term for seedling competitive effect
                                 prior(normal(0, 1), nlpar = "alphaAP")),
                       inits = "0",  
                       cores=4, 
                       chains=4,
                       iter=5000, 
                       thin=1,
                       refresh=100,
                       control = list(adapt_delta = 0.99, max_treedepth = 18))

annual.simple.poisson
fixef(annual.simple.poisson)
# annual.simple.poisson <- annual.simple   # backed up before running normal model to see if the issue is the log link
saveRDS(annual.simple.poisson, file="annual.simple.poisson.rds")

conditional_effects(annual.simple.poisson, points=TRUE)



## predict and plot ----

annual.simple.poisson<-readRDS("annual.simple.poisson.rds")
annual.simple.gaussian<-readRDS("annual.simple.gaussian.rds")


# predict over pm2
dat.new.annual <- expand.grid(
  seeded_a = max(annuals$seeded_a, na.rm=TRUE)
  # seeded_a = 0
  # seeded_a = mean(annuals$seeded_a, na.rm=TRUE) #0
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) # 0
  ,pm2 = seq(0,10, length.out=20)
  ,warmtrt = c("amb","warm")
)

pred.annual.gaussian <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual) 

pred.annual <- pred.annual.gaussian %>% 
  mutate(poisson = pred.annual.poisson$Estimate)
head(pred.annual); dim(pred.annual)
pred.annual <- pred.annual %>% rename(gaussian = Estimate) %>%
  pivot_longer(cols=c(gaussian, poisson), names_to="model", values_to = "percap_pred")
pred.annual

a <- pred.annual %>%  ggplot(aes(x = pm2, y = percap_pred, color=warmtrt)) + 
  geom_line(aes(linetype=model)) + 
  geom_jitter(data=annuals, aes(x=pm2, y=percap))+
  ylab("Per capita fecundity")+
  theme_pander()
a
# 

# do for seeded_s

# dat.new.annual <- expand.grid(
#   seeded_a = mean(annuals$seeded_a, na.rm=TRUE) #0
#   ,seeded_s= seq(0,max(annuals$seeded_s), length.out=100)
#   ,pm2 = mean(annuals$pm2, na.rm=TRUE) # 0
#   ,warmtrt = c("amb","warm")
# )
# 
# pred.annual <- 
#   as.data.frame(predict(annual.simple, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
#   cbind(dat.new.annual) 
# 
# b <- pred.annual %>%  ggplot(aes(x = seeded_s, y = Estimate, color=warmtrt)) + 
#   geom_line()+
#   geom_jitter(data=annuals, aes(x=seeded_s, y=percap))+
#   ylab("Per capita fecundity")+
#   theme_pander()
# # 


# predict over seeded_a

dat.new.annual <- expand.grid(
  seeded_a = seq(0,max(annuals$seeded_a), length.out=100)
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) #0
  ,pm2 = mean(annuals$pm2, na.rm=TRUE) # 0
  ,warmtrt = c("amb","warm")
)

pred.annual.gaussian <- 
  as.data.frame(predict(annual.simple.gaussian, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual) 

pred.annual.poisson <- 
  as.data.frame(predict(annual.simple.poisson, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual) 

pred.annual <- pred.annual.gaussian %>% 
  mutate(poisson = pred.annual.poisson$Estimate)
pred.annual <- pred.annual %>% rename(gaussian = Estimate) %>%
  pivot_longer(cols=c(gaussian, poisson), names_to="model", values_to = "percap_pred")
head(pred.annual); dim(pred.annual)


c <- pred.annual.gaussian %>%  ggplot(aes(x = seeded_a, y = percap_pred, color=warmtrt)) + 
  geom_line(aes()) +
  geom_jitter(data=annuals, aes(x=seeded_a, y=percap))+
  ylab("Per capita fecundity")+
  theme_pander()
#
c
a+c
# a/b/c
# a+b+c
ggsave(filename="annuals.poisson.predict.pdf", width=12, height=5)


## Get parameters ----

get_variables(annual.simple.gaussian)

fit.annuals.table.gaussian <- annual.simple.gaussian %>%
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
               names_sep = "_") %>% 
  group_by(param, treatment) %>% 
  # median_qi(.width = c(.95)) 
  mean_qi(.width = c(.95)) 

fit.annuals.table.gaussian


param.amb.gaus <- filter(fit.annuals.table.gaussian, treatment=="amb")
param.warm.gaus <- filter(fit.annuals.table.gaussian, treatment=="warm")

# try plotting with parameter estimates ----
# annual.simple <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a + alphaAP*pm2 + alphaAS*seeded_s)
# bev1 <- function(x, lam, AA,AS,N_S,AP,N_P) {lam  / (1+AA*x + AS*N_S + AP*N_P)} # Beverton-Holt
bev1 <- function(x, lam, AA,AP,N_P) {lam  / (1+AA*x + AP*N_P)} # Beverton-Holt

 ggplot(annuals, aes(x = seeded_a, y = percap, color=warmtrt)) + 
  geom_jitter()+
  geom_function(fun=bev1, color='red', size=1.5, args=c("lam"=100*param.amb.gaus$value[3], "AA"=param.amb.gaus$value[1], 
                                                          # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                          "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  # try with fixed values for comparison
  #geom_function(fun=bev1, color='orange', size=1.5, args=c("lam"=100, "AA"=.5,
  #                                                        # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
  #                                                        "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  #)) +
  geom_function(fun=bev1, color='blue', size=1.5, args=c("lam"=100*param.warm.gaus$value[3], "AA"=param.warm.gaus$value[1], 
                                                        # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                        "AP"=param.warm.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  # geom_function(fun=bev1, color='blue', size=1.5, args=c("lam"=45, "AA"=param.amb.gaus$value[1], 
  #                                                         # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
  #                                                         "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  # )) +
  ylab("Per capita fecundity")#+
  theme_pander()
cc
# 
# Looks like estimate for AA is way off

pred.annual %>%  ggplot(aes(x = pm2, y = percap_pred, color=warmtrt)) + 
  geom_line(aes(linetype=model)) + 
  geom_jitter(data=annuals, aes(x=pm2, y=percap))+
  geom_function(fun=bev1, color='red', size=1.5, args=c("lam"=100*param.amb.gaus$value[3], "AA"=param.amb.gaus$value[1], 
                                                        # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                        "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  geom_function(fun=bev1, color='blue', size=1.5, args=c("lam"=100*param.warm.gaus$value[3], "AA"=param.warm.gaus$value[1], 
                                                         # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                         "AP"=param.warm.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  ylab("Per capita fecundity")+
  theme_pander()
#


annual.simple.gaussian
plot(annual.simple.gaussian)



## !! Version with single covariate ----


### gaussian fit ----
## 1.1 Simple (seeds in:out) annuals percapita and scaled

# scale seeded_a
annuals.backup <- annuals
annuals <- mutate(annuals, seeded_a=seeded_a/1000)
hist(annuals$seeded_a)

annual.gaussian.1var <- brm(bf(percap ~ lambdaA*100 / (1 + alphaAA*seeded_a), # + alphaAS*seeded_s 
                                 lambdaA ~ warmtrt, #+ (1|block),
                                 alphaAA ~ warmtrt, #+ (1|block), #+ alphaAS
                                 nl=TRUE), 
                              data = annuals,
                              family = gaussian, #poisson, 
                              prior = c(prior(normal(1, .3), nlpar = "lambdaA"), 
                                        prior(normal(0, .1), nlpar = "alphaAA")),
                                        # prior(normal(0, .1), nlpar = "alphaAS"),  #added term for seedling competitive effect
                                        # prior(normal(0, .1), nlpar = "alphaAP")),
                              inits = "0",  
                              cores=4, 
                              chains=4,
                              iter=5000, 
                              thin=1,
                              refresh=100,
                              control = list(adapt_delta = 0.99, max_treedepth = 18))

annual.gaussian.1var
saveRDS(annual.gaussian.1var, file="annual.gaussian.1var.rds")

### poisson fit ----
## 1.1 Simple (seeds in:out) annuals percapita and scaled
annual.poisson.1var <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a  ), # + alphaAP*pm2+ alphaAS*seeded_s
                                lambdaA ~ warmtrt + (1|block),
                                alphaAA ~ warmtrt + (1|block), # + alphaAS
                                nl=TRUE), 
                             data = annuals,
                             family = poisson, #poisson distribution based on lina's suggestion.
                             prior = c(prior(normal(1, 1), nlpar = "lambdaA"), 
                                       prior(normal(0, 1), nlpar = "alphaAA")),
                                       # prior(normal(0, 1), nlpar = "alphaAS"),  #added term for seedling competitive effect
                                       # prior(normal(0, 1), nlpar = "alphaAP")),
                             inits = "0",  
                             cores=4, 
                             chains=4,
                             iter=5000, 
                             thin=1,
                             refresh=100,
                             control = list(adapt_delta = 0.99, max_treedepth = 18))

annual.poisson.1var
fixef(annual.poisson.1var)
# annual.simple.poisson <- annual.simple   # backed up before running normal model to see if the issue is the log link
saveRDS(annual.poisson.1var, file="annual.poisson.1var.rds")

conditional_effects(annual.simple.poisson, points=TRUE)



## predict and plot ----

annual.poisson.1var <-readRDS("annual.poisson.1var.rds")
annual.gaussian.1var <-readRDS("annual.gaussian.1var.rds")



# predict over seeded_a

dat.new.annual <- expand.grid(
  seeded_a = seq(0,max(annuals$seeded_a), length.out=100)
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) #0
  # ,pm2 = mean(annuals$pm2, na.rm=TRUE) # 0
  ,warmtrt = c("amb","warm")
)

pred.annual.gaussian <- 
  as.data.frame(predict(annual.gaussian.1var, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual) 

pred.annual.poisson <- 
  as.data.frame(predict(annual.poisson.1var, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual) 

pred.annual <- pred.annual.gaussian %>% 
  mutate(poisson = pred.annual.poisson$Estimate)
pred.annual <- pred.annual %>% rename(gaussian = Estimate) %>%
  pivot_longer(cols=c(gaussian, poisson), names_to="model", values_to = "percap_pred")
head(pred.annual); dim(pred.annual)


c <- pred.annual %>%  ggplot(aes(x = seeded_a, y = percap_pred, color=warmtrt)) + 
  geom_line(aes(linetype=model)) +
  geom_jitter(data=annuals, aes(x=seeded_a, y=percap))+
  ylab("Per capita fecundity")+
  theme_pander()
#
c


## Get parameters ----

get_variables(annual.gaussian.1var)

fit.annuals.table.gaussian <- annual.gaussian.1var %>%
  spread_draws(`b_.*`, regex = TRUE) %>% 
  mutate(
    lam_amb=b_lambdaA_Intercept,
    lam_warm=b_lambdaA_Intercept + b_lambdaA_warmtrtwarm,
    
    alphaAA_amb=b_alphaAA_Intercept,
    alphaAA_warm=b_alphaAA_Intercept + b_alphaAA_warmtrtwarm,
    
    # alphaAP_amb=b_alphaAP_Intercept,
    # alphaAP_warm=b_alphaAP_Intercept + b_alphaAP_warmtrtwarm,
    
    # alphaAS_amb=b_alphaAS_Intercept,
    # alphaAS_warm=b_alphaAS_Intercept + b_alphaAS_warmtrtwarm
  ) %>%
  dplyr::select(-contains("b_")) %>% 
  pivot_longer(-c(.chain, .iteration, .draw), 
               names_to = c("param", "treatment"),
               names_sep = "_") %>% 
  group_by(param, treatment) %>% 
  # median_qi(.width = c(.95)) 
  mean_qi(.width = c(.95)) 

fit.annuals.table.gaussian


param.amb.gaus <- filter(fit.annuals.table.gaussian, treatment=="amb")
param.warm.gaus <- filter(fit.annuals.table.gaussian, treatment=="warm")

# try plotting with parameter estimates ----
# annual.simple <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a + alphaAP*pm2 + alphaAS*seeded_s)
# bev1 <- function(x, lam, AA,AS,N_S,AP,N_P) {lam  / (1+AA*x + AS*N_S + AP*N_P)} # Beverton-Holt
bev1 <- function(x, lam, AA) {lam  / (1+AA*x )} # Beverton-Holt

cc <- annuals %>%  ggplot(aes(x = seeded_a, y = percap, color=warmtrt)) + 
  geom_jitter()+
  geom_function(fun=bev1, color='red', size=1.5, args=c("lam"=100*param.amb.gaus$value[2], "AA"=param.amb.gaus$value[1]
                                                        # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                        # "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  # try with fixed values for comparison
  geom_function(fun=bev1, color='orange', size=1.5, args=c("lam"=100, "AA"=.5
                                                           # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                           # "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  # geom_function(fun=bev1, color='blue', size=1.5, args=c("lam"=45, "AA"=param.amb.gaus$value[1], 
  #                                                         # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
  #                                                         "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  # )) +
  ylab("Per capita fecundity")+
  theme_pander()
cc
# 
# Looks like estimate for AA is way off



## !!!! Version with only Ambient, single covariate ----


### gaussian fit ----
## 1.1 Simple (seeds in:out) annuals percapita and scaled

# scale seeded_a
hist(annuals$seeded_a)

annuals_amb <- filter(annuals, warmtrt == "amb")

ggplot(annuals_amb, aes(x=seeded_a, y=percap))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~block)


ggpairs(select(annuals_amb, seeded_a, percap))

annual.1var.amb.gaussian <- brm(bf(percap ~ lambdaA*100 / (1 + alphaAA*seeded_a), # + alphaAS*seeded_s 
                                 lambdaA ~ 1#  (1|block),
                                 ,alphaAA ~ 1#(1|block), #+ alphaAS
                                 ,nl=TRUE), 
                              data = annuals_amb,
                              family = gaussian, #poisson, 
                              prior = c(prior(normal(1, 1), nlpar = "lambdaA"), 
                                        prior(normal(0, .1), nlpar = "alphaAA")),
                              # prior(normal(0, .1), nlpar = "alphaAS"),  #added term for seedling competitive effect
                              # prior(normal(0, .1), nlpar = "alphaAP")),
                              inits = "0",  
                              cores=4, 
                              chains=4,
                              iter=5000, 
                              thin=1,
                              refresh=100,
                              control = list(adapt_delta = 0.99, max_treedepth = 18))

annual.1var.amb.gaussian
saveRDS(annual.1var.amb.gaussian, file="annual.1var.amb.gaussian.rds")

### poisson fit ----
## 1.1 Simple (seeds in:out) annuals percapita and scaled
annual.1var.amb.poisson <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a  ), # + alphaAP*pm2+ alphaAS*seeded_s
                                lambdaA ~  (1|block),
                                alphaAA ~  (1|block), # + alphaAS
                                nl=TRUE), 
                             data = annuals_amb,
                             family = poisson, 
                             prior = c(prior(normal(1, 1), nlpar = "lambdaA"), 
                                       prior(normal(0, 1), nlpar = "alphaAA")),
                             # prior(normal(0, 1), nlpar = "alphaAS"),  #added term for seedling competitive effect
                             # prior(normal(0, 1), nlpar = "alphaAP")),
                             inits = "0",  
                             cores=4, 
                             chains=4,
                             iter=5000, 
                             thin=1,
                             refresh=100,
                             control = list(adapt_delta = 0.99, max_treedepth = 18))

annual.1var.amb.poisson
fixef(annual.1var.amb.poisson)
# annual.simple.poisson <- annual.simple   # backed up before running normal model to see if the issue is the log link
saveRDS(annual.1var.amb.poisson, file="annual.1var.amb.poisson.rds")

conditional_effects(annual.simple.poisson, points=TRUE)



## predict and plot ----

annual.1var.amb.poisson <-readRDS("annual.1var.amb.poisson.rds")
annual.1var.amb.gaussian <-readRDS("annual.1var.amb.gaussian.rds")



# predict over seeded_a

dat.new.annual <- expand.grid(
  seeded_a = seq(0,max(annuals$seeded_a), length.out=100)
  # ,seeded_s= mean(annuals$seeded_s, na.rm=TRUE) #0
  # ,pm2 = mean(annuals$pm2, na.rm=TRUE) # 0
  # ,warmtrt = c("amb","warm")
)

pred.1var.amb.gaussian <- 
  as.data.frame(predict(annual.1var.amb.gaussian, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual) 

pred.1var.amb.poisson <- 
  as.data.frame(predict(annual.1var.amb.poisson, newdata = dat.new.annual, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new.annual) 

pred.annual <- pred.1var.amb.gaussian %>% 
  mutate(poisson = pred.1var.amb.poisson$Estimate)
pred.annual <- pred.annual %>% rename(gaussian = Estimate) %>%
  pivot_longer(cols=c(gaussian, poisson), names_to="model", values_to = "percap_pred")
head(pred.annual); dim(pred.annual)


c <- pred.annual %>%  ggplot(aes(x = seeded_a, y = percap_pred)) + 
  geom_line(aes(linetype=model)) +
  geom_jitter(data=annuals_amb, aes(x=seeded_a, y=percap))+
  ylab("Per capita fecundity")+
  theme_pander()
#
c


## Get parameters ----

get_variables(annual.1var.amb.gaussian)

fit.annuals.table.gaussian <- annual.1var.amb.gaussian %>%
  spread_draws(`b_.*`, regex = TRUE) %>% 
  mutate(
    lam_amb=b_lambdaA_Intercept,
    alphaAA_amb=b_alphaAA_Intercept,

    # alphaAP_amb=b_alphaAP_Intercept,
    # alphaAP_warm=b_alphaAP_Intercept + b_alphaAP_warmtrtwarm,
    
    # alphaAS_amb=b_alphaAS_Intercept,
    # alphaAS_warm=b_alphaAS_Intercept + b_alphaAS_warmtrtwarm
  ) %>%
  dplyr::select(-contains("b_")) %>% 
  pivot_longer(-c(.chain, .iteration, .draw), 
               names_to = c("param")) %>% 
  group_by(param) %>% 
  # median_qi(.width = c(.95)) 
  mean_qi(.width = c(.95)) 

fit.annuals.table.gaussian


param.amb.gaus <- fit.annuals.table.gaussian

# try plotting with parameter estimates ----
# annual.simple <- brm(bf(as.integer(percap) ~ lambdaA*100 / (1 + alphaAA*seeded_a + alphaAP*pm2 + alphaAS*seeded_s)
# bev1 <- function(x, lam, AA,AS,N_S,AP,N_P) {lam  / (1+AA*x + AS*N_S + AP*N_P)} # Beverton-Holt
bev1 <- function(x, lam, AA) {lam  / (1+AA*x )} # Beverton-Holt

cc <- annuals_amb %>%  ggplot(aes(x = seeded_a, y = percap, color=warmtrt)) + 
  geom_jitter()+
  geom_function(fun=bev1, color='red', size=1.5, args=c("lam"=100*param.amb.gaus$value[2], "AA"=param.amb.gaus$value[1]
                                                        # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                        # "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  # try with fixed values for comparison
  geom_function(fun=bev1, color='orange', size=1.5, args=c("lam"=100, "AA"=1
                                                           # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
                                                           # "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  )) +
  # geom_function(fun=bev1, color='blue', size=1.5, args=c("lam"=45, "AA"=param.amb.gaus$value[1], 
  #                                                         # "AS"=param.amb$value[3], "N_S"=mean(annuals$seeded_s),
  #                                                         "AP"=param.amb.gaus$value[2], "N_P"=mean(annuals$pm2)
  # )) +
  ylab("Per capita fecundity")+
  theme_pander()
cc
# 
#  estimate for AA is still way off


param.amb.gaus

pred.annual %>%  ggplot(aes(x = seeded_a, y = percap_pred)) + 
  geom_line(aes(linetype=model)) +
  geom_jitter(data=annuals_amb, aes(x=seeded_a, y=percap))+
  geom_function(fun=bev1, color='orange', size=1.5, args=c("lam"=38, "AA"=.0215)) +
  ylab("Per capita fecundity")+
  theme_pander()



