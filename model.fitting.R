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

# This script is to estimate the parameters necessary to be able to perform invasion analysis and additional 
# simulation of an annual (Lolium multiflorum) and perennial (Festuca roemeri) in competition at various densities.
#  Here I estimate the inherent growth rates and competition coefficients in both warmed and ambient plots.

### Read in needed data ----
### I need to model EACH 'species' parameters in each warming treatment
### Annuals, perennial seedlings and perennial adults. 


### Density is measured as cover and counts/m2. It is a representation of the competitive context in each given plot.
### I need to correct counts to density per meter squared for some plots.  Full plots are 1m2 and 50:50 plots are .66m2. 
### This is for inclusion in formulas to estimate parameters including:
  # lambda_a, lambda_p, lambda_s and all combinations of alpha excluding alpha_xs.
density_spring20 <- right_join(plotkey, dplyr::select(vegplot2020, plotid, time, per_a, per_p, per_s, count_a, count_p, count_s))%>%
  mutate(am2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_a, count_a/.66))%>%
  mutate(pm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_p, count_p/.66))%>%
  mutate(sm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_s, count_s/.66))
#density_spring21 <- ### TO BE GATHERED IN JUNE 2021
#am2, pm2 and sm2 are the relevant data here. 
  
### Fecundity
# to actually run fecundity parameter estimation I need to join fecundity data with density data above
  
# annual seed production per capita (how do combine both of these?)
# phytometers
fecundity_phyt_a <- filter(fecundity, type=="a")

#plot level
fecundity_plot_a<-fecundity_plot

# plug seed production per capita
#phytometers
fecundity_phyt_p <-filter(fecundity, type=="p")
  


### seedling summer survival: 
# seedling counts in spring and fall
# this includes phytometers and background within each plot
# I could pull out phytometers to see if size, immediate neighbors matter
sumsur2020s<-right_join(plotkey, seedling_sumsur2020)
#sumsur2021s <- ### TO BE GATHERED JUNE AND SEPTEMBER 2021
#spring20_s and fall20_s.g are the relevant data here


#annual germination and spring survival
# cannot include phytometers! they are multiple seeded to ensure data. 
sprsur2020a<-filter(spr_sur2020, seeded_a!=0)%>%
  dplyr::select(-seeded_s, -spring20_s)%>%
  mutate(sprsur_a=spring20_a/seeded_a)
  
#seedling germination and seed survival
sprsur2020s <-filter(spr_sur2020, seeded_s!=0)%>%
  dplyr::select(-seeded_a, -spring20_a)%>%
  mutate(sprsur_s=spring20_s/seeded_s)




### Parameters ----
## Annuals Lambda and alphas: Seed production.  This is affected by competitors alphas and densities.  
## Annual germination + spring survival: Binomial distribution only dependent on seeds in. NOT dependent on competition (to start?)
## Perennial seedling germination + spring survival: Same as above.  Not dependent on competition. 
## Adult summer survival: Same as above.  Different time frame. 
## Adult lambda and alphas: like annuals above. 
## Perennial seedling summer survival: Essentially the lambda term for seedlings.  This is fit with densities and competitive effects. 
    #This is (only..?) where we fit competitive effects of others on seedlings. 
## Per capita effects of seedlings on adults: I dont completely understand this one. 
    # We need it for the full model but it seems to be based on the strength of annual competition somehow

#################
### ANNUALS ----

### Lambda and Alphas ----
dat<-left_join(fecundity_phyt_a, density_spring20)%>%
  mutate(out_a=seeds, density_a=am2, density_p=pm2)%>%
  select(1, 6, 7, 8, out_a, density_a, density_p)

#visualizations (annual fecundity as a function of perennial density (1), and annual density (2))
ggplot(dat, aes(x=out_a, y=density_p, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)

ggplot(dat, aes(x=density_a, y=out_a, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)

# check replication within blocks
table(dat$warmtrt, dat$block)

## NOTES / CHANGES: 
# parameter names can't have underscore or dots
# it's not liking the gamma priors (not clear to me why), so changed to normal

# parameter estimation:
annual_lambda <- brm(bf(out_a ~ lambdaA / (1+alphaAA*density_a + alphaAP*density_p), 
                        lambdaA ~ warmtrt + (1|block), 
                        alphaAA ~  warmtrt + (1|block), 
                        alphaAP ~  warmtrt + (1|block), nl=TRUE),
                     data = dat,
                     prior = c(prior(normal(300, 50), lb=0, nlpar = "lambdaA"), 
                               # prior(gamma(3, .01), nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaAA"),
                               prior(normal(0, .1), nlpar = "alphaAP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=5000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 12))
annual_lambda
plot(annual_lambda)
fixef(annual_lambda)
conditional_effects(annual_lambda)

### annual germination/spring survival ----
dat.surv.a <- sprsur2020a<-filter(spr_sur2020, seeded_a!=0)%>%  # just naming the data something different 
  dplyr::select(-seeded_s, -spring20_s)%>%
  mutate(sprsur_a=spring20_a/seeded_a)

dat.surv.a
# quick graphs of variables going into model
ggplot(dat.surv.a, aes(x=seeded_a, y=spring20_a)) +
  geom_jitter(aes(color=warmtrt), width=100)

ggplot(dat.surv.a, aes(x=warmtrt, y=spring20_a/seeded_a)) +
  geom_boxplot()+
  # geom_point(color='blue') +
  geom_jitter(width=.1) # 

# check replication within blocks
table(dat.surv$warmtrt, dat.surv$block) # 3 per block


###### fit simple annual binomial model
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

nd <- tibble(warmtrt = c("amb","warm"), block=NA, seeded_a=1)
n_iter<-100

fitted(annual_sprsur,newdata = nd, summary = TRUE ) %>%
  as_tibble() 

# Graph probabilities ~ treatment

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

##########################################
### Adult Perennial Parameter Fitting ----
### Adult Perennial Lambda/Alphas ----
datp<-left_join(fecundity_phyt_p, density_spring20)%>%
  mutate(out_p=seeds, density_a=am2, density_p=pm2)%>%
  select(1, 6, 7, 8, out_p, density_a, density_p)

#visualizations (adult perennial fecundity as a function of perennial density (1), and annual density (2))
ggplot(datp, aes(x=density_p, y=out_p, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)

ggplot(datp, aes(x=density_a, y=out_p, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)

perennial_lambda <- brm(bf(out_p ~ lambdaP / (1+alphaPA*density_a + alphaPP*density_p), 
                        lambdaP ~ warmtrt + (1|block), 
                        alphaPA ~  warmtrt + (1|block), 
                        alphaPP ~  warmtrt + (1|block), nl=TRUE),
                     data = datp,
                     prior = c(prior(normal(300, 50), lb=0, nlpar = "lambdaP"), 
                               # prior(gamma(3, .01), nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaPA"),
                               prior(normal(0, .1), nlpar = "alphaPP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=5000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 18))

perennial_lambda
plot(perennial_lambda)
fixef(perennial_lambda)
conditional_effects(perennial_lambda)

### Adult summer survival (DONT DO IT FOR NOW - 100% survival when not gophered)
# fit simple binomial model
#annual_sprsur <- brm(bf(FALLP|trials(SPRINGP) ~ warmtrt + (1|block)), 
#                     data = dat.surv,
 #                    family=binomial,
  #                   inits = "0",  
   #                  cores=4, 
    #                 chains=4,
     #                iter=5000, 
      #               thin=5,
       #              control = list(adapt_delta = 0.99, max_treedepth = 18))


#TRY WITH BERNOULLI

# mean value from mordecai is .88, i had 100%
# priors a~Gamma(8, 1) b~Gamma(.01, 1)


### Seedling Perennial Parameter Fitting ----

# Seedling Summer Survival ----
dats<-left_join(select(seedling_sumsur2020, -time), density_spring20)%>%
  mutate(density_a=am2, density_p=pm2, density_s=pm2)%>%
  select(1, 8, 9, 10, fall20_s.g, spring20_s, density_a, density_p, density_s)


#visualizations (annual fecundity as a function of perennial density (1), and annual density (2))
ggplot(dats, aes(x=density_a, y=fall20_s.g/spring20_s, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)

ggplot(dats, aes(x=density_p, y=fall20_s.g/spring20_s, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)

seedling_sumsur <- brm(bf(fall20_s.g|trials(spring20_s) ~ sumsurS / (1+alphaSA*density_a + alphaSS*density_s + alphaSP*density_p), 
                        sumsurS ~ warmtrt + (1|block), 
                        alphaSA ~  warmtrt + (1|block), 
                        alphaSP ~  warmtrt + (1|block), 
                        alphaSS ~  warmtrt + (1|block), nl=TRUE),
                       family=binomial,
                     data = dats,
                     prior = c(prior(normal(1, 1), lb=0, nlpar = "sumsurS"), 
                               # prior(gamma(3, .01), nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaSA"),
                               prior(normal(0, .1), nlpar = "alphaSS"),
                               prior(normal(0, .1), nlpar = "alphaSP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=5000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 18))



### Seedling Spring Survival              - is .34 in mordecai -----
dat.surv <- sprsur2020s<-filter(spr_sur2020, seeded_s!=0)%>%  # just naming the data something different 
  dplyr::select(-seeded_a, -spring20_a)%>%
  mutate(sprsur_s=spring20_s/seeded_s)

dat.surv

# quick graphs of variables going into model
ggplot(dat.surv, aes(x=seeded_s, y=spring20_s)) +
  geom_jitter(aes(color=warmtrt), width=100)

ggplot(dat.surv, aes(x=warmtrt, y=spring20_s/seeded_s)) +
  geom_boxplot()+
  # geom_point(color='blue') +
  geom_jitter(width=.1) # 

# check replication within blocks
table(dat.surv$warmtrt, dat.surv$block) # 3 per block

# fit simple binomial model
seedling_sprsur <- brm(bf(spring20_s|trials(seeded_s) ~ warmtrt + (1|block)), 
                     data = dat.surv,
                     family=binomial,
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=5000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 18))

seedling_sprsur
summary(seedling_sprsur)
plot(seedling_sprsur)
fixef(seedling_sprsur)
conditional_effects(seedling_sprsur)


### Final Params: seedling competitive effects
# value from mordecai is .5, oddly much higher than other competitive effects in the model...?
alpha_as=theta*alpha_aa
alpha_ps=theta*alpha_pa

theta~Uniform(0, 1)
