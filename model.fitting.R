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


#############################
##### ORIGINAL ANNUAL LAMBDA ----
dat<-left_join(fecundity_phyt_a, density_spring20)%>%
  mutate(out_a=seeds, density_a=am2, density_p=pm2)%>%
  select(1, 6, 7, 8, out_a, density_a, density_p)

###visualizations ----
ggplot(dat, aes(x=density_p, y=out_a, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=1)+
   geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  #stat_smooth(method = "nls",
       #       formula = y ~ a/(1+b*x),
      #       method.args = list(start = list(a = 100, b = .1)),
        #      se = FALSE)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("Na(t+1)")+xlab("Np(t))")

ggplot(dat, aes(x=density_a, y=out_a, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=1)+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("Na(t+1)")+xlab("Na(t))")

# check replication within blocks
table(dat$warmtrt, dat$block)

## BRM fit ----
annual_lambda <- brm(bf(out_a ~ lambdaA / (1+alphaAA*density_a + alphaAP*density_p), 
                        lambdaA ~ warmtrt + (1|block), 
                        alphaAA ~  warmtrt + (1|block), 
                        alphaAP ~  warmtrt + (1|block), nl=TRUE),
                     data = dat,
                     prior = c(prior(normal(200, 50), lb=0, nlpar = "lambdaA"), 
                               # prior(gamma(3, .01), nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaAA"),
                               prior(normal(0, .1), nlpar = "alphaAP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=10000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 16))
annual_lambda
plot(annual_lambda)
fixef(annual_lambda)
conditional_effects(annual_lambda)


###################################
##### ORIGINALANNUAL SPRING SURVIVAL ----
dat.surv.a <- sprsur2020a<-filter(spr_sur2020, seeded_a!=0)%>%  # just naming the data something different 
  dplyr::select(-seeded_s, -spring20_s)%>%
  mutate(sprsur_a=spring20_a/seeded_a)
# quick graphs of variables going into model ----
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

nd <- tibble(warmtrt = c("amb","warm"), block=NA, seeded_a=1)
n_iter<-100

fitted(annual_sprsur,newdata = nd, summary = TRUE ) %>%
  as_tibble() 

# Graph probabilities ~ treatment ----

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



###############################
## SIMPLIFIED ANNUALS MODEL ----
## data ----
annuals<-left_join(fecundity_plot_a, dplyr::select(dat.surv.a, 1,9))
annuals<-left_join(annuals, plotkey)%>%
  mutate(seeded_a=ifelse(comptrt=="none"|comptrt=='seedling perennials', 8, ifelse(comptrt=="adult perennials", 16, seeded_a)))%>%
  mutate(plotseeds=seeds*count_a)
annuals<-left_join(annuals, density_spring20)

#visualize ----
ggplot(annuals, aes(x=seeded_a, y=log(plotseeds), color=warmtrt)) +
  geom_point(aes(shape=comptrt), width=1)+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
 scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("Na(t+1)")+xlab("Na(t))")

### BRM fits ----
  annual.simple <- brm(bf(plotseeds ~ (lambdaA*seeded_a) / (1+alphaAA*seeded_a + alphaAP*pm2), 
                          lambdaA ~ warmtrt + (1|block), 
                          alphaAA ~  warmtrt + (1|block), 
                          alphaAP ~  warmtrt + (1|block), nl=TRUE),
                       data = annuals,
                       prior = c(prior(normal(800, 200), lb=0, nlpar = "lambdaA"), 
                                 # prior(gamma(3, .01), nlpar = "lambdaA"), 
                                 prior(normal(0, .1), nlpar = "alphaAA"),
                                 prior(normal(0, .1), nlpar = "alphaAP")),
                       inits = "0",  
                       cores=4, 
                       chains=4,
                       iter=10000, 
                       thin=5,
                       control = list(adapt_delta = 0.99, max_treedepth = 20))

plot(annual.simple)

### ambient only
annual.ambient <- brm(bf(plotseeds ~ (lambdaA*seeded_a) / (1+alphaAA*seeded_a + alphaAP*pm2), 
                        lambdaA ~ (1|block),
                        alphaAA ~  (1|block), 
                        alphaAP ~  (1|block), nl=TRUE),
                     data = subset(annuals, warmtrt=="amb"),
                     prior = c(prior(normal(800, 200), lb=0, nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaAA"),
                               prior(normal(0, .1), nlpar = "alphaAP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=10000, 
                     control = list(adapt_delta = 0.99, max_treedepth = 18))

### ambient only logged
# example from hallett: 
# m1A <- as.formula(log(AVseedout +1) ~ log(ag*(AVseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))

annual.ambient <- brm(bf(log(1+plotseeds) ~ log((seeded_a+1)*exp(log(lambdaA) - log((1+alphaAA*(seeded_a+1) + alphaAP*(pm2+1))))), 
                         lambdaA ~ (1|block),
                         alphaAA ~  (1|block), 
                         alphaAP ~  (1|block), nl=TRUE),
                      data = subset(annuals, warmtrt=="amb"),
                      prior = c(prior(normal(8, 2), lb=0, nlpar = "lambdaA"), 
                                prior(normal(0, .1), nlpar = "alphaAA"),
                                prior(normal(0, .1), nlpar = "alphaAP")),
                      inits = "0",  
                      cores=4, 
                      chains=4,
                      iter=10000, 
                      control = list(adapt_delta = 0.99, max_treedepth = 18))
################################
### ADULT PERENNIAL FECUNDITY ----
datp<-left_join(fecundity_phyt_p, density_spring20)%>%
  mutate(out_p=seeds, density_a=am2, density_p=pm2)%>%
  select(1, 6, 7, 8, out_p, density_a, density_p)
datp<-left_join(datp, dplyr::select(annuals, plotid, seeded_a))%>%
  mutate(seeded_a=ifelse(is.na(seeded_a), 8, seeded_a))

#visualizations----
ggplot(datp, aes(x=density_p, y=out_p, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("Np(t)")+ylab("ns(t+1)")

ggplot(datp, aes(x=seeded_a, y=out_p, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("ns(t+1)")+xlab("Na(t)")

#Fit BRM ----
perennial_lambda <- brm(bf(out_p ~ lambdaP / (1+alphaPA*seeded_a + alphaPP*density_p), 
                        lambdaP ~ warmtrt + (1|block), 
                        alphaPA ~  warmtrt + (1|block), 
                        alphaPP ~  warmtrt + (1|block), nl=TRUE),
                     data = subset(datp, warmtrt=="amb"),
                     prior = c(prior(normal(5000, 500), lb=0, nlpar = "lambdaP"), 
                               # prior(gamma(3, .01), nlpar = "lambdaA"), 
                               prior(normal(0, .1), nlpar = "alphaPA"),
                               prior(normal(0, .1), nlpar = "alphaPP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=10000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 18))

savedPL<-perennial_lambda
perennial_lambda
plot(perennial_lambda)
fixef(perennial_lambda)
conditional_effects(perennial_lambda)
#########################################
### ORIGINAL SEEDLING SUMMER SURVIVAL ----
dats<-left_join(select(seedling_sumsur2020, -time), density_spring20)%>%
  mutate(density_a=am2, density_p=pm2, density_s=pm2)%>%
  select(1, 8, 9, 10, fall20_s.g, spring20_s, density_a, density_p, density_s)%>%
  mutate(fall20_s.g=as.integer(fall20_s.g))

#visualizations----
ggplot(dats, aes(x=density_a, y=fall20_s.g/spring20_s, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ xlab("Na(t)")+ylab("Ss")

ggplot(dats, aes(x=density_p, y=fall20_s.g/spring20_s, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)+
  scale_colour_manual(values = c("dodgerblue", "darkred")) + xlab("Np(t)")+ylab("Ss")

### BRM fits ----
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

plot(seedling_sumsur)
summary(seedling_sumsur)


#######################################
### ORIGINAL SEEDLING SPRING SURVIVAL - is .34 in mordecai -----
dat.surv <- sprsur2020s<-filter(spr_sur2020, seeded_s!=0)%>%  # just naming the data something different 
  dplyr::select(-seeded_a, -spring20_a)%>%
  mutate(sprsur_s=spring20_s/seeded_s)

dat.surv

# quick graphs of variables going into model
ggplot(dat.surv, aes(x=seeded_s, y=spring20_s)) +
  geom_jitter(aes(color=warmtrt), width=100)+
  geom_smooth(aes(color=warmtrt), method="lm", se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ xlab("ns(t)")+ylab("Ns(t)")

ggplot(dat.surv, aes(x=warmtrt, y=spring20_s/seeded_s)) +
  geom_boxplot()+
  # geom_point(color='blue') +
  geom_jitter(width=.1) +
  scale_colour_manual(values = c("dodgerblue", "darkred"))

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





#############################################
## SIMPLIFIED SEEDLING (SEEDS IN, ADULTS OUT) ----
seedlings<-left_join(dats, dplyr::select(dat.surv, plotid, seeded_s))%>%
  mutate(seeded_s=ifelse(comptrt=="none"|comptrt=='annuals', 8, ifelse(comptrt=="adult perennials", 16, seeded_s)))
seedlings<-left_join(seedlings, dplyr::select(annuals, plotid, seeded_a))%>%
  mutate(seeded_a=ifelse(is.na(seeded_a), 8, seeded_a))

#vizualisation ----
ggplot(seedlings, aes(x=seeded_s, y=fall20_s.g, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+ylab("Np(t+1)")+xlab("Ns(t))")


#BRM fits ----

#normal----
seedling.simple <- brm(bf(fall20_s.g ~ (lambdaS*seeded_s) / (1+alphaSA*seeded_a + alphaSP*density_p+ alphaSS*seeded_s), 
                        lambdaS ~ warmtrt + (1|block), 
                        alphaSA ~  warmtrt + (1|block), 
                        alphaSP ~  warmtrt + (1|block),
                        alphaSS ~  warmtrt + (1|block),
                        nl=TRUE),
                     data = seedlings,
                     prior = c(prior(normal(.001, .0002), lb=0, nlpar = "lambdaS"), 
                               prior(normal(0, .1), nlpar = "alphaSA"),
                               prior(normal(0, .1), nlpar = "alphaSS"),
                               prior(normal(0, .1), nlpar = "alphaSP")),
                     inits = "0",  
                     cores=4, 
                     chains=4,
                     iter=100000, 
                     thin=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 20))

plot(seedling.simple)

#binomial ambient only----
seedling.ambient <- brm(bf(fall20_s.g|trials(seeded_s) ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                          lambdaS ~ (1|block), 
                          alphaSA ~  (1|block), 
                          alphaSP ~  (1|block), 
                          alphaSS ~  (1|block), nl=TRUE),
                       family=binomial,
                       data = subset(seedlings, warmtrt=="amb"),
                       prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaS"), 
                                 prior(normal(0, .1), nlpar = "alphaSA"),
                                 prior(normal(0, .1), nlpar = "alphaSS"),
                                 prior(normal(0, .1), nlpar = "alphaSP")),
                       inits = "0",  
                       cores=4, 
                       chains=4,
                       iter=5000, 
                       thin=5,
                       control = list(adapt_delta = 0.99, max_treedepth = 18))

#binomial warmed only ----
seedling.ambient <- brm(bf(fall20_s.g|trials(seeded_s) ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                           lambdaS ~ (1|block), 
                           alphaSA ~  (1|block), 
                           alphaSP ~  (1|block), 
                           alphaSS ~  (1|block), nl=TRUE),
                        family=binomial,
                        data = subset(seedlings, warmtrt=="warm"),
                        prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaS"), 
                                  prior(normal(0, .1), nlpar = "alphaSA"),
                                  prior(normal(0, .1), nlpar = "alphaSS"),
                                  prior(normal(0, .1), nlpar = "alphaSP")),
                        inits = "0",  
                        cores=4, 
                        chains=4,
                        iter=5000, 
                        thin=5,
                        control = list(adapt_delta = 0.99, max_treedepth = 18))



### Final Params: seedling competitive effects ----
# value from mordecai is .5, oddly much higher than other competitive effects in the model...?
alpha_as=theta*alpha_aa
alpha_ps=theta*alpha_pa

theta~Uniform(0, 1)