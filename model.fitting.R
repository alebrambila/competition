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
# ~~ I haven't actually made warmed and unwarmed models yet, just one version for each 'species' ~~

## Questions:
# How do I take into account random effects in the models: i.e. block, or eventually, year. 
# What is the second function in 'formulas'?  How do I incorporate it?
# How to do the binomial functions for survival.
# Are counts or cover more appropriate for competition context.  
  #"Per Capita' implies counts, but with variation in size, cover might be a better representation of competition strength.
  # In that case the alphas wouldn't be per-capita, but per-unit of cover (%)
# What is my unit of replication.  If it's plots in each treatment its very low (n=4), 
  #If it's individuals that survive/reproduce/etc. it's much more, but its generally not how my data is set up
  # I realized I had this question because I have individual fecundity/survival data for phytometers, but also 
      # a bunch more that I can infer from counts and plot level seed production data.  Can I merge these?

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
density_spring21 <- ### TO BE GATHERED IN JUNE 2021
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
sumsur2021s <- ### TO BE GATHERED JUNE AND SEPTEMBER 2021
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

### Annual Parameter Fitting ----

### Estimating Lambdas and Alphas
# out_a = measured seed output (per capita):  DATA!
# lambda_a = annual per capita seed production in the absence of competition: PARAMETER!
# alpha_aa = per capita competitive effect of annuals on annuals: PARAMETER!
# density_a = density of annuals (cover): DATA!
#alpha_ap = per capita competitive effect of adult perennials on annuals: PARAMETER
# density_p = density of adult perennials: DATA!

#There is a second function that I'm not sure what to do with: Y ~ Normal(log(out_a)) where Y is the actual per capita seed output.  
  # there is also a parameter for error term for this formula: 1/t, t~1/Gamma(.001, .001)


dat<-left_join(fecundity_phyt_a, density_spring20)%>%
  mutate(out_a=seeds, density_a=am2, density_p=pm2)%>%
  select(1, 6, 7, 8, out_a, density_a, density_p)

# quick graphs of variables going into model
ggplot(dat, aes(y=out_a)) +
  geom_jitter(aes(x=density_a, color=warmtrt))
ggplot(dat, aes(y=out_a)) +
  geom_jitter(aes(x=density_p, color=warmtrt))
ggplot(dat, aes(x=density_a, y=density_p)) +
  geom_jitter(aes(color=warmtrt))


ggplot(dat, aes(x=density_a, y=out_a, color=warmtrt)) +
  geom_jitter()+
  # geom_smooth(method = 'lm',formula = y ~ x + I(x^2))
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 300, b = .1)),
              se = FALSE)

# check replication within blocks
table(dat$warmtrt, dat$block)

# BRMS Version of Lambda/Alpha
## NOTES / CHANGES: 
# parameter names can't have underscore or dots
# it's not liking the gamma priors (not clear to me why), so changed to normal

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
                     control = list(adapt_delta = 0.99, max_treedepth = 18))


### Estimating annual germination and spring survival (NO COMPETITIVE EFFECT?)

dat<-sprsur2020a<-filter(spr_sur2020, seeded_a!=0)%>%
  dplyr::select(-seeded_s, -spring20_s)%>%
  mutate(sprsur_a=spring20_a/seeded_a)


# setting up beta-binomial (vs regular binomial) as described here: https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html#the-beta-binomial-distribution
# need a custom family to account for overdispersion
beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"

stanvars <- stanvar(scode = stan_funs, block = "functions")


annual_sprsur <- brm(bf(spring20_a|vint(seeded_a) ~ warmtrt + (1|block)), 
    data = dat,
    family=beta_binomial2,
    prior = c(prior(gamma(1, 1), lb=0, nlpar = "a"), # not nlpar, not a nonlinear
              prior(gamma(1, 1), nlpar = "b"),
    inits = "0",  #list(lambda=100, aA=1, aB=1, aL=1, aV=1, aE=1),
    cores=4, 
    chains=4,
    iter=10000, 
    thin=2,
    control = list(adapt_delta = 0.99, max_treedepth = 18))



### Adult Perennial Parameter Fitting ----
### Adult Lambda/Alphas

# out_p = measured adult seed output (per capita):  DATA!
# lambda_p = adult per capita seed production in the absence of competition: PARAMETER!
# alpha_pa = per capita competitive effect of annuals on adults: PARAMETER!
# density_a = density of annuals (cover): DATA!
#alpha_pp = per capita competitive effect of adult perennials on perennials: PARAMETER
# density_p = density of adult perennials: DATA!

#Same second function: Y ~ Normal(log(out_p)) where Y is the actual per capita seed output, with: sigma=1/t, t~1/Gamma(.001, .001)

#Least Squares Version
bevholt_adultlambda <- as.formula(out_p=lambda_p/(1+alpha_pa*density_a + 
                                                     alpha_pp*density_p))

m1_adultlambda <- nlsLM(bev_holt, start = list(lambda_a=200, alpha_pa=5, alpha_pp=5),
                         lower = c(0, -Inf, -Inf), 
                         control=nls.lm.control(maxiter=40000), trace=T,
                         data = dat)


# BRMS Version of Lambda/Alpha
m3_vumy <- brm(bf(out_p=lambda_p/(1+alpha_pa*density_a + 
                                    alpha_pp*density_p), 
                  lambda_p ~ 1,
                  alpha_pa ~ 1,
                  alpha_pp ~ 1,
                  nl=TRUE),
               data = dat,
               prior = c(prior(gamma(200, 1), lb=0, nlpar = "lambda_p"), 
                         prior(gamma(1, 9), nlpar = "alpha_pa"),
                         prior(gamma(1, 9), nlpar = "alpha_pp")),
               inits = "0",  #list(lambda=100, aA=1, aB=1, aL=1, aV=1, aE=1),
               cores=4, 
               chains=4,
               iter=10000, 
               thin=2,
               control = list(adapt_delta = 0.99, max_treedepth = 18))


### Adult summer survival
summer_p ~ Beta-Binomial( a, b, spring_p)
# mean value from mordecai is .88, i had 100%
# priors a~Gamma(8, 1) b~Gamma(.01, 1)


### Seedling Perennial Parameter Fitting ----

# out_s = measured perennial seedling survival:  DATA! (mean=.29 for mordecai)
# summersurvival_s = seedling per capita summer survival in the absence of competition: PARAMETER!
# alpha_sa = per capita competitive effect of annuals on adults: PARAMETER!
# density_a = density of annuals (cover): DATA!
#alpha_sp = per capita competitive effect of adult perennials on seedlings: PARAMETER
# density_p = density of adult perennials: DATA!
#alpha_ss = per capita competitive effect of seedlings on seedlings: PARAMETER
# density_s = density of seedlings: DATA!

#DIFFERENT second function: Y ~ Normal(log(out_p)) where Y is the actual # of surviving seedlings, 
    # with: Y ~ Beta-Binomial(a, a(1-out_s)/out_s, density_x) a~Gamma(1, 1)

#Least Squares Version
bevholt_seedling_summer <- as.formula(out_s=summersurvival_s/(1+alpha_sa*density_a + 
                                                                alpha_ss*density_s+
                                                                alpha_sp*density_p))

m1_summersurvival_s <- nlsLM(bev_holt, start = list(summersurvival_s=0.3, alpha_sa=5, alpha_ss=5, alpha_sp=5),
                        lower = c(0, -Inf, -Inf, -Inf), 
                        control=nls.lm.control(maxiter=40000), trace=T,
                        data = dat)


# BRMS Version of Seedling Summer Survival
m3_vumy <- brm(bf(out_s=summersurvival_s/(1+alpha_sa*density_a + 
                                            alpha_ss*density_s+
                                            alpha_sp*density_p), 
                  summersurvival_s ~ 1,
                  alpha_sa ~ 1,
                  alpha_ss ~ 1,
                  alpha_sp ~ 1,
                  nl=TRUE),
               data = dat,
               prior = c(prior(beta(1, 1), lb=0, nlpar = "summersurvival_s"), 
                         prior(gamma(1, 9), nlpar = "alpha_sa"),
                         prior(gamma(1, 9), nlpar = "alpha_ss"),
                         prior(gamma(1, 9), nlpar = "alpha_sp")),
               inits = "0",  #list(lambda=100, aA=1, aB=1, aL=1, aV=1, aE=1),
               cores=4, 
               chains=4,
               iter=10000, 
               thin=2,
               control = list(adapt_delta = 0.99, max_treedepth = 18))


### Seedling Spring Survival              - is .34 in mordecai
bevholt_seedlingspring <- as.formula(spring_s~BetaBinomial(a, b, seedsin_s)) # not sure how to do this binomial part
priors: a, b, Gamma(1, 1)

### Final Params: seedling competitive effects
# value from mordecai is .5, oddly much higher than other competitive effects in the model...?
alpha_as=theta*alpha_aa
alpha_ps=theta*alpha_pa

theta~Uniform(0, 1)
