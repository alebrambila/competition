library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source('data.cleaning.R')

#ANNUALS
annuals2<-annuals%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         percap=as.integer(percap))%>%
  select(warmtrt, percap, starting_pm2, seeded_am2)

N <- length(annuals2$percap)
percap <- annuals2$percap
warmtrt <- annuals2$warmtrt
starting_pm2 <- annuals2$starting_pm2
seeded_am2 <- annuals2$seeded_am2

annuals_datavector <- c("N", "percap", "seeded_am2","starting_pm2", "warmtrt")

initials <- list(lambdaA_amb=4.5, lambdaA_slope=-0.5, 
                 alphaAA_amb=-1, alphaAA_slope=0, 
                 alphaAP_amb=-1, alphaAP_slope=0)

initials1<- list(initials, initials, initials)

annuals_fit <- stan(file="annuals_model.stan", 
                    data=annuals_datavector,
                    iter=5000,
                    chains = 3,
                   # thin = 1,
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                    init = initials1)

traceplot(annuals_fit, pars="lambdaA_amb")
pairs(annuals_fit, pars=c("lambdaA_amb", "lambdaA_slope", "alphaAA_amb", "alphaAA_slope", "alphaAP_amb", "alphaAP_slope"))
# pairs(no_dist_seeds_brho_hi_hi, pars = c('lambda_int', 'lambda_slope')

### Save posterior distributions to file
save(annuals_fit, file = "annuals_fit032122.rdata")
#annuals_fit<-load("annuals_fit032122.rdata")

## Look at resulting estimated parameter distributions
stan_dens(annuals_fit, pars = c("lambdaA_amb", "lambdaA_slope", "alphaAA_amb", "alphaAA_slope", "alphaAP_amb", "alphaAP_slope"))


#ADULT PERENNIALS

adults<-dat_p%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         fecundity=as.integer(fecundity))%>%
  select(warmtrt, fecundity, starting_pm2, seeded_am2)

N <- length(adults$fecundity)
fecundity <- adults$fecundity
warmtrt <- adults$warmtrt
starting_pm2 <- adults$starting_pm2
seeded_am2 <- adults$seeded_am2

adults_datavector <- c("N", "fecundity", "seeded_am2","starting_pm2", "warmtrt")

initials <- list(lambdaP_amb=10, lambdaP_slope=-0.5, 
                 alphaPA_amb=-1, alphaPA_slope=0, 
                 alphaPP_amb=-1, alphaPP_slope=0)

initials1<- list(initials, initials, initials)

adults_fit <- stan(file="adults_model.stan", 
                    data=adults_datavector,
                    iter=3000,
                    chains = 3,
                    # thin = 1,
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                    init = initials1)

traceplot(adults_fit, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))
pairs(adults_fit, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))
# pairs(no_dist_seeds_brho_hi_hi, pars = c('lambda_int', 'lambda_slope')

### Save posterior distributions to file
save(adults_fit, file = "adults_fit032122.rdata")
#annuals_fit<-load("adults_fit032122.rdata")

## Look at resulting estimated parameter distributions
stan_dens(adults_fit, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))



#SEEDLING PERENNIALS

seedlings2<-seedlings%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         survival=as.integer((fall.g/seeded_sm2*10000)))%>% #survival is a proportion, but i need integers so scale up
  select(warmtrt, survival, starting_pm2, seeded_am2, seeded_sm2)

N <- length(seedlings2$survival)
survival <- seedlings2$survival
warmtrt <- seedlings2$warmtrt
starting_pm2 <- seedlings2$starting_pm2
seeded_am2 <- seedlings2$seeded_am2
seeded_sm2 <- seedlings2$seeded_sm2

seedlings_datavector <- c("N", "survival", "seeded_am2","seeded_sm2", "starting_pm2", "warmtrt")

initials <- list(lambdaS_amb=7, lambdaS_slope=0, 
                 alphaSA_amb=-1, alphaSA_slope=0, 
                 alphaSS_amb=-1, alphaSS_slope=0,
                 alphaSP_amb=-1, alphaSP_slope=0)

initials1<- list(initials, initials, initials)

seedlings_fit <- stan(file="seedlings_model.stan", 
                   data=seedlings_datavector,
                   iter=3000,
                   chains = 3,
                   # thin = 1,
                   control = list(adapt_delta = 0.9, max_treedepth = 10),
                   init = initials1)

traceplot(seedlings_fit, pars=c("lambdaS_amb", "lambdaS_slope", "alphaSA_amb", "alphaSA_slope", "alphaSP_amb", "alphaSP_slope", "alphaSS_amb", "alphaSS_slope"))
pairs(seedlings_fit, pars=c("lambdaS_amb", "lambdaS_slope", "alphaSA_amb", "alphaSA_slope", "alphaSP_amb", "alphaSP_slope", "alphaSS_amb", "alphaSS_slope"))
# pairs(no_dist_seeds_brho_hi_hi, pars = c('lambda_int', 'lambda_slope')

### Save posterior distributions to file
save(seedlings_fit, file = "seedlings_fit032122.rdata")
#annuals_fit<-load("adults_fit032122.rdata")

## Look at resulting estimated parameter distributions
stan_dens(seedlings_fit, pars=c("lambdaS_amb", "lambdaS_slope", "alphaSA_amb", "alphaSA_slope", "alphaSP_amb", "alphaSP_slope", "alphaSS_amb", "alphaSS_slope"))
#remember, survival (lambdaS) is scaled up by 10,000X so a value of exp(4.4)/10,000 is actually 0.8% survival, or ~ 1 in 100



## Extract all parameter estimates
annuals_estimates <- rstan::extract(annuals_fit)
adults_estimates <- rstan::extract(adults_fit)
seedlings_estimates <- rstan::extract(seedlings_fit)

