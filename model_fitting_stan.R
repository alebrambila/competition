library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source('data.cleaning.R')

#ANNUALS ----
annuals2020<-subset(annuals, time==2020)%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         percap=as.integer(percap))%>%
  select(warmtrt, percap, starting_pm2, seeded_am2)

annuals2021<-subset(annuals, time==2021)%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         percap=as.integer(percap))%>%
  select(warmtrt, percap, starting_pm2, seeded_am2)

ggplot(annuals2020, aes(y=seeded_am2, x=starting_pm2, color=log(percap))) +geom_jitter()
ggplot(annuals, aes(x=seeded_am2, color=starting_pm2, y=log(percap))) +geom_point()

# alternate versions with perennial phytometers dropped to zero in no competition treatment
# because they are far away

annuals2020<-subset(annuals, time==2020)%>%
  mutate(starting_pm2=ifelse(comptrt=="none"|comptrt=="seedling perennials", 0, starting_pm2 ))%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         percap=as.integer(percap))%>%
  select(warmtrt, percap, starting_pm2, seeded_am2)

annuals2021<-subset(annuals, time==2021)%>%
  mutate(starting_pm2=ifelse(comptrt=="none"|comptrt=="seedling perennials", 0, starting_pm2 ))%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         percap=as.integer(percap))%>%
  select(warmtrt, percap, starting_pm2, seeded_am2)


#check prior by finiding max and mean percap of the highest treatment
annuals_lambdaprior<-annuals%>%
  filter(starting_pm2<3&seeded_am2<100)%>%
  select(percap)%>%
  log()

#2020 data
N <- length(annuals2020$percap)
percap <- annuals2020$percap
warmtrt <- annuals2020$warmtrt
starting_pm2 <- annuals2020$starting_pm2
seeded_am2 <- annuals2020$seeded_am2
annuals_datavector <- c("N", "percap", "seeded_am2","starting_pm2", "warmtrt")



initials <- list(lambdaA_amb=4.5, lambdaA_slope=-0.5, 
                 alphaAA_amb=-1, alphaAA_slope=0, 
                 alphaAP_amb=-1, alphaAP_slope=0)

initials1<- list(initials, initials, initials)

annuals_fit2020 <- stan(file="annuals_model.stan", 
                    data=annuals_datavector,
                    iter=5000,
                    chains = 3,
                   # thin = 1,
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                    init = initials1)


#2021data
N <- length(annuals2021$percap)
percap <- annuals2021$percap
warmtrt <- annuals2021$warmtrt
starting_pm2 <- annuals2021$starting_pm2
seeded_am2 <- annuals2021$seeded_am2
annuals_datavector <- c("N", "percap", "seeded_am2","starting_pm2", "warmtrt")

annuals_fit2021 <- stan(file="annuals_model.stan", 
                        data=annuals_datavector,
                        iter=5000,
                        chains = 3,
                        # thin = 1,
                        control = list(adapt_delta = 0.9, max_treedepth = 10),
                        init = initials1)

#traceplot(annuals_fit, pars="lambdaA_amb")
#pairs(annuals_fit, pars=c("lambdaA_amb", "lambdaA_slope", "alphaAA_amb", "alphaAA_slope", "alphaAP_amb", "alphaAP_slope"))
# pairs(no_dist_seeds_brho_hi_hi, pars = c('lambda_int', 'lambda_slope')

### Save posterior distributions to file
#save(annuals_fit, file = "annuals_fit032122.rdata")
#annuals_fit<-load("annuals_fit032122.rdata")

## Look at resulting estimated parameter distributions
stan_dens(annuals_fit2020, pars = c("lambdaA_amb", "lambdaA_slope", "alphaAA_amb", "alphaAA_slope", "alphaAP_amb", "alphaAP_slope"))
stan_dens(annuals_fit2021, pars = c("lambdaA_amb", "lambdaA_slope", "alphaAA_amb", "alphaAA_slope", "alphaAP_amb", "alphaAP_slope"))

annuals_estimates2020 <-as.data.frame(annuals_fit2020)%>%
  mutate(lambdaA_warm=lambdaA_amb+lambdaA_slope, 
         alphaAA_warm=alphaAA_amb+alphaAA_slope,
         alphaAP_warm=alphaAP_amb+alphaAP_slope)%>%
  select(-lambdaA_slope, -alphaAA_slope, -alphaAP_slope)%>%
  gather(param, value, 1:7)%>%
  mutate(exp=exp(value))%>%
  group_by(param)%>%
  summarize(raw=mean(value), exp=mean(exp))

annuals_estimates2021 <-as.data.frame(annuals_fit2021)%>%
  mutate(lambdaA_warm=lambdaA_amb+lambdaA_slope, 
         alphaAA_warm=alphaAA_amb+alphaAA_slope,
         alphaAP_warm=alphaAP_amb+alphaAP_slope)%>%
  select(-lambdaA_slope, -alphaAA_slope, -alphaAP_slope)%>%
  gather(param, value, 1:7)%>%
  mutate(exp=exp(value))%>%
  group_by(param)%>%
  summarize(raw=mean(value), exp=mean(exp))


## plot against data
post<-as.data.frame(get_posterior_mean(annuals_fit2021))

bev1 <- function(x, lam, alpha_a) {(lam / (1+alpha_a*x))} 
bev1.exp <- function(x, lam, alpha) {exp(lam / (1+alpha*x))} 

#against pm2
ggplot(subset(annuals, time==2020))+
  geom_jitter(data=filter(annuals, warmtrt == "amb"), aes(x=starting_pm2, y=percap), color="blue1", width=.5)+
 geom_function(fun=bev1, color='blue1', size=1.5, 
               args=c("lam"=exp(post$`mean-chain:1`[1]), "alpha_a"=exp(post$`mean-chain:1`[5])),linetype='dashed') +
  geom_jitter(data=filter(annuals, warmtrt == "warm"), aes(x=starting_pm2, y=percap), color="red1", width=.5) +
  geom_function(fun=bev1, color='red1', size=1.5, 
                args=c("lam"=exp(post$`mean-chain:1`[1]+post$`mean-chain:1`[2]), "alpha_a"=exp(post$`mean-chain:1`[5]+post$`mean-chain:1`[2])),linetype='dashed')# +
  geom_function(fun=bev1.exp, color='blue1', size=1.5, 
              args=c("lam"=(post$`mean-chain:1`[1]), "alpha"=(post$`mean-chain:1`[5])))+
  geom_function(fun=bev1.exp, color='red1', size=1.5, 
              args=c("lam"=(post$`mean-chain:1`[1]+post$`mean-chain:1`[2]), "alpha"=(post$`mean-chain:1`[5]+post$`mean-chain:1`[2])))




#ADULT PERENNIALS ----

#check prior by finiding max and mean percap of the highest treatment
adults_lambdaprior<-adults%>%
  filter(starting_pm2<3&seeded_am2<100)%>%
  select(percap)%>%
  log()


adults<-dat_p%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         fecundity=as.integer(fecundity))%>%
  select(warmtrt, fecundity, starting_pm2, seeded_am2, time, comptrt)

#check prior by finiding max and mean percap of the highest treatment
ggplot(adults, aes(x=starting_pm2, y=seeded_am2, color=log(fecundity)))+geom_jitter()

adults_lambdaprior<-adults%>%
  filter(starting_pm2<3&seeded_am2<100)%>%
  select(fecundity)%>%
  log()

adults2020<-filter(adults, time==2020)
adults2021<-filter(adults, time==2021)

# alternate versions with perennial phytometers dropped to zero in no competition treatment
# because they are far away

adults2020<-filter(adults, time==2020)%>%
  mutate(seeded_am2=ifelse(comptrt=="none"|comptrt=="seedling perennials", 0, seeded_am2 ))%>%
  select(warmtrt, fecundity, starting_pm2, seeded_am2)

adults2021<-filter(adults, time==2021)%>%
  mutate(seeded_am2=ifelse(comptrt=="none"|comptrt=="seedling perennials", 0, seeded_am2 ))%>%
  select(warmtrt, fecundity, starting_pm2, seeded_am2)



#2020
N <- length(adults2020$fecundity)
fecundity <- adults2020$fecundity
warmtrt <- adults2020$warmtrt
starting_pm2 <- adults2020$starting_pm2
seeded_am2 <- adults2020$seeded_am2

#2021
N <- length(adults2021$fecundity)
fecundity <- adults2021$fecundity
warmtrt <- adults2021$warmtrt
starting_pm2 <- adults2021$starting_pm2
seeded_am2 <- adults2021$seeded_am2

adults_datavector <- c("N", "fecundity", "seeded_am2","starting_pm2", "warmtrt")

initials <- list(lambdaP_amb=10, lambdaP_slope=-0.5, 
                 alphaPA_amb=-1, alphaPA_slope=0, 
                 alphaPP_amb=-1, alphaPP_slope=0)

initials1<- list(initials, initials, initials)

adults_fit2021 <- stan(file="adults_model.stan", 
                    data=adults_datavector,
                    iter=3000,
                    chains = 3,
                    # thin = 1,
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                    init = initials1)

traceplot(adults_fit2020, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))
traceplot(adults_fit2021, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))

pairs(adults_fit2020, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))
pairs(adults_fit2021, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))

### Save posterior distributions to file
save(adults_fit, file = "adults_fit032122.rdata")
#annuals_fit<-load("adults_fit032122.rdata")

## Look at resulting estimated parameter distributions
stan_dens(adults_fit2020, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))
stan_dens(adults_fit2021, pars=c("lambdaP_amb", "lambdaP_slope", "alphaPA_amb", "alphaPA_slope", "alphaPP_amb", "alphaPP_slope"))


adults_estimates2020 <-as.data.frame(adults_fit2020)%>%
  mutate(lambdaP_warm=lambdaP_amb+lambdaP_slope, 
         alphaPA_warm=alphaPA_amb+alphaPA_slope,
         alphaPP_warm=alphaPP_amb+alphaPP_slope)%>%
  select(-lambdaP_slope, -alphaPA_slope, -alphaPP_slope)%>%
  gather(param, value, 1:7)%>%
  mutate(exp=exp(value))%>%
  group_by(param)%>%
  summarize(raw=mean(value), exp=mean(exp))

adults_estimates2021 <-as.data.frame(adults_fit2021)%>%
  mutate(lambdaP_warm=lambdaP_amb+lambdaP_slope, 
         alphaPA_warm=alphaPA_amb+alphaPA_slope,
         alphaPP_warm=alphaPP_amb+alphaPP_slope)%>%
  select(-lambdaP_slope, -alphaPA_slope, -alphaPP_slope)%>%
  gather(param, value, 1:7)%>%
  mutate(exp=exp(value))%>%
  group_by(param)%>%
  summarize(raw=mean(value), exp=mean(exp))


#SEEDLING PERENNIALS

seedlings2020<-filter(seedlings, time==2020)%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         survivors=as.integer(fall.g)) %>%#make integer
  select(warmtrt, survivors, starting_pm2, seeded_am2, seeded_sm2)

seedlings2021<-filter(seedlings, time==2021)%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         survivors=as.integer(fall.g)) %>%#make integer
  select(warmtrt, survivors, starting_pm2, seeded_am2, seeded_sm2)

#alternate modified
seedlings2020<-filter(seedlings, time==2020)%>%
  mutate(seeded_am2=ifelse(comptrt=="none"|comptrt=="adult perennials", 0, seeded_am2 ))%>%
  mutate(starting_pm2=ifelse(comptrt=="none"|comptrt=="annuals", 0, starting_pm2 ))%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         survivors=as.integer(fall.g)) %>%#make integer
  select(warmtrt, survivors, starting_pm2, seeded_am2, seeded_sm2)

seedlings2021<-filter(seedlings, time==2021)%>%
  mutate(seeded_am2=ifelse(comptrt=="none"|comptrt=="adult perennials", 0, seeded_am2 ))%>%
  mutate(starting_pm2=ifelse(comptrt=="none"|comptrt=="annuals", 0, starting_pm2 ))%>%
  mutate(warmtrt = ifelse(warmtrt == 'warm', 1, 0), 
         survivors=as.integer(fall.g)) %>%#make integer
  select(warmtrt, survivors, starting_pm2, seeded_am2, seeded_sm2)


#what is a reasonable prior
#check prior by finiding max and mean percap of the highest treatment
annuals_lambdaprior<-seedlings%>%
  filter(starting_pm2<3&seeded_am2<100)%>%
  mutate(survival=fall.g/seeded_sm2)%>%
  select(survival)%>%
  log()
#mean: -4.7
#max: -3.6



#2020
N <- length(seedlings2020$survivors)
survivors <- seedlings2020$survivors
warmtrt <- seedlings2020$warmtrt
starting_pm2 <- seedlings2020$starting_pm2
seeded_am2 <- seedlings2020$seeded_am2
seeded_sm2 <- seedlings2020$seeded_sm2
#2021
N <- length(seedlings2021$survivors)
survivors <- seedlings2021$survivors
warmtrt <- seedlings2021$warmtrt
starting_pm2 <- seedlings2021$starting_pm2
seeded_am2 <- seedlings2021$seeded_am2
seeded_sm2 <- seedlings2021$seeded_sm2

seedlings_datavector <- c("N", "survivors", "seeded_am2","seeded_sm2", "starting_pm2", "warmtrt")

initials <- list(lambdaS_amb=7, lambdaS_slope=0, 
                 alphaSA_amb=-1, alphaSA_slope=0, 
                 alphaSS_amb=-1, alphaSS_slope=0,
                 alphaSP_amb=-1, alphaSP_slope=0)

initials1<- list(initials, initials, initials)

seedlings_fit2021 <- stan(file="seedlings_model.stan", 
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
stan_dens(seedlings_fit2020, pars=c("survivalS_amb", "survivalS_slope", "alphaSA_amb", "alphaSA_slope", "alphaSP_amb", "alphaSP_slope", "alphaSS_amb", "alphaSS_slope"))
stan_dens(seedlings_fit2021, pars=c("survivalS_amb", "survivalS_slope", "alphaSA_amb", "alphaSA_slope", "alphaSP_amb", "alphaSP_slope", "alphaSS_amb", "alphaSS_slope"))

seedlings_estimates2020 <-as.data.frame(seedlings_fit2020)%>%
  mutate(survivalS_warm=survivalS_amb+survivalS_slope, 
         alphaSA_warm=alphaSA_amb+alphaSA_slope,
         alphaSP_warm=alphaSP_amb+alphaSP_slope,
         alphaSS_warm=alphaSS_amb+alphaSS_slope)%>%
  select(-survivalS_slope, -alphaSA_slope, -alphaSP_slope, -alphaSS_slope)%>%
  gather(param, value, 1:9)%>%
  mutate(exp=exp(value))%>%
  group_by(param)%>%
  summarize(raw=mean(value), exp=mean(exp))%>%
  filter(param!="lp__")

seedlings_estimates2021 <-as.data.frame(seedlings_fit2021)%>%
  mutate(survivalS_warm=survivalS_amb+survivalS_slope, 
         alphaSA_warm=alphaSA_amb+alphaSA_slope,
         alphaSP_warm=alphaSP_amb+alphaSP_slope,
         alphaSS_warm=alphaSS_amb+alphaSS_slope)%>%
  select(-survivalS_slope, -alphaSA_slope, -alphaSP_slope, -alphaSS_slope)%>%
  gather(param, value, 1:9)%>%
  mutate(exp=exp(value))%>%
  group_by(param)%>%
  summarize(raw=mean(value), exp=mean(exp))%>%
  filter(param!="lp__")

## Extract all parameter estimates
annuals_estimates <- rstan::extract(annuals_fit)
adults_estimates <- rstan::extract(adults_fit)
seedlings_estimates <- rstan::extract(seedlings_fit)

