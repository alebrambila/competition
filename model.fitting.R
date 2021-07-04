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

# JD: look at data
head(annuals); dim(annuals)
# how many plots / block?  answer: 12
annuals %>% group_by(block) %>%
  summarize(n.plots = n_distinct(plotid))
table(annuals$warmtrt, annuals$block)

# per capita
ggplot(subset(annuals), aes(x=pm2, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F) +xlab("perennial adults")+ylab("annual percapita seeds out")

ggplot(subset(annuals), aes(x=seeded_a, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+xlab("annual seeds in")+ylab("annual percapita seeds out")

#added in seedling perennials (not sure it's going to be helpful)
ggplot(subset(annuals, seeded_a>150), aes(x=seeded_s, y=percap, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt), width=.5)+scale_colour_manual(values = c("dodgerblue", "darkred"))+
  geom_smooth(method = 'lm',formula = y ~ x , se=F)+xlab("perennial seeds in")+ylab("annual percapita seeds out")

#ggplot(annuals, aes(x=pm2, y=seeded_a, color=percap)) +
#  geom_jitter(aes(shape=comptrt), height=50)+
#  scale_color_gradient(low = "blue", high = "red", na.value = NA)

### BRM fits ----
## 1.1 Simple (seeds in:out) annuals percapita and scaled
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

################################
### 2. ADULT PERENNIAL FECUNDITY ----
#data ----
datp<-left_join(fecundity_phyt_p, density_spring20)%>%
  mutate(out_p=seeds, density_a=am2, density_p=pm2, density_s=sm2)%>%
  select(1, 6, 7, 8, out_p, density_a, density_p, density_s)
datp<-left_join(datp, dplyr::select(annuals, plotid, seeded_a))%>%
  mutate(seeded_a=ifelse(is.na(seeded_a), 8, seeded_a))

#visualizations----
ggplot(datp, aes(x=density_p, y=out_p, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+
  labs(x="adult perennial density", y="adult perennial fecundity")

ggplot(datp, aes(x=seeded_a, y=out_p, color=warmtrt)) +
  geom_jitter(aes(shape=comptrt))+
  geom_smooth(method = 'lm',formula = y ~ x + I(x^2), se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+
  labs(x="annual seeds in", y="adult perennial fecundity")

table(datp$warmtrt, datp$block)
table(datp$warmtrt, datp$comptrt)

datp$warmtrt <- as.factor(datp$warmtrt)

#Fit BRM----
## try simple model first to get ballpark param estimates
#perennial_lambda0 <- brm(bf(out_p ~ lambdaP*5000 / (1+alphaPA*seeded_a + alphaPP*density_p), 
#                           lambdaP ~ 1, 
#                           alphaPA ~  1, 
#                           alphaPP ~  1, nl=TRUE),
#                        data = subset(datp),
#                        prior = c(prior(normal(1, 1), lb=0, nlpar = "lambdaP"), 
#                                  prior(normal(0, 1), nlpar = "alphaPA"), # broadedned the priors
#                                  prior(normal(0, 1), nlpar = "alphaPP")),
#                        inits = "0",  
#                        cores=3, 
#                        chains=3,
#                        #refresh=100,
#                        iter=5000, 
#                        thin=2,
#                        control = list(adapt_delta = 0.99, max_treedepth = 18))
#perennial_lambda0
#plot(perennial_lambda0)
#conditional_effects(perennial_lambda0)

### 2.1 Perennial fecundity scaled
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
#plw1<-perennial.lambda
#stancode(perennial.lambda)
#perennial.lambda
#plot(perennial.lambda)
#fixef(perennial.lambda)
#conditional_effects(readPL)%>% plot(points=T)
#conditional_effects(readPL, effects = "seeded_a:warmtrt")%>% plot(points=T)
#conditional_effects(readPL, effects = "density_p:warmtrt")%>% plot(points=T)

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
#data ----

dat.sprsurv.0 <- left_join(plotkey, spr_sur2020)%>%  # just naming the data something different 
  dplyr::select(-seeded_a, -spring20_a)
dat.sprsurv.00<-left_join(dat.sprsurv.0, dplyr::select(annuals, plotid, seeded_a))%>%
  mutate(seeded_a=ifelse(comptrt=="none"|comptrt=='seedling perennials', 4/0.06019467, ifelse(comptrt=="adult perennials", 8/0.06019467, ifelse(comptrt=="seedlings+adults", 0, seeded_a))))%>% # germination correction factor for phytometers.  how many stems I had in spring divided by germination factor to give me how many seeds were added.   
    mutate(seeded_s=ifelse(comptrt=="none"|comptrt=='annuals', 4/0.06019467, ifelse(comptrt=="adult perennials", 8/0.06019467, ifelse(comptrt=="annuals+adults", 0, seeded_s))))%>% # germination correction factor for phytometers.  how many stems I had in spring divided by germination factor to give me how many seeds were added.   
      mutate(seeded_s=ifelse(comptrt=='seedling perennials', 4500, ifelse(comptrt=="seedlings+adults", 1500, seeded_s))) %>% # germination correction factor for phytometers.  how many stems I had in spring divided by germination factor to give me how many seeds were added.   
  dplyr::select(1, 2, 3, 4, 9, 10, 11)
dat.sprsurv.000<-left_join(dat.sprsurv.00, select(datp, plotid, density_p))%>%
  mutate(density_p=ifelse(is.na(density_p), 0, density_p))
dat.sprsurv<-left_join(dat.sprsurv.000, select(density_spring20, plotid, count_s))%>%
  mutate(spring_20s=count_s)%>%
  select(-count_s)
dat.sprsurv<-unique(dat.sprsurv)

dat.sumsurv<-left_join(select(dat.sprsurv, -spring20_s), select(seedling_sumsur2020, plotid, fall20_s, fall20_s.g, spring20_s,spring20_s.g, gopher_ss_correction))%>%
  mutate(fall20_s.g=as.integer(fall20_s.g), spring20_s.g=as.integer(spring20_s.g))%>%
  #filter(!is.na(spring20_s))%>%
  mutate(seeded_s.g=seeded_s*(1-gopher_ss_correction), seeded_a.g=seeded_a*(1-gopher_ss_correction)) # adjust seed addition for down for gophers 
dat.sumsurv<-unique(dat.sumsurv)

#vizualisation ----
#three way competition for seeds:adults, and subdivided into spring survival, and summer survival
# to see only plots where seedlings made up a significant component of the community I replaced
# data with subset(data, seeded_s>150).  this can be removed to see full data but proportions may be misleading

#seeds:adults vs. seeds
ggplot(subset(dat.sumsurv, seeded_s>100), aes(x=seeded_s.g, y=fall20_s/seeded_s.g, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("seeded perennials")+ylab("seed:adult survival rate")

#seeds:adults vs. perennials
ggplot(subset(dat.sumsurv, seeded_s>100), aes(x=density_p, y=fall20_s/seeded_s.g, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("adult perennial density")+ylab("seed:adult survival rate")

#seeds:adults vs. annuals
# cool - it looks like perennials facilitate themselves, but annuals inhibit
ggplot(subset(dat.sumsurv, seeded_s>150), aes(x=seeded_a, y=fall20_s/seeded_s.g, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("seeded annuals")+ylab("seed:adult survival rate")


#spring survival vs. seeds
ggplot(subset(dat.sprsurv,seeded_s!=0&seeded_s>150), aes(x=seeded_s, y=spring_20s/seeded_s, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("seeded perennials")+ylab("spring survival rate")

#spring survival vs. perennials
ggplot(subset(dat.sprsurv,seeded_s!=0&seeded_s>150), aes(x=density_p, y=spring_20s/seeded_s, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("adult perennials")+ylab("spring survival rate")

#spring survival vs. annnuals
ggplot(subset(dat.sprsurv,seeded_s!=0&seeded_s>150), aes(x=seeded_a, y=spring_20s/seeded_s, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("annual seeds")+ylab("spring survival rate")

#summer survival vs. self
ggplot(subset(dat.sumsurv, seeded_s>150), aes(x=spring20_s.g, y=fall20_s/spring20_s.g, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("perennial stems")+ylab("summer survival rate")

#summer survival vs. adults perennials
ggplot(subset(dat.sumsurv, seeded_s>150), aes(x=density_p, y=fall20_s/spring20_s.g, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("adult perennial density")+ylab("summer survival rate")

#summer survival vs. annuals
ggplot(subset(dat.sumsurv, seeded_s>150), aes(x=seeded_a, y=fall20_s/spring20_s.g, color=warmtrt)) +
  geom_point(aes(shape=comptrt))+
  geom_smooth(method = 'lm', se=F)+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+xlab("seeded annuals")+ylab("sumer survival rate")



### 3.1 Simple (seeds in:adults out) perennial seedlings binomial ----
seedling.binomial.a<- brm(bf(fall20_s|trials(seeded_s.g) ~ lambdaS / (1+alphaSA*seeded_a + alphaSS*seeded_s + alphaSP*density_p), 
                           lambdaS +alphaSA +alphaSP+alphaSS~ 1+ (1|block), nl=TRUE),
                        family=binomial,
                        data = subset(dat.sumsurv, seeded_s>150&warmtrt=="amb"),   #running this with limited dataset as in teh figures above (only in seedling comptrts)
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

