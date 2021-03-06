library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

### ### ### ### ### ### ### ### 
## Bev-Holt ----
### ### ### ### ### ### ### ### 

bev1 <- function(x, lam, alpha) {lam / (1+alpha*x)} # Beverton-Holt

# Set arbitrary parameters
N<-100 # number of data points
lam.true<-20 # true lambda (fitness, no competition)
alpha.true<-.1 # competition coefficient

# simulate data
set.seed(101)
dat = tibble(
  density_data = runif(N, 0, 30) # competitor density
  ,fitness_sim = bev1(density_data, lam.true, alpha.true) # fitness, using BH
  ,fitness_data = fitness_sim + rnorm(N,0, 2) # add some noise
)
dat

# Simple Plot data
ggplot(dat, aes(x= density_data, y=fitness_data)) +
  geom_point() +
  geom_smooth()+
  stat_smooth(method = "nls",
            formula = y ~ a/(1+b*x), # BH
            method.args = list(start = list(a = 500, b = .1)),
            se = FALSE, 
            color='orange')
  

### ### ### ### ### ### ### ### 
# fit gaussian Bev-Holt ----
### ### ### ### ### ### ### ### 

fit1 <- brm(
  bf(fitness_data ~ lambda / (1+ alpha * density_data) 
     , alpha ~ 1 
     , lambda ~ 1 # treatment + (1|block)
     , nl=TRUE),
  data = dat, 
  prior = c(prior(gamma(20, 1), lb=0, nlpar = "lambda"), 
            prior(normal(0, 1), nlpar = "alpha")
  )
  ,inits = "0" 
  ,cores=4 
  ,chains=4
  ,iter=5000 
  # ,control = list(adapt_delta = 0.98)
)

plot(fit1)
fit1
conditional_effects(fit1) # quick plot

# Retrieve and plot coefficients
get_variables(fit1)

param.table <- fit1 %>% 
  gather_draws(`b_.*`, regex = TRUE) %>% 
  median_qi(.width = c(.95))

param.table

# plot fit to data
dat %>%  
  modelr::data_grid(density_data = modelr::seq_range(density_data, n = 101)) %>%
  add_fitted_draws(fit1, n = 100) %>%
  ggplot(aes(x = density_data, y = .value)) +
  # geom_line(aes(y = .value, group = .draw), alpha = 0.1, size = 0.5, color = "blue") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .8), alpha = 1) +  # , show.legend=FALSE
  scale_fill_brewer(palette = "Blues")+
  geom_point(data = dat, aes(x = density_data, y = fitness_data))+
  geom_function(fun=bev1, color='red', size=1.5, args=c("lam"=param.table$.value[2], "alpha"=param.table$.value[1]), alpha=.5) + # this function plots fitted line using estimated parameter values from gaussian model
  theme_minimal()

# GOOD:  complete overlap between predictions using model and that using parameter estimates


# fit Poisson model ----
fit1.pois <- brm(
  bf(as.integer(fitness_data) ~ lambda / (1+ alpha * density_data) # Main brms model formulation; syntax is like lmer usually, but for nonlinear models need to write out equations using parameters
     , alpha ~ 1 # For nonlinear models, need formulas for parameters (can model them as functions of other things; not doing that here, just an intercept)
     , lambda ~ 1 # treatment + (1|block)
     , nl=TRUE),
  data = dat, # dataframe with data; needs to include the data referred to in model formula
  prior = c(prior(gamma(20, 1), lb=0, nlpar = "lambda"), # Nonlinear models require priors for the fixed effects
            prior(normal(0, 1), nlpar = "alpha")
            #, prior(student_t(3, 0, 1.5), class="sigma")
  )
  ,family=poisson("log")
  ,inits = "0" # don't have to set initial values, but can
  ,cores=4 
  ,chains=4
  ,iter=5000 
  # ,control = list(adapt_delta = 0.98)
)

fit1.pois


param.table <- fit1.pois %>% 
  gather_draws(`b_.*`, regex = TRUE) %>% 
  group_by(.variable) %>%
  summarize(
    median=median(.value),
    param.mean=mean(.value),
    param.stdev=sd(.value),
    pval= round(sapply(pnorm(0,param.mean,param.stdev), function(x) min(x, 1-x)), 4)
  )
param.table



# plot fit to data
dat %>%  
  modelr::data_grid(density_data = modelr::seq_range(density_data, n = 101)) %>%
  add_fitted_draws(fit1.pois, n = 100, scale="response") %>%
  ggplot(aes(x = density_data, y = .value)) +
  # geom_line(aes(y = .value, group = .draw), alpha = 0.1, size = 0.5, color = "blue") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .8), alpha = 0.25) +  # , show.legend=FALSE
  scale_fill_brewer(palette = "Blues")+
  geom_point(data = dat, aes(x = density_data, y = fitness_data))+
   geom_function(fun=bev1, color='red', size=1.5, args=c("lam"=param.table$median[2], "alpha"=param.table$median[1]))+
#  geom_function(fun=bev1, color='orange', size=1.5, args=c("lam"=param.table$median[2], "alpha"=exp(param.table$median[1])))+
#  geom_function(fun=bev1, color='yellow', size=1.5, args=c("lam"=exp(param.table$median[2]), "alpha"=.12))+
  theme_minimal()

# Here's the issue: model is good, predicts data well.   But can't just use the parameter estimates from the poisson model, at least not the alpha; the lambda looks ok with a exp() transformation.  
# How to transform the alpha then?



# other way to make predictions and plot ----
dat.new <- expand.grid(
  density_data = seq(0,30, length.out=50)
)

pred.gaussian <- 
  as.data.frame(fitted(fit1, newdata = dat.new, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new) 

pred.poisson <- 
  as.data.frame(fitted(fit1.pois, newdata = dat.new, allow_new_levels=TRUE, probs=c(.05,.5,.95)))  %>%
  cbind(dat.new) 

pred.all <- pred.gaussian %>% 
  mutate(poisson = pred.poisson$Estimate)
head(pred.all); dim(pred.all)
pred.all <- pred.all %>% rename(gaussian = Estimate) %>%
  pivot_longer(cols=c(gaussian, poisson), names_to="model", values_to = "fitness")
pred.all

a <- pred.all %>%  ggplot(aes(x = density_data, y = fitness)) + 
  geom_line(aes(linetype=model)) + 
  geom_jitter(data=dat, aes(x=density_data, y=fitness_data))+
  ylab("fitness")+
  theme_minimal()
a



# trying to use fitted to recover values?   Doesn't work
vals <- fitted(fit1.pois, nlpar="alpha", scale="linear")
head(vals)
vals <- fitted(fit1.pois, nlpar="alpha", scale="response")
head(vals)
vals <- fitted(fit1.pois, scale="linear")
head(vals)
vals <- fitted(fit1.pois, scale="response")
head(vals)



## make predicted vs observed graph -----

fe <- as.data.frame(fixef(fit1.pois))

# make comparison using predict function
predict.data <- predict(fit1.pois, dat)
dd <- cbind(dat, predict.data)
ggplot(dd, aes(x=fitness_data, y=Estimate)) +
  geom_point()+
  geom_abline(slope=1, intercept=1)


# make it using parameter estimates (Doesn't work)
predict.data <- bev1(dat$density_data, lam=fe$Estimate[2],alpha=fe$Estimate[1] )
dd <- cbind(dat, predict.data)
ggplot(dd, aes(x=fitness_data, y=predict.data)) +
  geom_point()+
  geom_abline(slope=1, intercept=1)

# need to exp() it
predict.data <- exp(predict.data)
dd <- cbind(dat, predict.data)
ggplot(dd, aes(x=fitness_data, y=predict.data)) +
  geom_point()+
  geom_abline(slope=1, intercept=1)
