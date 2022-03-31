//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

//annual.simple.gaussian <- brm(bf(percap ~  lambdaA*100 / (1 + alphaAA*seeded_am2 + alphaAP*starting_pm2), # + alphaAS*seeded_s #switched seeded_a to seeded_am2 to get density, not sure if to switch pm2 to starting_pm2
//                                 lambdaA + alphaAA + alphaAP ~ warmtrt + (1|time),
//                                 # + time,
//                                 nl=TRUE), 
//                              data = annuals,
//                              family = gaussian,prior = c(prior(normal(0, 1), nlpar = "lambdaA"),
//                                        #  set_prior("uniform(0, 2)", class="Intercept",  nlpar = "lambdaA"), #an attempt to just set a lower bound for intercept to be positive
//                                        prior(normal(0, .1), nlpar = "alphaAA"),
//                                        prior(normal(0, .1), nlpar = "alphaAP")),

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
}

