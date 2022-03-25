//Annual population growth model

// The data is all of length 'N'
// Data inputs include 'percap' per-capita seed production of annuals at the plot level
// seeded_am2, the number of annual seeds added per plot per meter squared
// starting_pm2, the number of adult perennials initially planted in each plot per meter squared
// warwmtrt, whether a plot was warmed or ambient
// time, whether data is from 2020 or 2021
data {
  int<lower = 1> N;
  int<lower=0> percap[N]; // turn percap into an integer by rounding
  vector[N] seeded_am2;
  vector[N] starting_pm2;
  int warmtrt[N]; //change to 0 vs. 1
  //int time[N];
  // int Plot[N];int<lower = 1>P;

}

// Our model accepts three parameters
// lambdaA,intrinsic growth rate (seeds-in, seeds-out) of annuals
// alphaAA, conspecific competition annuals
// alphaAP, adult perennial competitive effect on annuals
// amb refers to the ambient treatment, slope refers to marinal effect of warming
// amb + slope = warmed parameter
parameters {
  real lambdaA_amb;
  real lambdaA_slope;
  real alphaAA_amb;
  real alphaAA_slope;
  real alphaAP_amb;
  real alphaAP_slope;
//  real<lower=0> sigma; sigma has to do with random effects
// if we want to do separately for each treatment/year do like above with <lower=0>
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // intermediate objects
  vector[N] percap_hat; 
  vector[N] lambdaA_e; //for each data point, what is data for relevant environment
  vector[N] alphaAA_e;
  vector[N] alphaAP_e;

  // set priors
  lambdaA_amb ~ normal(4, .5);
  lambdaA_slope ~ normal(0, 1);
  alphaAA_amb ~ normal(-2, 1);
  alphaAA_slope ~ normal(0, 1);
  alphaAP_amb ~ normal(-2, 1);
  alphaAP_slope ~ normal(0, 1);  
  
  // define the model
  for(i in 1:N){
    lambdaA_e[i] = exp(lambdaA_amb + lambdaA_slope * warmtrt[i]);
    alphaAA_e[i] = exp(alphaAA_amb + alphaAA_slope * warmtrt[i]);
    alphaAP_e[i] = exp(alphaAP_amb + alphaAP_slope * warmtrt[i]);
    percap_hat[i] = lambdaA_e[i] / (1 + alphaAA_e[i] * seeded_am2[i] + alphaAP_e[i] * starting_pm2[i]);
  }
  percap ~ poisson(percap_hat);
}


