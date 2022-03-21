

// The data is all of length 'N'
// Data inputs include 'fecundity' per-capita seed production of adults at the individual level
// seeded_am2, the number of annual seeds added per plot per meter squared
// starting_pm2, the number of adult perennials initially planted in each plot per meter squared
// warwmtrt, whether a plot was warmed or ambient
// time, whether data is from 2020 or 2021
data {
  int<lower = 1> N;
  int<lower=0> fecundity[N]; 
  vector[N] seeded_am2;
  vector[N] starting_pm2;
  int warmtrt[N]; //change to 0 vs. 1
  //int time[N];
  // int Plot[N];int<lower = 1>P;

}

// Our model accepts three parameters
// lambdaP,fecundity (plugs-in, seeds-out) of adults
// alphaPP, conspecific competition adults
// alphaPAA, annual competitive effects on adults
// amb refers to the ambient treatment, slope refers to marinal effect of warming
// amb + slope = warmed parameter
parameters {
  real lambdaP_amb;
  real lambdaP_slope;
  real alphaPA_amb;
  real alphaPA_slope;
  real alphaPP_amb;
  real alphaPP_slope;
//  real<lower=0> sigma; sigma has to do with random effects
// if we want to do separately for each treatment/year do like above with <lower=0>
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // intermediate objects
  vector[N] fecundity_hat; 
  vector[N] lambdaP_e; //for each data point, what is data for relevant environment
  vector[N] alphaPA_e;
  vector[N] alphaPP_e;

  // set priors
  lambdaP_amb ~ normal(0, 1);
  lambdaP_slope ~ normal(0, 1);
  alphaPA_amb ~ normal(0, 1);
  alphaPA_slope ~ normal(0, 1);
  alphaPP_amb ~ normal(0, 1);
  alphaPP_slope ~ normal(0, 1);  
  
  // define the model
  for(i in 1:N){
    lambdaP_e[i] = exp(lambdaP_amb + lambdaP_slope * warmtrt[i]);
    alphaPA_e[i] = exp(alphaPA_amb + alphaPA_slope * warmtrt[i]);
    alphaPP_e[i] = exp(alphaPP_amb + alphaPP_slope * warmtrt[i]);
    fecundity_hat[i] = lambdaP_e[i] / (1 + alphaPA_e[i] * seeded_am2[i] + alphaPP_e[i] * starting_pm2[i]);
  }
  fecundity ~ poisson(fecundity_hat);
}


