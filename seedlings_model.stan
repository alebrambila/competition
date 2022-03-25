
// The data is all of length 'N'
// Data inputs include 'fecundity' per-capita seed production of adults at the individual level
// seeded_am2, the number of annual seeds added per plot per meter squared
// starting_pm2, the number of adult perennials initially planted in each plot per meter squared
// warwmtrt, whether a plot was warmed or ambient
// time, whether data is from 2020 or 2021
data {
  int<lower = 1> N;
  int<lower=0> survivors[N]; 
  vector[N] seeded_am2;
    vector[N] seeded_sm2;
  vector[N] starting_pm2;
  int warmtrt[N]; //change to 0 vs. 1
  //int time[N];
  // int Plot[N];int<lower = 1>P;

}

// Our model accepts three parameters
// lambdaS,survival (seeds-in, adults out) of seedlings
// alphaSS, conspecific competition seedlings
// alphaSA, annual competitive effects on seedlings, alphaSP, aadult competitive effect on seedlings
// amb refers to the ambient treatment, slope refers to marinal effect of warming
// amb + slope = warmed parameter
parameters {
  real survivalS_amb;
  real survivalS_slope;
  real alphaSA_amb;
  real alphaSA_slope;
  real alphaSP_amb;
  real alphaSP_slope;
  real alphaSS_amb;
  real alphaSS_slope;
//  real<lower=0> sigma; sigma has to do with random effects
// if we want to do separately for each treatment/year do like above with <lower=0>
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // intermediate objects
  vector[N] survivors_hat; 
  vector[N] survivalS_e; //for each data point, what is data for relevant environment
  vector[N] alphaSA_e;
  vector[N] alphaSP_e;
  vector[N] alphaSS_e;


  // set priors
  survivalS_amb ~ normal(0, 1);
  survivalS_slope ~ normal(0, 1);
  alphaSA_amb ~ normal(0, 1);
  alphaSA_slope ~ normal(0, 1);
  alphaSP_amb ~ normal(0, 1);
  alphaSP_slope ~ normal(0, 1);  
    alphaSS_amb ~ normal(0, 1);
  alphaSS_slope ~ normal(0, 1);  
  
  // define the model
  for(i in 1:N){
    survivalS_e[i] = exp(survivalS_amb + survivalS_slope * warmtrt[i]);
    alphaSA_e[i] = exp(alphaSA_amb + alphaSA_slope * warmtrt[i]);
    alphaSP_e[i] = exp(alphaSP_amb + alphaSP_slope * warmtrt[i]);
      alphaSS_e[i] = exp(alphaSP_amb + alphaSP_slope * warmtrt[i]);
    survivors_hat[i] = survivalS_e[i] * seeded_sm2 / (1 + alphaSA_e[i] * seeded_am2[i] + alphaSP_e[i] * starting_pm2[i] + alphaSS_e[i] * seeded_sm2[i]);
  }
  survival ~ poisson(survival_hat);
}


