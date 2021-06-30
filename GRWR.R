# ------------------------------------------------------------------------------------
# Functions for use in coexistence calcualtions

# Determine equilibrium conditions for lolium seeds - to run for only a single timestep
pop.equilibrium <- function (N0a, sa, ga, alphaAA, lambdaA) { #number of annuals, seed survival annuals, germination annuals, self alpha, annual lambda
  Na <- sa*(1-ga)*N0a + N0a*lambdaA*30/(1+alphaAA*N0a) #predicted lambda is scaled, match scaling here
  return(N)
}

# Determine equilibrium conditions for festuca adults
pop.equilibrium <- function (N0p, NOs, sp, alphaSS, alphaSP, lambdaS) { # n adults, n seedlings, self seedling alpha, adult on seedling alpha, seedling survival to adult
  Np <- sp*N0p + N0s*lambdaS/(1 + alphaSS*seeded_s + alphaSP*density_p) 
  return(N)
}

# Determine equilibrium conditions for festuca seeds
pop.equilibrium <- function (N0p, NOs, ss, gs, alphaSS, alphaSP, lambdaP) { # n adult n seedling, seedling survival, germination seedling, self and adult alphas, adult seed production
  Ns <- ss*(1-gs)*N0s + N0p*(lambdaP*5000)/(1 + alphaPP*N0p)
  return(N)
}

# invader population growth rate one time step forward
pop.invade <- function (N0, resident, s, g, a_inter, lambda) {
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_inter*resident)
  return(N)
}

# resident population growth rate one time step forward
pop.resident <- function (N0, resident, s, g, a_intra, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*resident + resident*(lambda*g)/(1+a_intra*resident+a_inter*N0)
  return(N)
}

# ------------------------------------------------------------------------------------
# run models
# determine resident equilibrium abundances and low density growth rates

annual <- readA
seedling <- readPS
adult <- read PA

## Set germination and survival fractions from the literature
sa <- 
ga <- 
sp <- 1 # 100% adults survive observed in field - x in literature
ss <- 
gs <- 

# use the timeseries of environmental conditions for environmental variability
# for avena
N0 <- 550
time <- length(rainsummary$raintype)
N_avena <- rep(NA, time)
N_avena[1] <- N0

for (t in 1:time) {
  params <- subset(avena, treatment==rainsummary$raintype[t])
  N_avena[t+1] <- pop.equilibrium(N0=N_avena[t], s=as, g=ag, a_intra=params$aiA, lambda=params$lambda)
}

# check output
plot(seq(1:(time+1)), N_avena, type="l")


# invade avena first
avena_invade <- rep (NA, 72)
erodium_resident <- rep (NA, 72)
temp <- 1
for (t in 50:time) {
  params <- subset(avena, treatment==rainsummary$raintype[t])
  params_resident <- subset(erodium, treatment==rainsummary$raintype[t])
  avena_invade[temp] <- pop.invade(N0=1, resident=N_erodium[t], s=as, g=ag, a_inter=params$aiE, lambda=params$lambda)
  
  # sanity check that the resident isn't affected
  erodium_new <- pop.resident(N0=1, resident=N_erodium[t], s=es, g=eg, 
                              a_intra=params_resident$aiE, a_inter=params_resident$aiA, 
                              lambda=params_resident$lambda)
  erodium_resident[temp] <- erodium_new/N_erodium[t]
  
  temp  <- temp + 1 
}


avena_invader <- log(avena_invade)
erodium_invader <- log(erodium_invade)

avena_r <- log(avena_resident)
erodium_r <- log(erodium_resident)