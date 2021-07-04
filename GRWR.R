
# ------------------------------------------------------------------------------------
# pull in pars, visualize

par<-c("b_lambdaA_Intercept", "b_alphaAA_Intercept", "b_alphaAP_Intercept", "b_alphaAS_Intercept")
par0<-c("b_lambdaA_Intercept", "b_lambdaA_warmtrtwarm", "b_alphaAA_Intercept", "b_alphaAP_Intercept", "b_alphaAS_Intercept")
par1<-c("b_lambdaA_Intercept", "b_lambdaA_warmtrtwarm", "b_alphaAA_Intercept","b_alphaAA_warmtrtwarm", "b_alphaAP_Intercept","b_alphaAP_warmtrtwarm", "b_alphaAS_Intercept", "b_alphaAS_warmtrtwarm")

annual0 <- posterior_summary(annual.model, pars=par0, probs = c(0.025, 0.975), robust = FALSE)
annual<-as.tibble(annual0)%>%
  mutate(parameters=rownames(annual0)) %>%
  separate(parameters, sep="_", into=c("b", "par", "trt"))

ggplot(annual, aes(x=trt, y=Estimate*30)) +geom_point() +
  geom_errorbar(aes(ymin=Estimate*30-Est.Error*30, ymax=Estimate*30+Est.Error*30))+
  facet_wrap(~par, scales="free")

annualv0<- posterior_summary(annual.model.var, pars=par1, probs = c(0.025, 0.975), robust = FALSE)
annualv<-as.tibble(annualv0)%>%
  mutate(parameters=rownames(annualv0)) %>%
  separate(parameters, sep="_", into=c("b", "par", "trt"))

ggplot(annualv, aes(x=trt, y=Estimate*30)) +geom_point() +
  geom_errorbar(aes(ymin=Estimate*30-Est.Error*30, ymax=Estimate*30+Est.Error*30))+
  facet_wrap(~par, scales="free")

annual.warm0 <- (posterior_summary(annual.warm.model, pars=par, probs = c(0.025, 0.975), robust = FALSE))
annual.warm<-as.tibble(annual.warm0)%>%
  mutate(parameters=rownames(annual.warm0)) %>%
  mutate(parameters=paste(parameters, "_warm"))
annual.amb0  <- posterior_summary( annual.amb.model, pars=par, probs = c(0.025, 0.975), robust = FALSE)
annual.amb<-as.tibble(annual.amb0)%>%
  mutate(parameters=rownames(annual.amb0)) %>%
  mutate(parameters=paste(parameters, "_amb"))
annual.bind<-rbind(annual.warm, annual.amb)%>%
  separate(parameters, sep="_", into=c("b", "par", "int", "trt"))

ggplot(annual.bind, aes(x=trt, y=Estimate*30)) +geom_point() +
  geom_errorbar(aes(ymin=Estimate*30-Est.Error*30, ymax=Estimate*30+Est.Error*30))+
  facet_wrap(~par, scales="free")

ggplot(annual.bind, aes(x=trt, y=Estimate*30)) +geom_point() +
  geom_errorbar(aes(ymin=Q2.5*30, ymax=Q97.5*30))+
  facet_wrap(~par, scales="free")



seedling.warm <- posterior_summary(annual.warm.model, probs = c(0.025, 0.975), robust = FALSE)
seedling.amb  <- posterior_summary( annual.amb.model, probs = c(0.025, 0.975), robust = FALSE)

adult.warm <- posterior_summary(annual.warm.model, probs = c(0.025, 0.975), robust = FALSE)
adult.amb  <- posterior_summary( annual.amb.model, probs = c(0.025, 0.975), robust = FALSE)

lambdaAa<- annual.amb[]
lambdaAw<- annual.warm[]
lambdaSa<- seedling.amb[]
lambdaSw<- seedling.warm[]
lambdaPa<- adult.amb []
lambdaPw<- adult.warm[]
  
alphaAAa<- annual.amb[]
alphaAAw<- annual.warm[]
alphaASa<- annual.amb[]
alphaASw<- annual.warm[]
alphaAPa<- annual.amb[]
alphaAPw<- annual.warm[]
  
alphaPAa<- adult.amb []
alphaPAw<- adult.warm []
alphaPPa<- adult.amb []
alphaPPw<- adult.warm []
  
alphaSAa<- seedling.amb[]
alphaSAw<- seedling.warm[]
alphaSSa<- seedling.amb[]
alphaSSw<- seedling.warm[]
alphaSPa<- seedling.amb[]
alphaSPw<- seedling.warm[]
  

# ------------------------------------------------------------------------------------
## Set germination and survival fractions from the literature
sa <- .11 # ghersa 1984 # costa maia 2009 minus 4% permanently dormant
ga <- .89 #ghersa 1984
ga <- 100 #gundel 2007, gundel 2006
ga <- .97 #lin 2018

sp <- 1 # 100% adults survive observed in field - x in literature
sp <-.99 # or 2.5% fiegner 2007
ss <- 
gs <- .63 # fiegner 2007
gs <- .90 # schmidt 1998
gs <- .31 # maret 2005
gs <- .50 # mackin 2021

# ------------------------------------------------------------------------------------
# Functions for use in coexistence calcualtions

#equilibriums 

# Determine equilibrium conditions for lolium seeds - to run for only a single timestep
annual.equilibrium <- function (N0a, sa, ga, alphaAA, lambdaA) { #number of annuals, seed survival annuals, germination annuals, self alpha, annual lambda
  Na <- sa*(1-ga)*N0a + N0a*lambdaA*100/(1+alphaAA*N0a) #predicted lambda is scaled, match scaling here 100, 30??
  return(N)
}

# Determine equilibrium conditions for festuca adults
adult.equilibrium <- function (N0p, N0s, sp, alphaSS, alphaSP, lambdaS) { # n adults, n seedlings, self seedling alpha, adult on seedling alpha, seedling survival to adult
  Np <- sp*N0p + N0s*lambdaS/(1 + alphaSS*N0s + alphaSP*N0p) 
  return(N)
}

# Determine equilibrium conditions for festuca seeds
seedling.equilibrium <- function (N0p, N0s, ss, gs, alphaPP, lambdaP) { # n adult n seedling, seedling survival, germination seedling, self and adult alphas, adult seed production
  Ns <- ss*(1-gs)*N0s + N0p*(lambdaP*5000)/(1 + alphaPP*N0p)
  return(N)
}



#annual invading perennials

# annual invader one time step forward
annual.invade <- function  (N0a, N0p, N0s, sa, ga, alphaAA, alphaAS, alphaAP, lambdaA){
  Na <- sa*(1-ga)*N0a + N0a*lambdaA*100/(1+alphaAS*N0s+alphaAP*N0p)
  return(N)
}

# seedling one step forward
seedling.resident <- function (N0p, N0a, ss, gs, alphaPP, alphaPA, lambdaP) {
  Ns <- ss*(1-gs)*N0s + N0p*(lambdaP*5000)/(1 + alphaPP*N0p+ alphaPA*N0a)  
  return(N)
}

# adult one step forward
adult.resident <- function (N0p, N0s, N0a, sp, alphaSS, alphaSP, alphaSA,lambdaS) {
  Np <- sp*N0p + N0s*lambdaS/(1 + alphaSS*N0s + alphaSP*N0p+ alphaSA*N0a) 
  return(N)
}


#
# ------------------------------------------------------------------------------------
# run equilibriums

time <- 1:50
tmax=max(time)

# annuals
a <- rep(0, length(time))

N0a=1
eq.annuals.a <- tibble(time, a)
eq.annuals.a[1,2] = as.numeric(N0a)


for (t in 1:tmax) {
  eq.annuals.a[t+1] <- pop.equilibrium(eq.annuals[t], sa, ga, alphaAAa, lambdaAa)
}

eq.annuals.w<- tibble(time, a)
eq.annuals.w[1] <- N0a

for (t in 1:tmax) {
  eq.annuals.w[t+1] <- pop.equilibrium(eq.annuals[t], sa, ga, alphaAAw, lambdaAw)
}


#perennials

N0s=1
N0p=1
s <- rep(0, length(time))
p <- rep(0, length(time))

#amb
eq.perennials.a<-tibble(time, p, s)
eq.perennials.a[1,3] <- N0s
eq.perennials.a[1,2] <- N0p

for (t in 1:tmax) {
  eq.perennials.a[t+1, 2] <- pop.equilibrium(eq.perennials.a[t, 2], eq.perennials.a[t, 3], sp, alphaSSa, lambdaSa)
  eq.perennials.a[t+1, 3] <- pop.equilibrium(eq.perennials.a[t, 2], eq.perennials.a[t, 3], ss, gs, alphaSSa, alphaSPa, lambdaPa)
}

#warm
eq.perennials.w<-tibble(time, p, s)
eq.perennials.w[1,3] <- N0s
eq.perennials.w[1,2] <- N0p

for (t in 1:tmax) {
  eq.perennials.w[t+1, 2] <- pop.equilibrium(eq.perennials.w[t, 2], eq.perennials.w[t, 3], sp, alphaSSw, lambdaSw)
  eq.perennials.w[t+1, 3] <- pop.equilibrium(eq.perennials.w[t, 2], eq.perennials.w[t, 3], ss, gs, alphaSSw, alphaSPw, lambdaPw)
}

# check outputs
plot lines for each 

# ------------------------------------------------------------------------------------
# invasion

# annual into perennials 
#ambient
annual.invasion.a <- tibble(time, a, s, p)
annual.invasion.a[1,2] <- 1 #annuals
annual.invasion.a[1,3] <- eq.perennials.a[tmax,2] #adults at eq
annual.invasion.a[1,4] <- eq.perennials.a[tmax,3] #seedlings at eq


for (t in tmax) {
  avena_invasion.a[t+1, 2] <- annual.invade(N0a=annual.invasion.a[t,2],
                                    N0p=annual.invasion.a[t,3],
                                    N0s=annual.invasion.a[t,4], 
                                    sa, ga, alphaAAa, alphaASa, alphaAPa, lambdaAa)
  avena_invasion.a[t+1, 3] <- seedling.resident(N0a=annual.invasion.a[t,2],
                                      N0p=annual.invasion.a[t,3],
                                      ss, gs, alphaPPa, alphaPAa, lambdaPa)
  avena_invasion.a[t+1, 4] <- adult.resident(N0a=annual.invasion.a[t,2],
                                      N0p=annual.invasion.a[t,3],
                                      N0s=annual.invasion.a[t,4], 
                                      sp, alphaSSa, alphaSPa, alphaSAa, lambdaSa)
}


#warmed
annual.invasion.w <- tibble(time, a, s, p)
annual.invasion.w[1,2] <- 1
annual.invasion.w[1,3] <- eq.perennials.w[tmax,2] #adults at eq
annual.invasion.w[1,4] <- eq.perennials.w[tmax,3] #seedlings at eq

for (t in tmax) {
  avena_invasion.w[t+1, 2] <- annual.invade(N0a=annual.invasion.a[t,2],
                                            N0p=annual.invasion.a[t,3],
                                            N0s=annual.invasion.a[t,4], 
                                            sa, ga, alphaAAw, alphaASw, alphaAPw, lambdaAw)
  avena_invasion.w[t+1, 3] <- seedling.resident(N0a=annual.invasion.a[t,2],
                                                N0p=annual.invasion.a[t,3],
                                                ss, gs, alphaPPw, alphaPAw, lambdaPw)
  avena_invasion.w[t+1, 4] <- adult.resident(N0a=annual.invasion.a[t,2],
                                             N0p=annual.invasion.a[t,3],
                                             N0s=annual.invasion.a[t,4], 
                                             sp, alphaSSw, alphaSPw, alphaSAw, lambdaSw)
}



# perennials into annuals

#ambient
perennial.invasion.a <- tibble(time, a, s, p)
perennial.invasion.a[1,2] <-  eq.annuals.a[tmax,2] #annuals at eq
perennial.invasion.a[1,3] <- 1
perennial.invasion.a[1,4] <- .1


for (t in tmax) {
  perennial.invasion.a[t+1, 2] <- annual.invade(N0a=annual.invasion.a[t,2],
                                            N0p=annual.invasion.a[t,3],
                                            N0s=annual.invasion.a[t,4], 
                                            sa, ga, alphaAAa, alphaASa, alphaAPa, lambdaAa)
  perennial.invasion.a[t+1, 3] <- seedling.resident(N0a=annual.invasion.a[t,2],
                                                N0p=annual.invasion.a[t,3],
                                                ss, gs, alphaPPa, alphaPAa, lambdaPa)
  perennial.invasion.a[t+1, 4] <- adult.resident(N0a=annual.invasion.a[t,2],
                                             N0p=annual.invasion.a[t,3],
                                             N0s=annual.invasion.a[t,4], 
                                             sp, alphaSSa, alphaSPa, alphaSAa, lambdaSa)
}


#warmed
perennial.invasion.w <- tibble(time, a, s, p)
perennial.invasion.w[1,2] <-  eq.annuals.w[tmax,2] #annuals at eq
perennial.invasion.w[1,3] <- 1
perennial.invasion.w[1,4] <- .1


for (t in tmax) {
  perennial.invasion.w[t+1, 2] <- annual.invade(N0a=annual.invasion.w[t,2],
                                                N0p=annual.invasion.w[t,3],
                                                N0s=annual.invasion.w[t,4], 
                                                sa, ga, alphaAAa, alphaASa, alphaAPa, lambdaAa)
  perennial.invasion.w[t+1, 3] <- seedling.resident(N0a=annual.invasion.w[t,2],
                                                    N0p=annual.invasion.w[t,3],
                                                    ss, gs, alphaPPa, alphaPAa, lambdaPa)
  perennial.invasion.w[t+1, 4] <- adult.resident(N0a=annual.invasion.w[t,2],
                                                 N0p=annual.invasion.w[t,3],
                                                 N0s=annual.invasion.w[t,4], 
                                                 sp, alphaSSa, alphaSPa, alphaSAa, lambdaSa)
  # pars contains the following parameters (subscript 1 is the annual, 2 and 3 are the perennial seedling and adult, respectively):
  # g1, g2 = germination fraction
  # lambda1, lambda3 = seed production in the absence of competition
  # alpha11, alpha13, etc. = competition coefficients
  # phi = factor scaling competitive effects of annuals to competitive effects of perennial seedlings
  # s2, s3 = over-summer survival

  grwr.a = lambdaA*500/(1 + alphaAS*pop.s + alphaAP*pop.p))
  grwr.p = .5*(adultsumsur + (4* lambdaP/(1 + alphaPA*pop.a) * seedsur/(1 + alphaSA*pop.a) + adult.sumsur^2)^0.5))

grwr.a = g1*lambda1/(1 + g2*alpha11*phi*res.p[1] + alpha13*res.p[2])
grwr.p = 0.5*(s3 + (4* lambda3/(1 + alpha31*res.a*g1) * g2*s2/(1 + alpha21*res.a*g1) + s3^2)^0.5)

#last log whatever to go from lambda to little r ()
avena_invader <- log(avena_invade)
erodium_invader <- log(erodium_invade)

avena_r <- log(avena_resident)
erodium_r <- log(erodium_resident)