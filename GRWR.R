
# ------------------------------------------------------------------------------------
# pull in pars, visualize

par<-c("b_lambdaA_Intercept", "b_alphaAA_Intercept", "b_alphaAP_Intercept", "b_alphaAS_Intercept")
par0<-c("b_lambdaA_Intercept", "b_lambdaA_warmtrtwarm", "b_alphaAA_Intercept", "b_alphaAP_Intercept", "b_alphaAS_Intercept")
par1<-c("b_lambdaA_Intercept", "b_lambdaA_warmtrtwarm", "b_alphaAA_Intercept","b_alphaAA_warmtrtwarm", "b_alphaAP_Intercept","b_alphaAP_warmtrtwarm", "b_alphaAS_Intercept", "b_alphaAS_warmtrtwarm")

annual0 <- posterior_summary(A, pars=par0, probs = c(0.025, 0.975), robust = FALSE)
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

annual.summary   <-posterior_summary(annual.simple.gaussian, probs = c(0.025, 0.975), robust = FALSE)
adult.summary    <-posterior_summary(adult.simple.gaussian, probs = c(0.025, 0.975), robust = FALSE)
seedling.summary <-posterior_summary(seedling.simple, probs = c(0.025, 0.975), robust = FALSE)

annual.simple.gaussian
adult.simple.gaussian
seedling.simple


lambdaAa<-  annual.summary[1,1] #previously 1.30
lambdaAw<- annual.summary[1,1]+annual.summary[2,1] #previously 1.54 
lambdaSa<- seedling.summary[1,1] # .019
lambdaSw<- seedling.summary[1,1]+seedling.summary[2,1] # .019
lambdaPa<- adult.summary[1,1] #1.46
lambdaPw<- adult.summary[1,1]+adult.summary[2,1] #1.73
  
alphaAAa<- annual.summary[3,1] #.12
alphaAAw<- annual.summary[3,1]+annual.summary[4,1] # 28
#alphaASa<- .08#annual.amb[] #not in this model anymore
#alphaASw<- .05#annual.warm[] #not in this model anymore
alphaAPa<- annual.summary[5,1] #.10
alphaAPw<- annual.summary[5,1]+annual.summary[6,1] #.10
  
alphaPAa<- adult.summary[3,1] #.0012
alphaPAw<-  adult.summary[3,1] + adult.summary[4,1] #.0022
alphaPPa<-  adult.summary[5,1]# .0 
alphaPPw<- adult.summary[5,1]+adult.summary[6,1] # .17
  
alphaSAa<- seedling.summary[1,1] #.032
alphaSAw<- seedling.summary[1,1]+seedling.summary[2,1] #.011
alphaSSa<- seedling.summary[5,1] # .0042
alphaSSw<- seedling.summary[5,1]+seedling.summary[6,1] #.0012
alphaSPa<- seedling.summary[3,1] #.097
alphaSPw<- seedling.summary[3,1]+seedling.summary[4,1] # .094
  

# ------------------------------------------------------------------------------------
## Set germination and survival fractions from the literature
sa <- .11 # ghersa 1984 # costa maia 2009 minus 4% permanently dormant
ga <- .89 #ghersa 1984
ga <- 100 #gundel 2007, gundel 2006
ga <- .97 #lin 2018

sp <- 1 # 100% adults survive observed in field - x in literature
sp <-.975 # or 2.5% fiegner 2007
ss <- .0
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
  return(Na)
}

# Determine equilibrium conditions for festuca adults
adult.equilibrium <- function (N0p, N0s, sp, alphaSS, alphaSP, lambdaS) { # n adults, n seedlings, self seedling alpha, adult on seedling alpha, seedling survival to adult
  Np <- sp*N0p + N0s*lambdaS/(1 + alphaSS*N0s + alphaSP*N0p) 
  return(Np)
}

# Determine equilibrium conditions for festuca seeds
seedling.equilibrium <- function (N0p, N0s, ss, gs, alphaPP, lambdaP) { # n adult n seedling, seedling survival, germination seedling, self and adult alphas, adult seed production
  Ns <- ss*(1-gs)*N0s + N0p*(lambdaP*5000)/(1 + alphaPP*N0p)
  return(Ns)
}


#competition growth:
# annuals
annual.invade <- function  (N0a, N0p, N0s, sa, ga, alphaAA, alphaAP, lambdaA){
  Na <- sa*(1-ga)*N0a + N0a*lambdaA*100/(1+alphaAP*N0p+alphaAA*N0a)
  return(Na)
}

# seedlings
seedling.resident <- function (N0p, N0a, N0s, ss, gs, alphaPP, alphaPA, lambdaP) {
  Ns <- ss*(1-gs)*N0s + N0p*(lambdaP*5000)/(1 + alphaPP*N0p+ alphaPA*N0a)  
  return(Ns)
}

# adults
adult.resident <- function (N0p, N0s, N0a, sp, alphaSS, alphaSP, alphaSA,lambdaS) {
  Np <- sp*N0p + N0s*lambdaS/(1 + alphaSS*N0s + alphaSP*N0p+ alphaSA*N0a) 
  return(Np)
}


#
# ------------------------------------------------------------------------------------
# run equilibriums

time <- 1:300
tmax=max(time)

# annuals
a <- rep(0, length(time))

N0a=1
eq.annuals.a <- tibble(time, a)
eq.annuals.a[1,2] = as.numeric(N0a)


for (t in 1:tmax) {
  eq.annuals.a[t+1,2] <- annual.equilibrium(eq.annuals.a[t,2], sa, ga, alphaAAa, lambdaAa)
}

eq.annuals.w<- tibble(time, a)
eq.annuals.w[1,2] <- N0a

for (t in 1:tmax) {
  eq.annuals.w[t+1,2] <- annual.equilibrium(eq.annuals.w[t,2], sa, ga, alphaAAw, lambdaAw)
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

for (t in 1:300) {
  eq.perennials.a[t+1, 2] <- adult.equilibrium(eq.perennials.a[t, 2], eq.perennials.a[t, 3], sp, alphaSSa, alphaSPa, lambdaSa)
  eq.perennials.a[t+1, 3] <- seedling.equilibrium(eq.perennials.a[t, 2], eq.perennials.a[t, 3], ss, gs, alphaPPa, lambdaPa)
}


#warm
eq.perennials.w<-tibble(time, p, s)
eq.perennials.w[1,3] <- N0s
eq.perennials.w[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.w[t+1, 2] <- adult.equilibrium(eq.perennials.w[t, 2], eq.perennials.w[t, 3], sp, alphaSSw, alphaSPa, lambdaSw)
  eq.perennials.w[t+1, 3] <- seedling.equilibrium(eq.perennials.w[t, 2], eq.perennials.w[t, 3], ss, gs, alphaPPw, lambdaPw)
}



# check outputs
annual.eq<-rbind(mutate(eq.annuals.a, temp="amb"), mutate(eq.annuals.w, temp="warm"))
perennial.eq<-rbind(mutate(eq.perennials.a, temp="amb"), mutate(eq.perennials.w, temp="warm"))%>%
  gather("type", "density", p, s)%>%
  mutate(type=ifelse(type=="p", "adult", "seedling"))

eqp1<-ggplot(annual.eq, aes(time, a, color=temp))+ geom_line()+ylab("annual seed density")+  scale_colour_manual(values = c("dodgerblue", "darkred"))
eqp2<-ggplot(perennial.eq, aes(time, (density), color=temp, linetype=type, shape=type))+ geom_line()+ylab("perennial density")+scale_y_continuous(trans='log10', labels = scales::comma)+  scale_colour_manual(values = c("dodgerblue", "darkred"))

ggarrange(eqp1, eqp2)

# ------------------------------------------------------------------------------------
# invasion

# annual into perennials 
#ambient
annual.invasion.a <- tibble(time, a, p, s)
annual.invasion.a[1,2] <- 1 #annuals
annual.invasion.a[1,3] <- eq.perennials.a[300,2] #adults at eq
annual.invasion.a[1,4] <- eq.perennials.a[300,3] #seedlings at eq


for (t in 1:tmax) {
  annual.invasion.a[t+1, 2] <- annual.invade(N0a=annual.invasion.a[t,2],
                                             N0p=annual.invasion.a[t,3],N0s=annual.invasion.a[t,4],
                                    sa, ga, alphaAAa, alphaAPa, lambdaAa)
  annual.invasion.a[t+1, 4] <- seedling.resident(N0a=annual.invasion.a[t,2],
                                      N0p=annual.invasion.a[t,3],N0s=annual.invasion.a[t,4],
                                      ss, gs, alphaPPa, alphaPAa, lambdaPa)
  annual.invasion.a[t+1, 3] <- adult.resident(N0a=annual.invasion.a[t,2],
                                      N0p=annual.invasion.a[t,3],
                                      N0s=annual.invasion.a[t,4], 
                                      sp, alphaSSa, alphaSPa, alphaSAa, lambdaSa)
  }

#warmed
annual.invasion.w <- tibble(time, a, p, s)
annual.invasion.w[1,2] <- 1
annual.invasion.w[1,3] <- eq.perennials.w[300,2] #adults at eq
annual.invasion.w[1,4] <- eq.perennials.w[300,3] #seedlings at eq

for (t in 1:tmax) {
  annual.invasion.w[t+1, 2] <- annual.invade(N0a=annual.invasion.w[t,2],
                                             N0p=annual.invasion.w[t,3],N0s=annual.invasion.w[t,4],
                                             sa, ga, alphaAAw, alphaAPw, lambdaAw)
  annual.invasion.w[t+1, 4] <- seedling.resident(N0a=annual.invasion.w[t,2],
                                                 N0p=annual.invasion.w[t,3],N0s=annual.invasion.w[t,4],
                                                 ss, gs, alphaPPw, alphaPAw, lambdaPw)
  annual.invasion.w[t+1, 3] <- adult.resident(N0a=annual.invasion.w[t,2],
                                              N0p=annual.invasion.w[t,3],
                                              N0s=annual.invasion.w[t,4], 
                                              sp, alphaSSw, alphaSPw, alphaSAw, lambdaSw)
}


# perennials into annuals

#ambient
perennial.invasion.a <- tibble(time, a, p, s)
perennial.invasion.a[1,2] <-  eq.annuals.a[tmax,2] #annuals at eq
perennial.invasion.a[1,3] <- .5
perennial.invasion.a[1,4] <- .5


for (t in 1:tmax) {
  perennial.invasion.a[t+1, 2] <- annual.invade(N0a=perennial.invasion.a[t,2],
                                                N0p=perennial.invasion.a[t,3],N0s=perennial.invasion.a[t,4],
                                                sa, ga, alphaAAa, alphaAPa, lambdaAa)
  perennial.invasion.a[t+1, 4] <- seedling.resident(N0a=perennial.invasion.a[t,2],
                                                    N0p=perennial.invasion.a[t,3],N0s=perennial.invasion.a[t,4],
                                                    ss, gs, alphaPPa, alphaPAa, lambdaPa)
  perennial.invasion.a[t+1, 3] <- adult.resident(N0a=perennial.invasion.a[t,2],
                                                N0p=perennial.invasion.a[t,3],
                                                N0s=perennial.invasion.a[t,4], 
                                                sp, alphaSSa, alphaSPa, alphaSAa, lambdaSa)
}

#warmed
perennial.invasion.w <- tibble(time, a, p, s)
perennial.invasion.w[1,2] <-  eq.annuals.w[tmax,2] #annuals at eq
perennial.invasion.w[1,3] <- .5
perennial.invasion.w[1,4] <- .5


for (t in 1:tmax) {
  perennial.invasion.w[t+1, 2] <- annual.invade(N0a=perennial.invasion.w[t,2],
                                                N0p=perennial.invasion.w[t,3],N0s=perennial.invasion.w[t,4],
                                                sa, ga, alphaAAw, alphaAPw, lambdaAw)
  perennial.invasion.w[t+1, 4] <- seedling.resident(N0a=perennial.invasion.w[t,2],
                                                    N0p=perennial.invasion.w[t,3],N0s=perennial.invasion.w[t,4],
                                                    ss, gs, alphaPPw, alphaPAw, lambdaPw)
  perennial.invasion.w[t+1, 3] <- adult.resident(N0a=perennial.invasion.w[t,2],
                                                 N0p=perennial.invasion.w[t,3],
                                                 N0s=perennial.invasion.w[t,4], 
                                                 sp, alphaSSw, alphaSPw, alphaSAw, lambdaSw)
}
  

#visualize

# check outputs
annual.inv<-rbind(mutate(annual.invasion.a, temp="amb"), mutate(annual.invasion.w, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))


perennial.inv<-rbind(mutate(perennial.invasion.a, temp="amb"), mutate(perennial.invasion.w, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))

inv1<-ggplot(subset(annual.inv, time<10), aes(time, density, color=temp, linetype=type, shape=type))+   scale_colour_manual(values = c("dodgerblue", "darkred"))+
geom_line()+ylab("density")+scale_y_continuous(trans='log10')
inv2<-ggplot(subset(perennial.inv, time<50), aes(time, density, color=temp, linetype=type, shape=type))+   scale_colour_manual(values = c("dodgerblue", "darkred"))+
geom_line()+ylab("density")+scale_y_continuous(trans='log10')

ggarrange(inv1, inv2)


  # pars contains the following parameters (subscript 1 is the annual, 2 and 3 are the perennial seedling and adult, respectively):
  # g1, g2 = germination fraction
  # lambda1, lambda3 = seed production in the absence of competition
  # alpha11, alpha13, etc. = competition coefficients
  # phi = factor scaling competitive effects of annuals to competitive effects of perennial seedlings
  # s2, s3 = over-summer survival


annualGRWRa <- log(annual.invasion.a[2,2])%>%
  mutate(invader="Invader: Lolium", trt="ambient")
annualGRWRw <- log(annual.invasion.w[2,2])%>%
  mutate(invader="Invader: Lolium", trt="warmed")
perennialGRWRa = log(.5*(sp + (4* lambdaPa*5000/(1 + alphaPAa*perennial.invasion.a[1,2]) * lambdaSa/(1 + alphaSAa*perennial.invasion.a[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="ambient")
perennialGRWRw = log(.5*(sp + (4* lambdaPw*5000/(1 + alphaPAw*perennial.invasion.w[1,2]) * lambdaSw/(1 + alphaSAw*perennial.invasion.w[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="warmed")


allGRWR<-rbind(annualGRWRa, annualGRWRw, perennialGRWRa, perennialGRWRw)%>%
  mutate(GRWR=a)%>%
  select(-a)

ggplot(allGRWR, aes(x=trt, y=GRWR)) +geom_bar(stat="identity") +facet_wrap(~invader, scales="free") +xlab("") +ylab("Invader Growth Rate (r) When Rare")

#method 2, need to create a 2x2 transition matrix with (seed survival term, seed production term )
#                                                      (seeds maturing term, adults surviving term)
library(popbio)

eigen.analysis(perennial.invasion.a[1:2,3:4])
eigen.analysis(perennial.invasion.w[1:2,3:4])

