
# ------------------------------------------------------------------------------------
# mean parameters ----

#fitsum.annuals
#fitsum.seedling
#fitsum.adult

#2020
lambdaAam20<-  annuals_estimates2020$exp[5]
lambdaAwm20<- annuals_estimates2020$exp[6]
survivalSam20<- seedlings_estimates2020$exp[7]
survivalSwm20<- seedlings_estimates2020$exp[8]
lambdaPam20<- adults_estimates2020$exp[5]
lambdaPwm20<- adults_estimates2020$exp[6] 
  
alphaAAam20<- annuals_estimates2020$exp[1]
alphaAAwm20<- annuals_estimates2020$exp[2]
alphaAPam20<- annuals_estimates2020$exp[3]
alphaAPwm20<- annuals_estimates2020$exp[4] 
  
alphaPAam20<- adults_estimates2020$exp[1]
alphaPAwm20<-  adults_estimates2020$exp[2]
alphaPPam20<-  adults_estimates2020$exp[3]
alphaPPwm20<- adults_estimates2020$exp[4]
  
alphaSAam20<- seedlings_estimates2020$exp[1]
alphaSAwm20<- seedlings_estimates2020$exp[2]
alphaSSam20<- seedlings_estimates2020$exp[5]
alphaSSwm20<- seedlings_estimates2020$exp[6]
alphaSPam20<- seedlings_estimates2020$exp[3]
alphaSPwm20<- seedlings_estimates2020$exp[4]

#2021
lambdaAam21<-  annuals_estimates2021$exp[5]
lambdaAwm21<- annuals_estimates2021$exp[6]
survivalSam21<- seedlings_estimates2021$exp[7]
survivalSwm21<- seedlings_estimates2021$exp[8]
lambdaPam21<- adults_estimates2021$exp[5]
lambdaPwm21<- adults_estimates2021$exp[6] 

alphaAAam21<- annuals_estimates2021$exp[1]
alphaAAwm21<- annuals_estimates2021$exp[2]
alphaAPam21<- annuals_estimates2021$exp[3]
alphaAPwm21<- annuals_estimates2021$exp[4] 

alphaPAam21<- adults_estimates2021$exp[1]
alphaPAwm21<-  adults_estimates2021$exp[2]
alphaPPam21<-  adults_estimates2021$exp[3]
alphaPPwm21<- adults_estimates2021$exp[4]

alphaSAam21<- seedlings_estimates2021$exp[1]
alphaSAwm21<- seedlings_estimates2021$exp[2]
alphaSSam21<- seedlings_estimates2021$exp[5]
alphaSSwm21<- seedlings_estimates2021$exp[6]
alphaSPam21<- seedlings_estimates2021$exp[3]
alphaSPwm21<- seedlings_estimates2021$exp[4]

#all chains (posterior predictive plots) ----
#annuals_fit2020, adults_fit2020, seedlings_fit2020

#2020
lambdaAam20<-  head(rstan::extract(annuals_fit2020)$lambdaA_amb)
dimnames(as.array(annuals_fit2020)$lambdaA_amb)
as.array(annuals_fit2020)[[lambdaA_amb]]

lambdaAwm20<- annuals_estimates2020$exp[6]
survivalSam20<- seedlings_estimates2020$exp[7]
survivalSwm20<- seedlings_estimates2020$exp[8]
lambdaPam20<- adults_estimates2020$exp[5]
lambdaPwm20<- adults_estimates2020$exp[6] 

allfits.annuals1<-allfits.annuals%>%
  filter(.iteration>1000)%>%
  mutate(id=paste(.chain, .iteration, sep="_"))%>%
  select(id, param, treatment, value)
allfits.seedling1<-allfits.seedling%>%
  mutate(id=paste(.chain, .iteration, sep="_"))%>%
  filter(.iteration>625)%>%
  select(id, param, treatment, value)
allfits.adult1<-allfits.adult%>%
  mutate(id=paste(.chain, .iteration, sep="_"))%>%
  filter(.iteration>625)%>%
  select(id, param, treatment, value)

lambdaAa<- subset(allfits.annuals1, param=="lam"&treatment=="ambient")
lambdaAw<- subset(allfits.annuals1, param=="lam"&treatment=="warmed")
alphaAAa<- subset(allfits.annuals1, param=="alphaAA"&treatment=="ambient")
alphaAAw<- subset(allfits.annuals1, param=="alphaAA"&treatment=="warmed")
alphaAPa<- subset(allfits.annuals1, param=="alphaAP"&treatment=="ambient")
alphaAPw<- subset(allfits.annuals1, param=="alphaAP"&treatment=="warmed")

lambdaPa<- subset(allfits.adult1, param=="lam"&treatment=="ambient")
lambdaPw<- subset(allfits.adult1, param=="lam"&treatment=="warmed")
alphaPAa<- subset(allfits.adult1, param=="alphaPA"&treatment=="ambient")
alphaPAw<-  subset(allfits.adult1, param=="alphaPA"&treatment=="warmed")
alphaPPa<-  subset(allfits.adult1, param=="alphaPP"&treatment=="ambient")
alphaPPw<- subset(allfits.adult1, param=="alphaPP"&treatment=="warmed")

lambdaSa<- subset(allfits.seedling1, param=="lam"&treatment=="ambient")
lambdaSw<- subset(allfits.seedling1, param=="lam"&treatment=="warmed")
alphaSAa<- subset(allfits.seedling1, param=="alphaSA"&treatment=="ambient")
alphaSAw<- subset(allfits.seedling1, param=="alphaSA"&treatment=="warmed")
alphaSSa<- subset(allfits.seedling1, param=="alphaSS"&treatment=="ambient")
alphaSSw<- subset(allfits.seedling1, param=="alphaSS"&treatment=="warmed")
alphaSPa<- subset(allfits.seedling1, param=="alphaSP"&treatment=="ambient")
alphaSPw<- subset(allfits.seedling1, param=="alphaSP"&treatment=="warmed")


# ------------------------------------------------------------------------------------
## Set germination and survival fractions from the literature
sa <-.11# .11 # ghersa 1984 # costa maia 2009 minus 4% permanently dormant
ga <- .89 #ghersa 1984
ga <- 100 #gundel 2007, gundel 2006
ga <- .97 #lin 2018

sp <- 1 # 100% adults survive observed in field - x in literature
sp <-.975 # or 2.5% fiegner 2007
ss <- 0
gs <- .63 # fiegner 2007
gs <- .90 # schmidt 1998
gs <- .31 # maret 2005
gs <- .50 # mackin 2021

# ------------------------------------------------------------------------------------
# Functions for use in coexistence calcualtions

#equilibriums 

# Determine equilibrium conditions for lolium seeds - to run for only a single timestep
annual.equilibrium <- function (N0a, sa, ga, alphaAA, lambdaA) { #number of annuals, seed survival annuals, germination annuals, self alpha, annual lambda
  Na <- sa*(1-ga)*N0a + N0a*lambdaA/(1+alphaAA*N0a) #predicted lambda is scaled, match scaling here 100, 30??
  return(Na)
}

# Determine equilibrium conditions for festuca adults
adult.equilibrium <- function (N0p, N0s, sp, alphaSS, alphaSP, lambdaS) { # n adults, n seedlings, self seedling alpha, adult on seedling alpha, seedling survival to adult
  Np <- sp*N0p + N0s*lambdaS/(1 + alphaSS*N0s + alphaSP*N0p) 
  return(Np)
}

# Determine equilibrium conditions for festuca seeds
seedling.equilibrium <- function (N0p, N0s, ss, gs, alphaPP, lambdaP) { # n adult n seedling, seedling survival, germination seedling, self and adult alphas, adult seed production
  Ns <- ss*(1-gs)*N0s + N0p*(lambdaP)/(1 + alphaPP*N0p)
  return(Ns)
}


#competition growth:
# annuals
annual.invade <- function  (N0a, N0p, N0s, sa, ga, alphaAA, alphaAP, lambdaA){
  Na <- sa*(1-ga)*N0a + N0a*lambdaA/(1+alphaAP*N0p+alphaAA*N0a)
  return(Na)
}

# seedlings
seedling.resident <- function (N0p, N0a, N0s, ss, gs, alphaPP, alphaPA, lambdaP) {
  Ns <- ss*(1-gs)*N0s + N0p*(lambdaP)/(1 + alphaPP*N0p+ alphaPA*N0a)  
  return(Ns)
}

# adults
adult.resident <- function (N0p, N0s, N0a, sp, alphaSS, alphaSP, alphaSA,lambdaS) {
  Np <- sp*N0p + N0s*lambdaS/(1 + alphaSS*N0s + alphaSP*N0p+ alphaSA*N0a) 
  return(Np)
}

# ------------------------------------------------------------------------------------
# run equilibriums (MEAN)

time <- 1:300
tmax=max(time)

# annuals, mean parameters
a <- rep(0, length(time))

N0a=1

#2020 ambient
eq.annuals.a20 <- tibble(time, a)
eq.annuals.a20[1,2] = as.numeric(N0a)

for (t in 1:tmax) {
  eq.annuals.a20[t+1,2] <- annual.equilibrium(eq.annuals.a20[t,2], sa, ga, alphaAAam20, lambdaAam20)
}
#2021 ambient
eq.annuals.a21 <- tibble(time, a)
eq.annuals.a21[1,2] = as.numeric(N0a)

for (t in 1:tmax) {
  eq.annuals.a21[t+1,2] <- annual.equilibrium(eq.annuals.a21[t,2], sa, ga, alphaAAam21, lambdaAam21)
}

#2020 warmed
eq.annuals.w20<- tibble(time, a)
eq.annuals.w20[1,2] <- N0a

for (t in 1:tmax) {
  eq.annuals.w20[t+1,2] <- annual.equilibrium(eq.annuals.w20[t,2], sa, ga, alphaAAwm20, lambdaAwm20)
}

#2021 warmed
eq.annuals.w21<- tibble(time, a)
eq.annuals.w21[1,2] <- N0a

for (t in 1:tmax) {
  eq.annuals.w21[t+1,2] <- annual.equilibrium(eq.annuals.w21[t,2], sa, ga, alphaAAwm21, lambdaAwm21)
}


#perennials

N0s=1
N0p=1
s <- rep(0, length(time))
p <- rep(0, length(time))

#2020 amb
eq.perennials.a20<-tibble(time, p, s)
eq.perennials.a20[1,3] <- N0s
eq.perennials.a20[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.a20[t+1, 2] <- adult.equilibrium(eq.perennials.a20[t, 2], eq.perennials.a20[t, 3], sp, alphaSSam20, alphaSPam20, survivalSam20)
  eq.perennials.a20[t+1, 3] <- seedling.equilibrium(eq.perennials.a20[t, 2], eq.perennials.a20[t, 3], ss, gs, alphaPPam20, lambdaPam20)
}

#2021 amb
eq.perennials.a21<-tibble(time, p, s)
eq.perennials.a21[1,3] <- N0s
eq.perennials.a21[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.a21[t+1, 2] <- adult.equilibrium(eq.perennials.a21[t, 2], eq.perennials.a21[t, 3], sp, alphaSSam21, alphaSPam21, survivalSam21)
  eq.perennials.a21[t+1, 3] <- seedling.equilibrium(eq.perennials.a21[t, 2], eq.perennials.a21[t, 3], ss, gs, alphaPPam21, lambdaPam21)
}


#2020 warm
eq.perennials.w20<-tibble(time, p, s)
eq.perennials.w20[1,3] <- N0s
eq.perennials.w20[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.w20[t+1, 2] <- adult.equilibrium(eq.perennials.w20[t, 2], eq.perennials.w20[t, 3], sp, alphaSSwm20, alphaSPam20, survivalSwm20)
  eq.perennials.w20[t+1, 3] <- seedling.equilibrium(eq.perennials.w20[t, 2], eq.perennials.w20[t, 3], ss, gs, alphaPPwm20, lambdaPwm20)
}

#2021 warm
eq.perennials.w21<-tibble(time, p, s)
eq.perennials.w21[1,3] <- N0s
eq.perennials.w21[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.w21[t+1, 2] <- adult.equilibrium(eq.perennials.w21[t, 2], eq.perennials.w21[t, 3], sp, alphaSSwm21, alphaSPam21, survivalSwm21)
  eq.perennials.w21[t+1, 3] <- seedling.equilibrium(eq.perennials.w21[t, 2], eq.perennials.w21[t, 3], ss, gs, alphaPPwm21, lambdaPwm21)
}


# check outputs
annual.eq<-rbind(mutate(eq.annuals.a20, temp="amb", year=20), mutate(eq.annuals.w20, temp="warm", year=20), mutate(eq.annuals.a21, temp="amb", year=21), mutate(eq.annuals.w21, temp="warm", year=21))
perennial.eq<-rbind(mutate(eq.perennials.a20, temp="amb", year=20), mutate(eq.perennials.w20, temp="warm", year=20), mutate(eq.perennials.a21, temp="amb", year=21), mutate(eq.perennials.w21, temp="warm", year=21))%>%
  gather("type", "density", p, s)%>%
  mutate(type=ifelse(type=="p", "adult", "seedling"))

eqp1<-ggplot(annual.eq, aes(time, a, color=temp))+ geom_line()+ylab("annual seed density")+  scale_colour_manual(values = c("dodgerblue", "darkred")) +facet_wrap(~year)
eqp2<-ggplot(perennial.eq, aes(time, (density), color=temp, linetype=as.factor(type) ))+ geom_line()+ylab("perennial density")+scale_y_continuous(trans='log10', labels = scales::comma)+  scale_colour_manual(values = c("dodgerblue", "darkred"))+facet_wrap(~year)

ggarrange(eqp1, eqp2)

# ------------------------------------------------------------------------------------
# run equilibriums add 100 chains for each to visualization
warmed.annual.chains<-select(lambdaAw, 1:3)%>%mutate(eq=0)
ambient.annual.chains<-select(lambdaAa, 1:3)%>%mutate(eq=0)

warmed.seedling.chains<-select(lambdaSw, 1:3)%>%mutate(eq=0)
ambient.seedling.chains<-select(lambdaSa, 1:3)%>%mutate(eq=0)
  
warmed.adult.chains<-select(lambdaPw, 1:3)%>%mutate(eq=0)
ambient.adult.chains<-select(lambdaPa, 1:3)%>%mutate(eq=0)

#annuals
#for (i in 1:nrow(warmed.annual.chains)) {
##  a <- rep(0, length(time))
#  N0a=1
#  eq.annuals.a <- tibble(time, a)
#  eq.annuals.a[1,2] = as.numeric(N0a)
#  
#  for (t in 1:tmax) {
#    eq.annuals.a[t+1,2] <- annual.equilibrium(eq.annuals.a[t,2], sa, ga, alphaAAa[i, 4], lambdaAa[i, 4])
#  }
#  ambient.annual.chains[i, 4]<-eq.annuals.a[tmax, 2]
#  
#  eq.annuals.w<- tibble(time, a)
#  eq.annuals.w[1,2] <- N0a
#  
#  for (t in 1:tmax) {
#    eq.annuals.w[t+1,2] <- annual.equilibrium(eq.annuals.w[t,2], sa, ga, alphaAAw[i, 4], lambdaAw[i, 4])
 # }
#  warmed.annual.chains[i, 4]<-eq.annuals.w[tmax, 2]
#}


N0a=1
time <- 1:10
tmax=max(time)

warmed.annual.chains1<-select(lambdaAw, 1:3)
warmed.annual.chains1 <- cbind(warmed.annual.chains1, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0a)

ambient.annual.chains1<-select(lambdaAa, 1:3)
ambient.annual.chains1 <- cbind(ambient.annual.chains1, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0a)


for (i in 1:100) {
  for (t in 1:tmax) {
    ambient.annual.chains1[i, t+4] <- annual.equilibrium(ambient.annual.chains1[i, t+3], sa, ga, alphaAAa[i, 4], lambdaAa[i, 4])
  }
  for (t in 1:tmax) {
    warmed.annual.chains1[i, t+4] <- annual.equilibrium(warmed.annual.chains1[i, t+3], sa, ga, alphaAAw[i, 4], lambdaAw[i, 4])
  }
}


#perennials

tmax=300
warmed.perennial.seedling.chains1<-select(lambdaSw, 1:3)
warmed.perennial.seedling.chains1 <- cbind(warmed.perennial.seedling.chains1, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0s)

ambient.perennial.seedling.chains1<-select(lambdaSa, 1:3)
ambient.perennial.seedling.chains1 <- cbind(ambient.perennial.seedling.chains1, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0s)

warmed.perennial.adult.chains1<-select(lambdaPw, 1:3)
warmed.perennial.adult.chains1 <- cbind(warmed.perennial.adult.chains1, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0p)

ambient.perennial.adult.chains1<-select(lambdaPa, 1:3)
ambient.perennial.adult.chains1 <- cbind(ambient.perennial.adult.chains1, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0p)


for (i in 1:50) {  # this value is always 2000, no need to change it (just to get this to run I do)
  for (t in 1:tmax) {
    ambient.perennial.adult.chains1[i, t+4] <- adult.equilibrium(ambient.perennial.adult.chains1[i, t+3], ambient.perennial.seedling.chains1[i, t+3], sp, alphaSSa[i, 4], alphaSPa[i, 4], lambdaSa[i, 4])
    ambient.perennial.seedling.chains1[i, t+4] <- seedling.equilibrium(ambient.perennial.adult.chains1[i, t+3], ambient.perennial.seedling.chains1[i, t+3], ss, gs, alphaPPa[i, 4], lambdaPa[i, 4])
  }
  for (t in 1:tmax) {
    warmed.perennial.adult.chains1[i, t+4] <- adult.equilibrium(warmed.perennial.adult.chains1[i, t+3], warmed.perennial.seedling.chains1[i, t+3], sp, alphaSSw[i, 4], alphaSPw[i, 4], lambdaSw[i, 4])
    warmed.perennial.seedling.chains1[i, t+4] <- seedling.equilibrium(warmed.perennial.adult.chains1[i, t+3], warmed.perennial.seedling.chains1[i, t+3], ss, gs, alphaPPw[i, 4], lambdaPw[i, 4])
  }
}



# visualize chains
annual.chains<-  rbind(ambient.annual.chains1[1:100,],              warmed.annual.chains1[1:100,])%>%
  gather("time", "population", 4:14)%>% #this only runs to 14 because it only runs for 10 years
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(id, treatment))
seedling.chains<-rbind(ambient.perennial.seedling.chains1[1:50,],  warmed.perennial.seedling.chains1[1:50,])%>%
  gather("time", "population", 4:304) %>% # this runs to 304 because it takes 300 years to equilibrate
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(id, treatment))
adult.chains<-   rbind(ambient.perennial.adult.chains1[1:50,],     warmed.perennial.adult.chains1[1:50,])%>%
  gather("time", "population", 4:304)%>%
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(id, treatment))

eqannual<-ggplot(subset(annual.chains, population>0), aes(time, (population)))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("annual seed equilibrium density")+
  geom_line(data=subset(annual.eq, time<12), size=1.5, aes(x=time, y=a, color=temp))+
  scale_colour_manual(values = c("dodgerblue", "dodgerblue", "darkred", "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)
eqseedling<-ggplot(subset(seedling.chains, population>0), aes(time, population))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("perennial seed equilibrium density")+
  geom_line(data=subset(perennial.eq, type=="seedling"), size=1.5, aes(x=time, y=density, color=temp))+
  scale_colour_manual(values = c("dodgerblue","dodgerblue", "darkred", "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)
eqadult<-ggplot(subset(adult.chains, population>0), aes(time, (population)))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("perennial adult equilibrium density")+
  geom_line(data=subset(perennial.eq, type=="adult"), size=1.5, aes(x=time, y=density, color=temp))+
  scale_colour_manual(values = c("dodgerblue","dodgerblue", "darkred", "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)

ggarrange(eqannual, eqseedling, eqadult, nrow=1, ncol=3, common.legend = T)

# ------------------------------------------------------------------------------------
# invasion
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))
# annual into perennials 
#ambient
annual.invasion.a <- tibble(time, a, p, s)
annual.invasion.a[1,2] <- 1 #annuals
annual.invasion.a[1,3] <- eq.perennials.a[300,2] #adults at eq
annual.invasion.a[1,4] <- eq.perennials.a[300,3] #seedlings at eq


for (t in 1:tmax) {
  annual.invasion.a[t+1, 2] <- annual.invade(N0a=annual.invasion.a[t,2],
                                             N0p=annual.invasion.a[t,3],N0s=annual.invasion.a[t,4],
                                    sa, ga, alphaAAam, alphaAPam, lambdaAam)
  annual.invasion.a[t+1, 4] <- seedling.resident(N0a=annual.invasion.a[t,2],
                                      N0p=annual.invasion.a[t,3],N0s=annual.invasion.a[t,4],
                                      ss, gs, alphaPPam, alphaPAam, lambdaPam)
  annual.invasion.a[t+1, 3] <- adult.resident(N0a=annual.invasion.a[t,2],
                                      N0p=annual.invasion.a[t,3],
                                      N0s=annual.invasion.a[t,4], 
                                      sp, alphaSSam, alphaSPam, alphaSAam, lambdaSam)
  }

#warmed
annual.invasion.w <- tibble(time, a, p, s)
annual.invasion.w[1,2] <- 1
annual.invasion.w[1,3] <- eq.perennials.w[300,2] #adults at eq
annual.invasion.w[1,4] <- eq.perennials.w[300,3] #seedlings at eq

for (t in 1:tmax) {
  annual.invasion.w[t+1, 2] <- annual.invade(N0a=annual.invasion.w[t,2],
                                             N0p=annual.invasion.w[t,3],N0s=annual.invasion.w[t,4],
                                             sa, ga, alphaAAwm, alphaAPwm, lambdaAwm)
  annual.invasion.w[t+1, 4] <- seedling.resident(N0a=annual.invasion.w[t,2],
                                                 N0p=annual.invasion.w[t,3],N0s=annual.invasion.w[t,4],
                                                 ss, gs, alphaPPwm, alphaPAwm, lambdaPwm)
  annual.invasion.w[t+1, 3] <- adult.resident(N0a=annual.invasion.w[t,2],
                                              N0p=annual.invasion.w[t,3],
                                              N0s=annual.invasion.w[t,4], 
                                              sp, alphaSSwm, alphaSPwm, alphaSAwm, lambdaSwm)
}


# perennials into annuals

#ambient
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))

perennial.invasion.a <- tibble(time, a, p, s)
perennial.invasion.a[1,2] <-  eq.annuals.a[tmax,2] #annuals at eq
perennial.invasion.a[1,3] <- .5
perennial.invasion.a[1,4] <- .5


for (t in 1:tmax) {
  perennial.invasion.a[t+1, 2] <- annual.invade(N0a=perennial.invasion.a[t,2],
                                                N0p=perennial.invasion.a[t,3],N0s=perennial.invasion.a[t,4],
                                                sa, ga, alphaAAam, alphaAPam, lambdaAam)
  perennial.invasion.a[t+1, 4] <- seedling.resident(N0a=perennial.invasion.a[t,2],
                                                    N0p=perennial.invasion.a[t,3],N0s=perennial.invasion.a[t,4],
                                                    ss, gs, alphaPPam, alphaPAam, lambdaPam)
  perennial.invasion.a[t+1, 3] <- adult.resident(N0a=perennial.invasion.a[t,2],
                                                N0p=perennial.invasion.a[t,3],
                                                N0s=perennial.invasion.a[t,4], 
                                                sp, alphaSSam, alphaSPam, alphaSAam, lambdaSam)
}

#warmed
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))
perennial.invasion.w <- tibble(time, a, p, s)
perennial.invasion.w[1,2] <-  eq.annuals.w[tmax,2] #annuals at eq
perennial.invasion.w[1,3] <- .5
perennial.invasion.w[1,4] <- .5


for (t in 1:tmax) {
  perennial.invasion.w[t+1, 2] <- annual.invade(N0a=perennial.invasion.w[t,2],
                                                N0p=perennial.invasion.w[t,3],N0s=perennial.invasion.w[t,4],
                                                sa, ga, alphaAAwm, alphaAPwm, lambdaAwm)
  perennial.invasion.w[t+1, 4] <- seedling.resident(N0a=perennial.invasion.w[t,2],
                                                    N0p=perennial.invasion.w[t,3],N0s=perennial.invasion.w[t,4],
                                                    ss, gs, alphaPPwm, alphaPAwm, lambdaPwm)
  perennial.invasion.w[t+1, 3] <- adult.resident(N0a=perennial.invasion.w[t,2],
                                                 N0p=perennial.invasion.w[t,3],
                                                 N0s=perennial.invasion.w[t,4], 
                                                 sp, alphaSSwm, alphaSPwm, alphaSAwm, lambdaSwm)
}
  

#visualize

# check outputs
annual.inv<-rbind(mutate(annual.invasion.a, temp="amb"), mutate(annual.invasion.w, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))


perennial.inv<-rbind(mutate(perennial.invasion.a, temp="amb"), mutate(perennial.invasion.w, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))

inv1<-ggplot(subset(annual.inv, time<16), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
geom_line(size=1.25)+ylab("density")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))
inv2<-ggplot(subset(perennial.inv, time<16), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
geom_line(size=1.25)+ylab("density")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))

ggarrange(inv1, inv2, nrow=2, ncol=1)



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

perennialGRWRa = log(.5*(sp + (4* lambdaPam*5000/(1 + alphaPAam*perennial.invasion.a[1,2]) * lambdaSam/(1 + alphaSAam*perennial.invasion.a[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="ambient")
perennialGRWRw = log(.5*(sp + (4* lambdaPwm*5000/(1 + alphaPAwm*perennial.invasion.w[1,2]) * lambdaSwm/(1 + alphaSAwm*perennial.invasion.w[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="warmed")


allGRWR<-rbind(annualGRWRa, annualGRWRw, perennialGRWRa, perennialGRWRw)%>%
  mutate(GRWR=a)%>%
  select(-a)

#plot of mean paremter GRWR
ggplot(allGRWR, aes(x=trt, y=GRWR)) +geom_bar(stat="identity", aes(fill=trt)) +facet_wrap(~invader, scales="free") +xlab("") +ylab("Invader Growth Rate (r) When Rare") + scale_fill_manual(values = c("dodgerblue", "darkred"))


#GRWR with mean eq but all the parameters for invasion from there


#annual GRWR (chains)

GRWRchainAa<-select(lambdaAa[1:2000,], 1,3)
GRWRchainAa <- cbind(GRWRchainAa, t(rep(0,tmax-3)))%>%
  mutate(`1`=N0a)

GRWRchainAw<-select(lambdaAw[1:2000,], 1,3)
GRWRchainAw <- cbind(GRWRchainAw, t(rep(0,tmax-3)))%>%
  mutate(`1`=N0a)


#grwr algebra from mordecai, transitions adjusted for my model
for (i in 1:2000) {  # take 100 chains for invasion
  GRWRchainAa[i,4]<-log(lambdaAa[i,4]*100/(1 + alphaAPa[i,4]*annual.invasion.a[1,3])) 
  GRWRchainAw[i,4]<-log(lambdaAw[i,4]*100/(1 + alphaAPw[i,4]*annual.invasion.w[1,3])) 
}


GRWRchainA<-rbind(GRWRchainAa,GRWRchainAw)%>%
  mutate(GRWR=`2`)%>%
  select(id, treatment, GRWR)%>%
  mutate(type="annual")

GRWRchainAsum<-GRWRchainA%>%
  group_by(treatment,type)%>%
  summarize(meanGRWR=mean(GRWR), seGRWR=calcSE(GRWR))

annualGRWRchain<-ggplot(GRWRchainAsum, aes(x=treatment, y=meanGRWR))+
  #  scale_y_continuous(trans='log10')+
  ylim(-.03, 4.5)+
  geom_jitter(data=GRWRchainA, aes(x=treatment, y=GRWR), width=.2,size=.25)+geom_hline(yintercept = 0)+
  geom_point(stat="identity", position="dodge",aes(color=treatment), size=3)+
  scale_color_manual(values = c("dodgerblue","darkred"))+facet_wrap(~type)+ylab("")+
  geom_errorbar (aes(ymin=meanGRWR-seGRWR, ymax=meanGRWR+seGRWR, color=treatment), width=.2)+theme(text=element_text(size=16))


#perennials GRWR chains
GRWRchainPa<-select(lambdaPa[1:2000,], 1,3)
GRWRchainPa <- cbind(GRWRchainPa, t(rep(0,tmax-3)))%>%
  mutate(`1`=N0p)

GRWRchainPw<-select(lambdaPw[1:2000,], 1,3)
GRWRchainPw <- cbind(GRWRchainPw, t(rep(0,tmax-3)))%>%
  mutate(`1`=N0p)

#grwr algebra from mordecai, transitions adjusted for my model
for (i in 1:2000) {  # take 100 chains for invasion
  GRWRchainPa[i,4]<-log(.5*(sp + (4* lambdaPa[i,4]*5000/(1 + alphaPAa[i,4]*perennial.invasion.a[1,2]) * lambdaSa[i,4]/(1 + alphaSAa[i,4]*perennial.invasion.a[1,2]) + sp^2)^0.5))
  GRWRchainPw[i,4]<-log(.5*(sp + (4* lambdaPw[i,4]*5000/(1 + alphaPAw[i,4]*perennial.invasion.a[1,2]) * lambdaSw[i,4]/(1 + alphaSAw[i,4]*perennial.invasion.a[1,2]) + sp^2)^0.5))
  
}

GRWRchainP<-rbind(GRWRchainPa,GRWRchainPw)%>%
  mutate(GRWR=`2`)%>%
  select(id, treatment, GRWR)%>%
  mutate(type="perennial")

GRWRchainPsum<-GRWRchainP%>%
  group_by(treatment,type)%>%
  summarize(meanGRWR=mean(GRWR), seGRWR=calcSE(GRWR))

perennialGRWRchain<-ggplot(GRWRchainPsum, aes(x=treatment, y=meanGRWR))+
#  scale_y_continuous(trans='log10')+  
  ylim(-.03, 4.5)+
  geom_jitter(data=GRWRchainP, aes(x=treatment, y=GRWR), width=.2,size=.25)+geom_hline(yintercept = 0)+
  geom_point(stat="identity", position="dodge",aes(color=treatment), size=3)+
  scale_color_manual(values = c("dodgerblue","darkred"))+facet_wrap(~type)+ylab("GRWR")+
  geom_errorbar (aes(ymin=meanGRWR-seGRWR, ymax=meanGRWR+seGRWR, color=treatment), width=.2)+theme(text=element_text(size=16))
  

#visualize
ggarrange(perennialGRWRchain, annualGRWRchain, common.legend=T)

#compare with mean
ggplot(allGRWR, aes(x=trt, y=GRWR)) +geom_point(stat="identity", aes(color=trt)) +
  facet_wrap(~invader) + geom_hline(yintercept=0) +xlab("") +ylab("Invader Growth Rate (r) When Rare") + 
  scale_color_manual(values = c("dodgerblue", "darkred"))  +ylim(-.03, 4.5)+theme(text=element_text(size=16))


for (i in 1:100) {  # take 100 chains for invasion
  aPc[i, t+4] <- adult.equilibrium(aPc[i, t+3], aPs[i, t+3], sp, alphaSSa[i, 4], alphaSPa[i, 4], lambdaSa[i, 4])
  aSc[i, t+4] <- seedling.equilibrium(aPc[i, t+3], aPs[i, t+3], ss, gs, alphaPPa[i, 4], lambdaPa[i, 4])
}
for (t in 1:tmax) {
  wPc[i, t+4] <- adult.equilibrium(wPc[i, t+3], wSc[i, t+3], sp, alphaSSw[i, 4], alphaSPw[i, 4], lambdaSw[i, 4])
  wSc[i, t+4] <- seedling.equilibrium(wPc[i, t+3], wSc[i, t+3], ss, gs, alphaPPw[i, 4], lambdaPw[i, 4])
}

}





#method 2, need to create a 2x2 transition matrix with (seed survival term, seed production term )
#                                                      (seeds maturing term, adults surviving term)
library(popbio)

eigen.analysis(perennial.invasion.a[1:2,3:4])
eigen.analysis(perennial.invasion.w[1:2,3:4])


#### ----# ------------------------------------------------------------------------------------
# invasion sandbox
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))
# annual into perennials 
#ambient
annual.invasion.a <- tibble(time, a, p, s)
annual.invasion.a[1,2] <- 1 #annuals
annual.invasion.a[1,3] <- eq.perennials.a[300,2] #adults at eq
annual.invasion.a[1,4] <- eq.perennials.a[300,3] #seedlings at eq


for (t in 1:tmax) {
  annual.invasion.a[t+1, 2] <- annual.invade(N0a=annual.invasion.a[t,2],
                                             N0p=annual.invasion.a[t,3],N0s=annual.invasion.a[t,4],
                                             sa, ga, alphaAAam, alphaAPam, lambdaAam)
  annual.invasion.a[t+1, 4] <- seedling.resident(N0a=annual.invasion.a[t,2],
                                                 N0p=annual.invasion.a[t,3],N0s=annual.invasion.a[t,4],
                                                 ss, gs, alphaPPam, alphaPAam, lambdaPam)
  annual.invasion.a[t+1, 3] <- adult.resident(N0a=annual.invasion.a[t,2],
                                              N0p=annual.invasion.a[t,3],
                                              N0s=annual.invasion.a[t,4], 
                                              .75, alphaSSam, alphaSPam, alphaSAam, lambdaSam)
}

#warmed
annual.invasion.w <- tibble(time, a, p, s)
annual.invasion.w[1,2] <- 1
annual.invasion.w[1,3] <- eq.perennials.w[300,2] #adults at eq
annual.invasion.w[1,4] <- eq.perennials.w[300,3] #seedlings at eq

for (t in 1:tmax) {
  annual.invasion.w[t+1, 2] <- annual.invade(N0a=annual.invasion.w[t,2],
                                             N0p=annual.invasion.w[t,3],N0s=annual.invasion.w[t,4],
                                             sa, ga, alphaAAwm, alphaAPwm, lambdaAwm)
  annual.invasion.w[t+1, 4] <- seedling.resident(N0a=annual.invasion.w[t,2],
                                                 N0p=annual.invasion.w[t,3],N0s=annual.invasion.w[t,4],
                                                 ss, gs, alphaPPwm, alphaPAwm, lambdaPwm)
  annual.invasion.w[t+1, 3] <- adult.resident(N0a=annual.invasion.w[t,2],
                                              N0p=annual.invasion.w[t,3],
                                              N0s=annual.invasion.w[t,4], 
                                              .75, alphaSSwm, alphaSPwm, alphaSAwm, lambdaSwm)
}


# perennials into annuals

#ambient
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))

perennial.invasion.a <- tibble(time, a, p, s)
perennial.invasion.a[1,2] <-  eq.annuals.a[tmax,2] #annuals at eq
perennial.invasion.a[1,3] <- .5
perennial.invasion.a[1,4] <- .5


for (t in 1:tmax) {
  perennial.invasion.a[t+1, 2] <- annual.invade(N0a=perennial.invasion.a[t,2],
                                                N0p=perennial.invasion.a[t,3],N0s=perennial.invasion.a[t,4],
                                                sa, ga, alphaAAam, alphaAPam, lambdaAam)
  perennial.invasion.a[t+1, 4] <- seedling.resident(N0a=perennial.invasion.a[t,2],
                                                    N0p=perennial.invasion.a[t,3],N0s=perennial.invasion.a[t,4],
                                                    ss, gs, alphaPPam, alphaPAam, lambdaPam)
  perennial.invasion.a[t+1, 3] <- adult.resident(N0a=perennial.invasion.a[t,2],
                                                 N0p=perennial.invasion.a[t,3],
                                                 N0s=perennial.invasion.a[t,4], 
                                                 .75, alphaSSam, alphaSPam, alphaSAam, lambdaSam)
}

#warmed
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))
perennial.invasion.w <- tibble(time, a, p, s)
perennial.invasion.w[1,2] <-  eq.annuals.w[tmax,2] #annuals at eq
perennial.invasion.w[1,3] <- .5
perennial.invasion.w[1,4] <- .5


for (t in 1:tmax) {
  perennial.invasion.w[t+1, 2] <- annual.invade(N0a=perennial.invasion.w[t,2],
                                                N0p=perennial.invasion.w[t,3],N0s=perennial.invasion.w[t,4],
                                                sa, ga, alphaAAwm, alphaAPwm, lambdaAwm)
  perennial.invasion.w[t+1, 4] <- seedling.resident(N0a=perennial.invasion.w[t,2],
                                                    N0p=perennial.invasion.w[t,3],N0s=perennial.invasion.w[t,4],
                                                    ss, gs, alphaPPwm, alphaPAwm, lambdaPwm)
  perennial.invasion.w[t+1, 3] <- adult.resident(N0a=perennial.invasion.w[t,2],
                                                 N0p=perennial.invasion.w[t,3],
                                                 N0s=perennial.invasion.w[t,4], 
                                                 .75, alphaSSwm, alphaSPwm, alphaSAwm, lambdaSwm)
}


#visualize

# check outputs
annual.inv<-rbind(mutate(annual.invasion.a, temp="amb"), mutate(annual.invasion.w, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))


perennial.inv<-rbind(mutate(perennial.invasion.a, temp="amb"), mutate(perennial.invasion.w, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))

inv1<-ggplot(subset(annual.inv, time<30), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
  geom_line(size=1.25)+ylab("density")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))
inv2<-ggplot(subset(perennial.inv, time<30), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
  geom_line(size=1.25)+ylab("density")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))

ggarrange(inv1, inv2, nrow=2, ncol=1)


