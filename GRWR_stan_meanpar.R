# mean parameters ----

# new version. annual and seedling parameters will be averaged across years
# only allow adult perennial to vary

#annual
lambdaAam<- (annuals_estimates2020$exp[5]+annuals_estimates2021$exp[5])/2 
lambdaAwm<- (annuals_estimates2020$exp[6]+annuals_estimates2021$exp[6])/2 

alphaAAam<- ( annuals_estimates2020$exp[1]+ annuals_estimates2021$exp[1])/2
alphaAAwm<- ( annuals_estimates2020$exp[2]+annuals_estimates2021$exp[2])/2

alphaPAam<-( adults_estimates2020$exp[1]+adults_estimates2021$exp[1])/2
alphaPAwm<- ( adults_estimates2020$exp[2]+ adults_estimates2021$exp[2])/2

alphaSAam<- ( seedlings_estimates2020$exp[1]+ seedlings_estimates2021$exp[1])/2
alphaSAwm<- ( seedlings_estimates2020$exp[2]+ seedlings_estimates2021$exp[2])/2

#seedling
survivalSam<- (seedlings_estimates2020$exp[7]+ seedlings_estimates2021$exp[7])/2
survivalSwm<- (seedlings_estimates2020$exp[8]+ seedlings_estimates2021$exp[8])/2

alphaSSam<- (seedlings_estimates2020$exp[5]+seedlings_estimates2021$exp[5])/2
alphaSSwm<- (seedlings_estimates2020$exp[6]+ seedlings_estimates2021$exp[6])/2

#adult
alphaAPam20<- annuals_estimates2020$exp[3]
alphaAPwm20<- annuals_estimates2020$exp[4]
alphaAPam21<- annuals_estimates2021$exp[3]
alphaAPwm21<- annuals_estimates2021$exp[4]

lambdaPam20<- adults_estimates2020$exp[5]
lambdaPwm20<- adults_estimates2020$exp[6]
lambdaPam21<- adults_estimates2021$exp[5]
lambdaPwm21<- adults_estimates2021$exp[6] 

alphaPPam20<-  adults_estimates2020$exp[3]
alphaPPwm20<- adults_estimates2020$exp[4]
alphaPPam21<-  adults_estimates2021$exp[3]
alphaPPwm21<- adults_estimates2021$exp[4]

alphaSPam20<- seedlings_estimates2020$exp[3]
alphaSPwm20<- seedlings_estimates2020$exp[4]
alphaSPam21<- seedlings_estimates2021$exp[3]
alphaSPwm21<- seedlings_estimates2021$exp[4]

# ------------------------------------------------------------------------------------
## Set germination and survival fractions from the literature
sa <-0 # .11 # ghersa 1984 # costa maia 2009 minus 4% permanently dormant
ga <- .89  #ghersa 1984
#ga <- 100 #gundel 2007, gundel 2006
#ga <- .97 #lin 2018

#sp <- 1 # 100% adults survive observed in field - x in literature
sp <-.975*.95 # or 2.5% fiegner 2007
ss <- 0.1
gs <- .63 # fiegner 2007
#gs <- .90 # schmidt 1998
#gs <- .31 # maret 2005
#gs <- .50 # mackin 2021

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

# annuals
a <- rep(0, length(time))

N0a=1

#2020 ambient
eq.annuals.a20 <- tibble(time, a)
eq.annuals.a20[1,2] = as.numeric(N0a)

for (t in 1:tmax) {
  eq.annuals.a20[t+1,2] <- annual.equilibrium(eq.annuals.a20[t,2], sa, ga, alphaAAam20, lambdaAam20)
}
#2021 ambient (no longer different)
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
eq.annuals.w21[1,2] <- N0a#

for (t in 1:tmax) {
  eq.annuals.w21[t+1,2] <- annual.equilibrium(eq.annuals.w21[t,2], sa, ga, alphaAAwm21, lambdaAwm21)
}

#mean warmed
eq.annuals.w<- tibble(time, a)
eq.annuals.w[1,2] <- N0a

for (t in 1:tmax) {
  eq.annuals.w[t+1,2] <- annual.equilibrium(eq.annuals.w[t,2], sa, ga, alphaAAwm, lambdaAwm)
}

#mean amb
eq.annuals.a<- tibble(time, a)
eq.annuals.a[1,2] <- N0a#

for (t in 1:tmax) {
  eq.annuals.a[t+1,2] <- annual.equilibrium(eq.annuals.a[t,2], sa, ga, alphaAAam, lambdaAam)
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
  eq.perennials.a20[t+1, 2] <- adult.equilibrium(eq.perennials.a20[t, 2], eq.perennials.a20[t, 3], sp, alphaSSam, alphaSPam20, survivalSam)
  eq.perennials.a20[t+1, 3] <- seedling.equilibrium(eq.perennials.a20[t, 2], eq.perennials.a20[t, 3], ss, gs, alphaPPam20, lambdaPam20)
}


#2021 amb
eq.perennials.a21<-tibble(time, p, s)
eq.perennials.a21[1,3] <- N0s
eq.perennials.a21[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.a21[t+1, 2] <- adult.equilibrium(eq.perennials.a21[t, 2], eq.perennials.a21[t, 3], sp, alphaSSam, alphaSPam21, survivalSam)
  eq.perennials.a21[t+1, 3] <- seedling.equilibrium(eq.perennials.a21[t, 2], eq.perennials.a21[t, 3], ss, gs, alphaPPam21, lambdaPam21)
}



#2020 warm
eq.perennials.w20<-tibble(time, p, s)
eq.perennials.w20[1,3] <- N0s
eq.perennials.w20[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.w20[t+1, 2] <- adult.equilibrium(eq.perennials.w20[t, 2], eq.perennials.w20[t, 3], sp, alphaSSwm, alphaSPam20, survivalSwm)
  eq.perennials.w20[t+1, 3] <- seedling.equilibrium(eq.perennials.w20[t, 2], eq.perennials.w20[t, 3], ss, gs, alphaPPwm20, lambdaPwm20)
}


#2021 warm
eq.perennials.w21<-tibble(time, p, s)
eq.perennials.w21[1,3] <- N0s
eq.perennials.w21[1,2] <- N0p

for (t in 1:300) {
  eq.perennials.w21[t+1, 2] <- adult.equilibrium(eq.perennials.w21[t, 2], eq.perennials.w21[t, 3], sp, alphaSSwm, alphaSPam21, survivalSwm)
  eq.perennials.w21[t+1, 3] <- seedling.equilibrium(eq.perennials.w21[t, 2], eq.perennials.w21[t, 3], ss, gs, alphaPPwm21, lambdaPwm21)
}


# check outputs
annual.eq<-rbind(mutate(eq.annuals.a20, temp="amb", year=20), mutate(eq.annuals.w20, temp="warm", year=20), mutate(eq.annuals.a21, temp="amb", year=21), mutate(eq.annuals.w21, temp="warm", year=21), mutate(eq.annuals.a, temp="amb", year=0), mutate(eq.annuals.w, temp="warm", year=0))
perennial.eq<-rbind(mutate(eq.perennials.a20, temp="amb", year=20), mutate(eq.perennials.w20, temp="warm", year=20), mutate(eq.perennials.a21, temp="amb", year=21), mutate(eq.perennials.w21, temp="warm", year=21))%>%
  gather("type", "density", p, s)%>%
  mutate(type=ifelse(type=="p", "adult", "seedling"))

eqp1<-ggplot(filter(annual.eq, time<300), aes(time, log(a), color=temp))+ geom_line()+ylab("annual seed density")+  scale_colour_manual(values = c("dodgerblue", "darkred")) +facet_wrap(~year)
eqp2<-ggplot(filter(perennial.eq, time<300), aes(time, (density), color=temp, linetype=as.factor(type) ))+ geom_line()+ylab("perennial density")+scale_y_continuous(trans='log10', labels = scales::comma)+  scale_colour_manual(values = c("dodgerblue", "darkred"))+facet_wrap(~year)

ggarrange(eqp1, eqp2)


# ------------------------------------------------------------------------------------
# invasion
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))
# annual into perennials 
#2020
#ambient
annual.invasion.a20 <- tibble(time, a, p, s)
annual.invasion.a20[1,2] <- 1 #annuals
annual.invasion.a20[1,3] <- eq.perennials.a20[300,2] #adults at eq
annual.invasion.a20[1,4] <- eq.perennials.a20[300,3] #seedlings at eq


for (t in 1:tmax) {
  annual.invasion.a20[t+1, 2] <- annual.invade(N0a=annual.invasion.a20[t,2],
                                             N0p=annual.invasion.a20[t,3],N0s=annual.invasion.a20[t,4],
                                    sa, ga, alphaAAam, alphaAPam20, lambdaAam)
  annual.invasion.a20[t+1, 4] <- seedling.resident(N0a=annual.invasion.a20[t,2],
                                      N0p=annual.invasion.a20[t,3],N0s=annual.invasion.a20[t,4],
                                      ss, gs, alphaPPam20, alphaPAam, lambdaPam20)
  annual.invasion.a20[t+1, 3] <- adult.resident(N0a=annual.invasion.a20[t,2],
                                      N0p=annual.invasion.a20[t,3],
                                      N0s=annual.invasion.a20[t,4], 
                                      sp, alphaSSam, alphaSPam20, alphaSAam, survivalSam)
  }

#warmed
annual.invasion.w20 <- tibble(time, a, p, s)
annual.invasion.w20[1,2] <- 1
annual.invasion.w20[1,3] <- eq.perennials.w20[300,2] #adults at eq
annual.invasion.w20[1,4] <- eq.perennials.w20[300,3] #seedlings at eq

for (t in 1:tmax) {
  annual.invasion.w20[t+1, 2] <- annual.invade(N0a=annual.invasion.w20[t,2],
                                             N0p=annual.invasion.w20[t,3],N0s=annual.invasion.w20[t,4],
                                             sa, ga, alphaAAwm, alphaAPwm20, lambdaAwm)
  annual.invasion.w20[t+1, 4] <- seedling.resident(N0a=annual.invasion.w20[t,2],
                                                 N0p=annual.invasion.w20[t,3],N0s=annual.invasion.w20[t,4],
                                                 ss, gs, alphaPPwm20, alphaPAwm, lambdaPwm20)
  annual.invasion.w20[t+1, 3] <- adult.resident(N0a=annual.invasion.w20[t,2],
                                              N0p=annual.invasion.w20[t,3],
                                              N0s=annual.invasion.w20[t,4], 
                                              sp, alphaSSwm, alphaSPwm20, alphaSAwm, survivalSwm)
}

#2021
#ambient
annual.invasion.a21 <- tibble(time, a, p, s)
annual.invasion.a21[1,2] <- 1 #annuals
annual.invasion.a21[1,3] <- eq.perennials.a21[300,2] #adults at eq
annual.invasion.a21[1,4] <- eq.perennials.a21[300,3] #seedlings at eq


for (t in 1:tmax) {
  annual.invasion.a21[t+1, 2] <- annual.invade(N0a=annual.invasion.a21[t,2],
                                               N0p=annual.invasion.a21[t,3],N0s=annual.invasion.a21[t,4],
                                               sa, ga, alphaAAam, alphaAPam21, lambdaAam)
  annual.invasion.a21[t+1, 4] <- seedling.resident(N0a=annual.invasion.a21[t,2],
                                                   N0p=annual.invasion.a21[t,3],N0s=annual.invasion.a21[t,4],
                                                   ss, gs, alphaPPam21, alphaPAam, lambdaPam21)
  annual.invasion.a21[t+1, 3] <- adult.resident(N0a=annual.invasion.a21[t,2],
                                                N0p=annual.invasion.a21[t,3],
                                                N0s=annual.invasion.a21[t,4], 
                                                sp, alphaSSam, alphaSPam21, alphaSAam, survivalSam)
}

#warmed
annual.invasion.w21 <- tibble(time, a, p, s)
annual.invasion.w21[1,2] <- 1
annual.invasion.w21[1,3] <- eq.perennials.w21[300,2] #adults at eq
annual.invasion.w21[1,4] <- eq.perennials.w21[300,3] #seedlings at eq

for (t in 1:tmax) {
  annual.invasion.w21[t+1, 2] <- annual.invade(N0a=annual.invasion.w21[t,2],
                                               N0p=annual.invasion.w21[t,3],N0s=annual.invasion.w21[t,4],
                                               sa, ga, alphaAAwm, alphaAPwm21, lambdaAwm)
  annual.invasion.w21[t+1, 4] <- seedling.resident(N0a=annual.invasion.w21[t,2],
                                                   N0p=annual.invasion.w21[t,3],N0s=annual.invasion.w21[t,4],
                                                   ss, gs, alphaPPwm21, alphaPAwm, lambdaPwm21)
  annual.invasion.w21[t+1, 3] <- adult.resident(N0a=annual.invasion.w21[t,2],
                                                N0p=annual.invasion.w21[t,3],
                                                N0s=annual.invasion.w21[t,4], 
                                                sp, alphaSSwm, alphaSPwm21, alphaSAwm, survivalSwm)
}


# perennials into annuals

#ambient
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))

perennial.invasion.a20 <- tibble(time, a, p, s)
perennial.invasion.a20[1,2] <-  eq.annuals.a[tmax,2] #annuals at eq
perennial.invasion.a20[1,3] <- .0001 # set based on eq. fractions 
perennial.invasion.a20[1,4] <- .999


for (t in 1:tmax) {
  perennial.invasion.a20[t+1, 2] <- annual.invade(N0a=perennial.invasion.a20[t,2],
                                                N0p=perennial.invasion.a20[t,3],N0s=perennial.invasion.a20[t,4],
                                                sa, ga, alphaAAam, alphaAPam20, lambdaAam)
  perennial.invasion.a20[t+1, 4] <- seedling.resident(N0a=perennial.invasion.a20[t,2],
                                                    N0p=perennial.invasion.a20[t,3],N0s=perennial.invasion.a20[t,4],
                                                    ss, gs, alphaPPam20, alphaPAam, lambdaPam20)
  perennial.invasion.a20[t+1, 3] <- adult.resident(N0a=perennial.invasion.a20[t,2],
                                                N0p=perennial.invasion.a20[t,3],
                                                N0s=perennial.invasion.a20[t,4], 
                                                sp, alphaSSam, alphaSPam20, alphaSAam, survivalSam)
}

#warmed
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))
perennial.invasion.w20 <- tibble(time, a, p, s)
perennial.invasion.w20[1,2] <-  eq.annuals.w[tmax,2] #annuals at eq
perennial.invasion.w20[1,3] <- .0001
perennial.invasion.w20[1,4] <- .999


for (t in 1:tmax) {
  perennial.invasion.w20[t+1, 2] <- annual.invade(N0a=perennial.invasion.w20[t,2],
                                                N0p=perennial.invasion.w20[t,3],N0s=perennial.invasion.w20[t,4],
                                                sa, ga, alphaAAwm, alphaAPwm20, lambdaAwm)
  perennial.invasion.w20[t+1, 4] <- seedling.resident(N0a=perennial.invasion.w20[t,2],
                                                    N0p=perennial.invasion.w20[t,3],N0s=perennial.invasion.w20[t,4],
                                                    ss, gs, alphaPPwm20, alphaPAwm, lambdaPwm20)
  perennial.invasion.w20[t+1, 3] <- adult.resident(N0a=perennial.invasion.w20[t,2],
                                                 N0p=perennial.invasion.w20[t,3],
                                                 N0s=perennial.invasion.w20[t,4], 
                                                 sp, alphaSSwm, alphaSPwm20, alphaSAwm, survivalSwm)
}
  
#ambient
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))

perennial.invasion.a21 <- tibble(time, a, p, s)
perennial.invasion.a21[1,2] <-  eq.annuals.a[tmax,2] #annuals at eq
perennial.invasion.a21[1,3] <- .0001 # set based on eq. fractions 
perennial.invasion.a21[1,4] <- .999


for (t in 1:tmax) {
  perennial.invasion.a21[t+1, 2] <- annual.invade(N0a=perennial.invasion.a21[t,2],
                                                  N0p=perennial.invasion.a21[t,3],N0s=perennial.invasion.a21[t,4],
                                                  sa, ga, alphaAAam, alphaAPam21, lambdaAam)
  perennial.invasion.a21[t+1, 4] <- seedling.resident(N0a=perennial.invasion.a21[t,2],
                                                      N0p=perennial.invasion.a21[t,3],N0s=perennial.invasion.a21[t,4],
                                                      ss, gs, alphaPPam21, alphaPAam, lambdaPam21)
  perennial.invasion.a21[t+1, 3] <- adult.resident(N0a=perennial.invasion.a21[t,2],
                                                   N0p=perennial.invasion.a21[t,3],
                                                   N0s=perennial.invasion.a21[t,4], 
                                                   sp, alphaSSam, alphaSPam21, alphaSAam, survivalSam)
}
#2021
#warmed
time=1:300
a <- rep(0, length(time))
s <- rep(0, length(time))
p <- rep(0, length(time))
perennial.invasion.w21 <- tibble(time, a, p, s)
perennial.invasion.w21[1,2] <-  eq.annuals.w[tmax,2] #annuals at eq
perennial.invasion.w21[1,3] <- .0001
perennial.invasion.w21[1,4] <- .999


for (t in 1:tmax) {
  perennial.invasion.w21[t+1, 2] <- annual.invade(N0a=perennial.invasion.w21[t,2],
                                                  N0p=perennial.invasion.w21[t,3],N0s=perennial.invasion.w21[t,4],
                                                  sa, ga, alphaAAwm, alphaAPwm21, lambdaAwm)
  perennial.invasion.w21[t+1, 4] <- seedling.resident(N0a=perennial.invasion.w21[t,2],
                                                      N0p=perennial.invasion.w21[t,3],N0s=perennial.invasion.w21[t,4],
                                                      ss, gs, alphaPPwm21, alphaPAwm, lambdaPwm21)
  perennial.invasion.w21[t+1, 3] <- adult.resident(N0a=perennial.invasion.w21[t,2],
                                                   N0p=perennial.invasion.w21[t,3],
                                                   N0s=perennial.invasion.w21[t,4], 
                                                   sp, alphaSSwm, alphaSPwm21, alphaSAwm, survivalSwm)
}



#visualize

# check outputs
annual.inv20<-rbind(mutate(annual.invasion.a20, temp="amb"), mutate(annual.invasion.w20, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))
annual.inv21<-rbind(mutate(annual.invasion.a21, temp="amb"), mutate(annual.invasion.w21, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))

perennial.inv20<-rbind(mutate(perennial.invasion.a20, temp="amb"), mutate(perennial.invasion.w20, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))

perennial.inv21<-rbind(mutate(perennial.invasion.a21, temp="amb"), mutate(perennial.invasion.w21, temp="warm"))%>%
  gather("type", "density", a, p, s)%>%
  mutate(type=ifelse(type=="p", "adult", ifelse(type=="s", "seedling", "annual")))

inv120<-ggplot(subset(annual.inv20), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
  geom_line(size=1.25)+ylab("density (2020)")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))
inv121<-ggplot(subset(annual.inv21), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
  geom_line(size=1.25)+ylab("density (2021)")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))

inv220<-ggplot(subset(perennial.inv20), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
geom_line(size=1.25)+ylab("density (2020)")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))
inv221<-ggplot(subset(perennial.inv21), aes(time, density,  color=type, shape=type, linetype=type))+   scale_colour_manual(values = c("darkgreen", "limegreen", "darkgreen"))+
  geom_line(size=1.25)+ylab("density (2021)")+scale_y_continuous(trans='log10', labels = scales::comma)+facet_wrap(~temp)+ scale_linetype_manual(values=c(1, 1, 3))+
  theme(text=element_text(size=16))

ggarrange(inv120, inv220, inv121, inv221, nrow=2, ncol=2, common.legend = T)


#GRWR----
  # pars contains the following parameters (subscript 1 is the annual, 2 and 3 are the perennial seedling and adult, respectively):
  # g1, g2 = germination fraction
  # lambda1, lambda3 = seed production in the absence of competition
  # alpha11, alpha13, etc. = competition coefficients
  # phi = factor scaling competitive effects of annuals to competitive effects of perennial seedlings
  # s2, s3 = over-summer survival


annualGRWRa20 <- log(annual.invasion.a20[2,2])%>%
  mutate(invader="Invader: Lolium", trt="ambient", year=2020)
annualGRWRw20 <- log(annual.invasion.w20[2,2])%>%
  mutate(invader="Invader: Lolium", trt="warmed", year=2020)
annualGRWRa21 <- log(annual.invasion.a21[2,2])%>%
  mutate(invader="Invader: Lolium", trt="ambient", year=2021)
annualGRWRw21 <- log(annual.invasion.w21[2,2])%>%
  mutate(invader="Invader: Lolium", trt="warmed", year=2021)


perennialGRWRa20 = log(.5*(sp + (4* lambdaPam20/(1 + alphaPAam20*perennial.invasion.a20[1,2]) * survivalSam20/(1 + alphaSAam20*perennial.invasion.a20[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="ambient", year=2020)
perennialGRWRw20 = log(.5*(sp + (4* lambdaPwm20/(1 + alphaPAwm20*perennial.invasion.w20[1,2]) * survivalSwm20/(1 + alphaSAwm20*perennial.invasion.w20[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="warmed", year=2020)
perennialGRWRa21 = log(.5*(sp + (4* lambdaPam21/(1 + alphaPAam21*perennial.invasion.a21[1,2]) * survivalSam21/(1 + alphaSAam21*perennial.invasion.a21[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="ambient", year=2021)
perennialGRWRw21 = log(.5*(sp + (4* lambdaPwm21/(1 + alphaPAwm21*perennial.invasion.w21[1,2]) * survivalSwm21/(1 + alphaSAwm21*perennial.invasion.w21[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="warmed", year=2021)


allGRWR<-rbind(annualGRWRa20, annualGRWRw20, perennialGRWRa20, perennialGRWRw20, annualGRWRa21, annualGRWRw21, perennialGRWRa21, perennialGRWRw21)%>%
  mutate(GRWR=a)%>%
  select(-a)



#plot of mean paremter GRWR
ggplot(allGRWR, aes(x=trt, y=GRWR)) + geom_hline(yintercept=0)+geom_bar(stat="identity", aes(fill=trt)) +facet_wrap(year~invader, scales="free") +xlab("") +ylab("Invader Growth Rate (r) When Rare") + scale_fill_manual(values = c("dodgerblue", "darkred"))

#sensitivity .. doing this manually. Change each parameter down 5%, run everything then save it as a new object
baselineGRWR <- allGRWR%>%mutate(param="baseline") #done
lambdaaGRWR <- allGRWR%>%mutate(param="lambdaa") #done
lambdapGRWR <- allGRWR%>%mutate(param="lambdap") # done
lambdasGRWR <- allGRWR%>%mutate(param="lambdas") # done


alphaapGRWR <- allGRWR%>%mutate(param="aap") # done
alphaaaGRWR <- allGRWR%>%mutate(param="aaa") #done
alphappGRWR <- allGRWR%>%mutate(param="app") # done
alphapaGRWR <- allGRWR%>%mutate(param="apa") # done
alphassGRWR <- allGRWR%>%mutate(param="ass") # done
alphasaGRWR <- allGRWR%>%mutate(param="asa") # done
alphaspGRWR <- allGRWR%>%mutate(param="asp") # done

#gaGRWR <- allGRWR%>%mutate(param="ga") #
#gsGRWR <- allGRWR%>%mutate(param="gs") #
spGRWR <- allGRWR%>%mutate(param="sp") #
#ssGRWR <- allGRWR%>%mutate(param="ss") #
#saGRWR <- allGRWR%>%mutate(param="sa") #

sensitivity<-rbind(baselineGRWR, lambdaaGRWR,lambdapGRWR,lambdasGRWR, spGRWR,# ssGRWR, saGRWR, gaGRWR, gsGRWR,
                   alphaapGRWR, alphaaaGRWR, alphappGRWR, alphapaGRWR, alphassGRWR, alphasaGRWR, alphaspGRWR  )%>%
  spread(param, GRWR)%>%
  mutate(lambdaa=(lambdaa-baseline), lambdap=(lambdap-baseline), lambdas=(lambdas-baseline), 
         sp=(sp-baseline), # ga=(ga-baseline), gs=(gs-baseline),
       #  ss=(ss-baseline), sa=(sa-baseline), 
         aap=(aap-baseline), aaa=(aaa-baseline), app=(app-baseline), 
         apa=(apa-baseline), ass=(ass-baseline), asa=(asa-baseline), 
         asp=(asp-baseline) )%>%
  gather(param, GRWRdiff, 4:15)%>%
  filter(param!="baseline")%>%
  group_by(param, invader, trt)%>%
  mutate(GRWRdiff=mean(GRWRdiff))

sensitivity$param<-factor(sensitivity$param, levels = c("lambdap", "lambdas", "lambdaa", "sp", "app", "apa", "asp",  "ass", "asa", "aap", "aaa"))

ggplot(subset(sensitivity), aes(param, GRWRdiff))+
  geom_bar(stat="identity", aes(fill=trt), position="dodge")+
  facet_wrap(~invader, scales="free")+
  geom_hline(yintercept = 0)+ theme(axis.text.x=element_text(angle = -90, hjust = 0))+
scale_fill_manual(values = c("dodgerblue", "darkred"))

#final populations ----
annual.annual<-rbind((annual.invasion.a20[300,2])%>%
                      mutate(type="lolium", invader="Invader: Lolium", trt="ambient", year=2020),
                    (annual.invasion.w20[300,2])%>%
                      mutate(type="lolium", invader="Invader: Lolium", trt="warmed", year=2020),
                    (annual.invasion.a21[300,2])%>%
                      mutate(type="lolium", invader="Invader: Lolium", trt="ambient", year=2021),
                    (annual.invasion.w21[300,2])%>%
                      mutate(type="lolium", invader="Invader: Lolium", trt="warmed", year=2021))%>%
  mutate(eqpop=a)%>%
  select(-a)
                    
annual.seedling<-rbind((annual.invasion.a20[300,4])%>%
                      mutate(type="festuca seedling", invader="Invader: Lolium", trt="ambient", year=2020),
                    (annual.invasion.w20[300,4])%>%
                      mutate(type="festuca seedling", invader="Invader: Lolium", trt="warmed", year=2020),
                    (annual.invasion.a21[300,4])%>%
                      mutate(type="festuca seedling", invader="Invader: Lolium", trt="ambient", year=2021),
                    (annual.invasion.w21[300,4])%>%
                      mutate(type="festuca seedling", invader="Invader: Lolium", trt="warmed", year=2021))%>%
  mutate(eqpop=s)%>%
  select(-s)

annual.adult<-rbind((annual.invasion.a20[300,3])%>%
                       mutate(type="festuca adult", invader="Invader: Lolium", trt="ambient", year=2020),
                     (annual.invasion.w20[300,3])%>%
                       mutate(type="festuca adult", invader="Invader: Lolium", trt="warmed", year=2020),
                     (annual.invasion.a21[300,3])%>%
                       mutate(type="festuca adult", invader="Invader: Lolium", trt="ambient", year=2021),
                     (annual.invasion.w21[300,3])%>%
                       mutate(type="festuca adult", invader="Invader: Lolium", trt="warmed", year=2021))%>%
  mutate(eqpop=p)%>%
  select(-p)

annual_eqpop<-rbind(annual.annual, annual.adult, annual.seedling)





#plot of equilibrium plots when annuals are invadiang

nocomp<-rbind(
  mutate(annual.eq, type="lolium", density=a)%>%select(-a),
  mutate(perennial.eq, type=ifelse(type=="adult", "festuca adult", "festuca seedling")
))%>%
  filter(time==300)%>%
  mutate(density=ifelse(temp=="amb"&type=="lolium", 22586, density))%>%
  mutate(density=ifelse(temp=="warm"&type=="lolium", 14581, density))%>%
  filter(year!=0)
  
  

tog<-mutate(nocomp, year=ifelse(year==20, 2020, 2021))%>%
  mutate(trt=ifelse(temp=="amb", "ambient", "warmed"))
tog<-left_join(tog, annual_eqpop)%>%
  mutate(prop=eqpop/density)

r1<- ggplot(nocomp, aes(x=as.factor(year), y=density)) + 
  geom_bar(stat="identity", position="dodge", aes(fill=temp)) +
  facet_wrap(~type, scales="free") +xlab("") +
  ylab("Equilibrium population without competition") + 
  scale_fill_manual(values = c("dodgerblue", "darkred"))

r2 <- ggplot(annual_eqpop, aes(x=as.factor(year), y=eqpop)) + 
  geom_bar(stat="identity", position="dodge", aes(fill=trt)) +
  facet_wrap(~type, scales="free") +xlab("") +
  ylab("Equilibrium population with competition") + 
  scale_fill_manual(values = c("dodgerblue", "darkred"))

r3 <- ggplot(tog, aes(x=as.factor(year), y=prop)) + 
  geom_bar(stat="identity", position="dodge", aes(fill=trt)) +
  facet_wrap(~type, scales="free") +xlab("") +
  ylab("Proportion of equilibrium population with intra-specific competition ") + 
  scale_fill_manual(values = c("dodgerblue", "darkred"))


ggarrange(r1, r2, r3, nrow=3, ncol=1)

#plot of equilibrium plots when perennials are invadiang
#ggplot(perennial_eqpop,aes(x=as.factor(year), y=eqpop)) + 
#  geom_bar(stat="identity", position="dodge", aes(fill=trt)) +
#  facet_wrap(~type, scales="free") +xlab("") +
#  ylab("Equilibrium population under perennial invasion") + 
#  scale_fill_manual(values = c("dodgerblue", "darkred"))
s
perennialGRWRa20 = log(.5*(sp + (4* lambdaPam20/(1 + alphaPAam20*perennial.invasion.a20[1,2]) * survivalSam20/(1 + alphaSAam20*perennial.invasion.a20[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="ambient", year=2020)
perennialGRWRw20 = log(.5*(sp + (4* lambdaPwm20/(1 + alphaPAwm20*perennial.invasion.w20[1,2]) * survivalSwm20/(1 + alphaSAwm20*perennial.invasion.w20[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="warmed", year=2020)
perennialGRWRa21 = log(.5*(sp + (4* lambdaPam21/(1 + alphaPAam21*perennial.invasion.a21[1,2]) * survivalSam21/(1 + alphaSAam21*perennial.invasion.a21[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="ambient", year=2021)
perennialGRWRw21 = log(.5*(sp + (4* lambdaPwm21/(1 + alphaPAwm21*perennial.invasion.w21[1,2]) * survivalSwm21/(1 + alphaSAwm21*perennial.invasion.w21[1,2]) + sp^2)^0.5))%>%
  mutate(invader="Invader: Festuca", trt="warmed", year=2021)


allGRWR<-rbind(annualGRWRa20, annualGRWRw20, perennialGRWRa20, perennialGRWRw20, annualGRWRa21, annualGRWRw21, perennialGRWRa21, perennialGRWRw21)%>%
  mutate(GRWR=a)%>%
  select(-a)


#GRWR with mean eq but all the parameters for invasion from there

#GRWR with all the chains ----
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


#all chains (posterior predictive plots) ----

#2020
#annuals
iterations_a20<-(as.data.frame(rstan::extract(annuals_fit2020)))%>%
  mutate(lambdaA_warm=lambdaA_amb+lambdaA_slope, 
         alphaAA_warm=alphaAA_amb+alphaAA_slope,
         alphaAP_warm=alphaAP_amb+alphaAP_slope)%>%
  select(-lambdaA_slope, -alphaAA_slope, -alphaAP_slope)
iterations_a20<-rowid_to_column(iterations_a20)%>%
  gather(param, value, 2:8)%>%
  mutate(exp=exp(value))%>%
  select(-value)%>%
  spread(param, exp)

#adults
iterations_p20<-(as.data.frame(rstan::extract(adults_fit2020)))%>%
  mutate(lambdaP_warm=lambdaP_amb+lambdaP_slope, 
         alphaPA_warm=alphaPA_amb+alphaPA_slope,
         alphaPP_warm=alphaPP_amb+alphaPP_slope)%>%
  select(-lambdaP_slope, -alphaPA_slope, -alphaPP_slope)
iterations_p20<-rowid_to_column(iterations_p20)%>%
  gather(param, value, 2:8)%>%
  mutate(exp=exp(value))%>%
  select(-value)%>%
  spread(param, exp)

#seedlings
iterations_s20<-(as.data.frame(rstan::extract(seedlings_fit2020)))%>%
  mutate(survivalS_warm=survivalS_amb+survivalS_slope, 
         alphaSA_warm=alphaSA_amb+alphaSA_slope,
         alphaSP_warm=alphaSP_amb+alphaSP_slope,
         alphaSS_warm=alphaSS_amb+alphaSS_slope)%>%
  select(-survivalS_slope, -alphaSA_slope, -alphaSP_slope, -alphaSS_slope)
iterations_s20<-rowid_to_column(iterations_s20)%>%
  gather(param, value, 2:9)%>%
  mutate(exp=exp(value))%>%
  select(-value)%>%
  spread(param, exp)

#2021

#annuals
iterations_a21<-(as.data.frame(rstan::extract(annuals_fit2021)))%>%
  mutate(lambdaA_warm=lambdaA_amb+lambdaA_slope, 
         alphaAA_warm=alphaAA_amb+alphaAA_slope,
         alphaAP_warm=alphaAP_amb+alphaAP_slope)%>%
  select(-lambdaA_slope, -alphaAA_slope, -alphaAP_slope)
iterations_a21<-rowid_to_column(iterations_a21)%>%
  gather(param, value, 2:8)%>%
  mutate(exp=exp(value))%>%
  select(-value)%>%
  spread(param, exp)

#adults
iterations_p21<-(as.data.frame(rstan::extract(adults_fit2021)))%>%
  mutate(lambdaP_warm=lambdaP_amb+lambdaP_slope, 
         alphaPA_warm=alphaPA_amb+alphaPA_slope,
         alphaPP_warm=alphaPP_amb+alphaPP_slope)%>%
  select(-lambdaP_slope, -alphaPA_slope, -alphaPP_slope)
iterations_p21<-rowid_to_column(iterations_p21)%>%
  gather(param, value, 2:8)%>%
  mutate(exp=exp(value))%>%
  select(-value)%>%
  spread(param, exp)

#seedlings
iterations_s21<-(as.data.frame(rstan::extract(seedlings_fit2021)))%>%
  mutate(survivalS_warm=survivalS_amb+survivalS_slope, 
         alphaSA_warm=alphaSA_amb+alphaSA_slope,
         alphaSP_warm=alphaSP_amb+alphaSP_slope,
         alphaSS_warm=alphaSS_amb+alphaSS_slope)%>%
  select(-survivalS_slope, -alphaSA_slope, -alphaSP_slope, -alphaSS_slope)
iterations_s21<-rowid_to_column(iterations_s21)%>%
  gather(param, value, 2:9)%>%
  mutate(exp=exp(value))%>%
  select(-value)%>%
  spread(param, exp)


year1<-cbind(subset(iterations_a20, rowid<1001), subset(iterations_p20, rowid<1001), subset(iterations_s20, rowid<1001))%>%
  select(-8, -9, -16, -17, -24)%>%
  gather("param", "value", 2:21)%>%
  mutate(warmtrt=ifelse(grepl("warm", param), "warmed", "ambient"))%>%
  mutate(stage=ifelse(grepl("AA", param), "annual", NA))%>%
  mutate(stage=ifelse(grepl("AP", param), "annual", stage))%>%
  mutate(stage=ifelse(grepl("lambdaA", param), "annual", stage))%>%
  mutate(stage=ifelse(grepl("PA", param), "adult", stage))%>%
  mutate(stage=ifelse(grepl("PP", param), "adult", stage))%>%
  mutate(stage=ifelse(grepl("lambdaP", param), "adult", stage))%>%
  mutate(stage=ifelse(grepl("SA", param), "seedling", stage))%>%
  mutate(stage=ifelse(grepl("SP", param), "seedling", stage))%>%
  mutate(stage=ifelse(grepl("survivalS", param), "seedling", stage))%>%
  mutate(stage=ifelse(grepl("SS", param), "seedling", stage))%>%
  mutate(type=ifelse(grepl("alpha", param), "alpha", "lambda"))%>%
  mutate(type=ifelse(grepl("AA", param), "self_alpha", type))%>%
  mutate(type=ifelse(grepl("SS", param), "self_alpha", type))%>%
  mutate(type=ifelse(grepl("PP", param), "self_alpha", type))%>%
  mutate(type=ifelse(grepl("AP", param), "adult_alpha", type))%>%
  mutate(type=ifelse(grepl("SP", param), "adult_alpha", type))%>%
  mutate(type=ifelse(grepl("SA", param), "annual_alpha", type))%>%
  mutate(type=ifelse(grepl("PA", param), "annual_alpha", type))%>%
  mutate(year=2020)

year2<-cbind(subset(iterations_a21, rowid<1001), subset(iterations_p21, rowid<1001), subset(iterations_s21, rowid<1001))%>%
  select(-8, -9, -16, -17, -24)%>%
  gather("param", "value", 2:21)%>%
  mutate(warmtrt=ifelse(grepl("warm", param), "warmed", "ambient"))%>%
  mutate(stage=ifelse(grepl("AA", param), "annual", NA))%>%
  mutate(stage=ifelse(grepl("AP", param), "annual", stage))%>%
  mutate(stage=ifelse(grepl("lambdaA", param), "annual", stage))%>%
  mutate(stage=ifelse(grepl("PA", param), "adult", stage))%>%
  mutate(stage=ifelse(grepl("PP", param), "adult", stage))%>%
  mutate(stage=ifelse(grepl("lambdaP", param), "adult", stage))%>%
  mutate(stage=ifelse(grepl("SA", param), "seedling", stage))%>%
  mutate(stage=ifelse(grepl("SP", param), "seedling", stage))%>%
  mutate(stage=ifelse(grepl("survivalS", param), "seedling", stage))%>%
  mutate(stage=ifelse(grepl("SS", param), "seedling", stage))%>%
  mutate(type=ifelse(grepl("alpha", param), "alpha", "lambda"))%>%
  mutate(type=ifelse(grepl("AA", param), "self_alpha", type))%>%
  mutate(type=ifelse(grepl("SS", param), "self_alpha", type))%>%
  mutate(type=ifelse(grepl("PP", param), "self_alpha", type))%>%
  mutate(type=ifelse(grepl("AP", param), "adult_alpha", type))%>%
  mutate(type=ifelse(grepl("SP", param), "adult_alpha", type))%>%
  mutate(type=ifelse(grepl("SA", param), "annual_alpha", type))%>%
  mutate(type=ifelse(grepl("PA", param), "annual_alpha", type))%>%
  mutate(year=2021)

alliter<-rbind(year1, year2)


ggplot(subset(alliter, stage=="annual"), aes(as.factor(year), value, color=warmtrt)) +geom_boxplot()  +facet_wrap(~type, scale='free') +ylab("parameter value") +scale_color_manual(values=c("blue", "red"))
ggplot(subset(alliter, stage=="adult"), aes(as.factor(year), value, color=warmtrt)) +geom_boxplot()  +facet_wrap(~type, scale='free') +ylab("parameter value") +scale_color_manual(values=c("blue", "red"))
ggplot(subset(alliter, stage=="seedling"), aes(as.factor(year), value, color=warmtrt)) +geom_boxplot()  +facet_wrap(~type, scale='free') +ylab("parameter value") +scale_color_manual(values=c("blue", "red"))


lambdaAa20<-  iterations_a20$lambdaA_amb
lambdaAw20<- iterations_a20$lambdaA_warm
survivalSa20<- iterations_s20$survivalS_amb
survivalSw20<- iterations_s20$survivalS_warm
lambdaPa20<- iterations_p20$lambdaP_amb
lambdaPw20<- iterations_p20$lambdaP_warm

alphaAAa20<- iterations_a20$alphaAA_amb
alphaAAw20<- iterations_a20$alphaAA_warm
alphaAPa20<- iterations_a20$alphaAP_amb
alphaAPw20<- iterations_a20$alphaAP_warm

alphaPAa20<- iterations_p20$alphaPA_amb
alphaPAw20<-  iterations_p20$alphaPA_warm
alphaPPa20<- iterations_p20$alphaPP_amb
alphaPPw20<- iterations_p20$alphaPP_warm

alphaSAa20<- iterations_s20$alphaSA_amb
alphaSAw20<- iterations_s20$alphaSA_warm
alphaSSa20<- iterations_s20$alphaSS_amb
alphaSSw20<- iterations_s20$alphaSS_warm
alphaSPa20<- iterations_s20$alphaSP_amb
alphaSPw20<- iterations_s20$alphaSP_warm



lambdaAa21<-  iterations_a21$lambdaA_amb
lambdaAw21<- iterations_a21$lambdaA_warm
survivalSa21<- iterations_s21$survivalS_amb
survivalSw21<- iterations_s21$survivalS_warm
lambdaPa21<- iterations_p21$lambdaP_amb
lambdaPw21<- iterations_p21$lambdaP_warm

alphaAAa21<- iterations_a21$alphaAA_amb
alphaAAw21<- iterations_a21$alphaAA_warm
alphaAPa21<- iterations_a21$alphaAP_amb
alphaAPw21<- iterations_a21$alphaAP_warm

alphaPAa21<- iterations_p21$alphaPA_amb
alphaPAw21<-  iterations_p21$alphaPA_warm
alphaPPa21<- iterations_p21$alphaPP_amb
alphaPPw21<- iterations_p21$alphaPP_warm

alphaSAa21<- iterations_s21$alphaSA_amb
alphaSAw21<- iterations_s21$alphaSA_warm
alphaSSa21<- iterations_s21$alphaSS_amb
alphaSSw21<- iterations_s21$alphaSS_warm
alphaSPa21<- iterations_s21$alphaSP_amb
alphaSPw21<- iterations_s21$alphaSP_warm


# run equilibriums add 100 chains for each to visualization
#warmed.annual.chains<-as.tibble(lambdaAw20)%>%mutate(eq=0)
#ambient.annual.chains<-as.tibble(lambdaAa20)%>%mutate(eq=0)

#warmed.seedling.chains<-as.tibble(survivalSw20)%>%mutate(eq=0)
#ambient.seedling.chains<-as.tibble(survivalSa20)%>%mutate(eq=0)

#warmed.adult.chains<-as.tibble(lambdaPw20)%>%mutate(eq=0)
#ambient.adult.chains<-as.tibble(lambdaPa20)%>%mutate(eq=0)

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

#2020
warmed.annual.chains120<-as.tibble(lambdaAw20[1:100])
warmed.annual.chains120 <- cbind(warmed.annual.chains120, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0a)

ambient.annual.chains120<-as.tibble(lambdaAa20[1:100])
ambient.annual.chains120 <- cbind(ambient.annual.chains120, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0a)


for (i in 1:100) {
  for (t in 1:tmax) {
    ambient.annual.chains120[i, t+2] <- annual.equilibrium(ambient.annual.chains120[i, t+1], sa, ga, alphaAAa20[i], lambdaAa20[i])
  }
  for (t in 1:tmax) {
    warmed.annual.chains120[i, t+2] <- annual.equilibrium(warmed.annual.chains120[i, t+1], sa, ga, alphaAAw20[i], lambdaAw20[i])
  }
}

#2021
warmed.annual.chains121<-as.tibble(lambdaAw21[1:100])
warmed.annual.chains121 <- cbind(warmed.annual.chains121, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0a)

ambient.annual.chains121<-as.tibble(lambdaAa21[1:100])
ambient.annual.chains121 <- cbind(ambient.annual.chains121, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0a)


for (i in 1:100) {
  for (t in 1:tmax) {
    ambient.annual.chains121[i, t+2] <- annual.equilibrium(ambient.annual.chains121[i, t+1], sa, ga, alphaAAa21[i], lambdaAa21[i])
  }
  for (t in 1:tmax) {
    warmed.annual.chains121[i, t+2] <- annual.equilibrium(warmed.annual.chains121[i, t+1], sa, ga, alphaAAw21[i], lambdaAw21[i])
  }
}


#perennials

tmax=300

#set up empty dataframes with rows as iterations and columns as time
warmed.seedling.chains120<-as.tibble(survivalSw20[1:100])
warmed.seedling.chains120 <- cbind(warmed.seedling.chains120, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0s)

ambient.seedling.chains120<-as.tibble(survivalSa20[1:100])
ambient.seedling.chains120 <- cbind(ambient.seedling.chains120, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0s)

warmed.adult.chains120<-as.tibble(lambdaPw20[1:100])
warmed.adult.chains120 <- cbind(warmed.adult.chains120, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0p)

ambient.adult.chains120<-as.tibble(lambdaPa20[1:100])
ambient.adult.chains120 <- cbind(ambient.adult.chains120, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0p)

#populate data in the data frames
for (i in 1:100) {  # this value is always 2000, no need to change it (just to get this to run I do)
  for (t in 1:tmax) {
    ambient.adult.chains120[i, t+2] <- adult.equilibrium(ambient.adult.chains120[i, t+1], ambient.seedling.chains120[i, t+1], sp, alphaSSa20[i], alphaSPa20[i], survivalSa20[i])
    ambient.seedling.chains120[i, t+2] <- seedling.equilibrium(ambient.adult.chains120[i, t+1], ambient.seedling.chains120[i, t+1], ss, gs, alphaPPa20[i], lambdaPa20[i])
  }
  for (t in 1:tmax) {
    warmed.adult.chains120[i, t+2] <- adult.equilibrium(warmed.adult.chains120[i, t+1], warmed.seedling.chains120[i, t+1], sp, alphaSSw20[i], alphaSPw20[i], survivalSw20[i])
    warmed.seedling.chains120[i, t+2] <- seedling.equilibrium(warmed.adult.chains120[i,  t+1], warmed.seedling.chains120[i, t+1], ss, gs, alphaPPw20[i], lambdaPw20[i])
  }
}


#2021
warmed.seedling.chains121<-as.tibble(survivalSw21[1:100])
warmed.seedling.chains121 <- cbind(warmed.seedling.chains121, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0s)

ambient.seedling.chains121<-as.tibble(survivalSa21[1:100])
ambient.seedling.chains121 <- cbind(ambient.seedling.chains121, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0s)

warmed.adult.chains121<-as.tibble(lambdaPw21[1:100])
warmed.adult.chains121 <- cbind(warmed.adult.chains121, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0p)

ambient.adult.chains121<-as.tibble(lambdaPa21[1:100])
ambient.adult.chains121 <- cbind(ambient.adult.chains121, t(rep(0,tmax+1)))%>%
  mutate(`1`=N0p)

#populate data in the daata frames
for (i in 1:100) {  # this value is always 2100, no need to change it (just to get this to run I do)
  for (t in 1:tmax) {
    ambient.adult.chains121[i, t+2] <- adult.equilibrium(ambient.adult.chains121[i, t+1], ambient.seedling.chains121[i, t+1], sp, alphaSSa21[i], alphaSPa21[i], survivalSa21[i])
    ambient.seedling.chains121[i, t+2] <- seedling.equilibrium(ambient.adult.chains121[i, t+1], ambient.seedling.chains121[i, t+1], ss, gs, alphaPPa21[i], lambdaPa21[i])
  }
  for (t in 1:tmax) {
    warmed.adult.chains121[i, t+2] <- adult.equilibrium(warmed.adult.chains121[i, t+1], warmed.seedling.chains121[i, t+1], sp, alphaSSw21[i], alphaSPw21[i], survivalSw21[i])
    warmed.seedling.chains121[i, t+2] <- seedling.equilibrium(warmed.adult.chains121[i, t+1], warmed.seedling.chains121[i, t+1], ss, gs, alphaPPw21[i], lambdaPw21[i])
  }
}


ambient.annual.chains120<-rowid_to_column(ambient.annual.chains120)
warmed.annual.chains120<-rowid_to_column(warmed.annual.chains120)
ambient.annual.chains121<-rowid_to_column(ambient.annual.chains121)
warmed.annual.chains121<-rowid_to_column(warmed.annual.chains121)

ambient.seedling.chains120<-rowid_to_column(ambient.seedling.chains120)
warmed.seedling.chains120<-rowid_to_column(warmed.seedling.chains120)
ambient.seedling.chains121<-rowid_to_column(ambient.seedling.chains121)
warmed.seedling.chains121<-rowid_to_column(warmed.seedling.chains121)

ambient.adult.chains120<-rowid_to_column(ambient.adult.chains120)
warmed.adult.chains120<-rowid_to_column(warmed.adult.chains120)
ambient.adult.chains121<-rowid_to_column(ambient.adult.chains121)
warmed.adult.chains121<-rowid_to_column(warmed.adult.chains121)

# visualize chains
annual.chains20<-  rbind(mutate(ambient.annual.chains120, treatment="amb")[1:100,], mutate(warmed.annual.chains120, treatment="warm")[1:100,])%>%
  gather("time", "population", 3:13)%>% #this only runs to 14 because it only runs for 10 years
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(as.character(rowid), as.character(treatment), sep="_"))

seedling.chains20<-rbind(mutate(ambient.seedling.chains120, treatment="amb")[1:50,],  mutate(warmed.seedling.chains120, treatment="warm")[1:50,])%>%
  gather("time", "population", 3:303) %>% # this runs to 304 because it takes 300 years to equilibrate
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(as.character(rowid), as.character(treatment), sep="_"))

adult.chains20<-   rbind(mutate(ambient.adult.chains120, treatment="amb")[1:50,], mutate(warmed.adult.chains120, treatment="warm")[1:50,])%>%
  gather("time", "population", 3:303)%>%
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(as.character(rowid), as.character(treatment), sep="_"))


annual.chains21<-  rbind(mutate(ambient.annual.chains121, treatment="amb")[1:100,], mutate(warmed.annual.chains121, treatment="warm")[1:100,])%>%
  gather("time", "population", 3:13)%>% #this only runs to 14 because it only runs for 10 years
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(as.character(rowid), as.character(treatment), sep="_"))

seedling.chains21<-rbind(mutate(ambient.seedling.chains121, treatment="amb")[1:50,],  mutate(warmed.seedling.chains121, treatment="warm")[1:50,])%>%
  gather("time", "population", 3:303) %>% # this runs to 304 because it takes 300 years to equilibrate
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(as.character(rowid), as.character(treatment), sep="_"))

adult.chains21<-   rbind(mutate(ambient.adult.chains121, treatment="amb")[1:50,], mutate(warmed.adult.chains121, treatment="warm")[1:50,])%>%
  gather("time", "population", 3:303)%>%
  mutate(time=as.numeric(time))%>%
  mutate(id=paste(as.character(rowid), as.character(treatment), sep="_"))


eqannual20<-ggplot(subset(annual.chains20, population>0), aes(time, (population)))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("annual seed equilibrium density 2020")+
  #geom_line(data=subset(annual.eq, time<12), size=1.5, aes(x=time, y=a, color=temp))+ #?
  scale_colour_manual(values = c("dodgerblue",  "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)

eqannual21<-ggplot(subset(annual.chains21, population>0), aes(time, (population)))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("annual seed equilibrium density 2021")+
  #geom_line(data=subset(annual.eq, time<12), size=1.5, aes(x=time, y=a, color=temp))+ #?
  scale_colour_manual(values = c("dodgerblue",  "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)

eqseedling20<-ggplot(subset(seedling.chains20, population>0), aes(time, population))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("perennial seed equilibrium density 2020")+
  # geom_line(data=subset(perennial.eq, type=="seedling"), size=1.5, aes(x=time, y=density, color=temp))+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)

eqseedling21<-ggplot(subset(seedling.chains21, population>0), aes(time, population))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("perennial seed equilibrium density 2021")+
  # geom_line(data=subset(perennial.eq, type=="seedling"), size=1.5, aes(x=time, y=density, color=temp))+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)


eqadult20<-ggplot(subset(adult.chains20, population>0), aes(time, (population)))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("perennial adult equilibrium density 2020")+
  geom_line(data=subset(subset(perennial.eq,year==20) , type=="adult"), size=1.5, aes(x=time, y=density, color=temp))+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)

eqadult21<-ggplot(subset(adult.chains21, population>0), aes(time, (population)))+ 
  geom_line(size=.3, alpha=.15, aes(group=id, color=treatment))+ylab("perennial adult equilibrium density 2021")+
  geom_line(data=subset(subset(perennial.eq,year==21) ,  type=="adult"), size=1.5, aes(x=time, y=density, color=temp))+
  scale_colour_manual(values = c("dodgerblue", "darkred"))+scale_y_continuous(trans='log10', labels = scales::comma)

## population ove time
ggarrange(eqannual20, eqannual21, nrow=1, ncol=2, common.legend = T)
ggarrange(eqadult20, eqadult21, nrow=1, ncol=2, common.legend = T)
ggarrange(eqseedling20, eqseedling21, nrow=1, ncol=2, common.legend = T)






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


#visualize ----

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


#sensitivity test


