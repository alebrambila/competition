## predict and plot COUNTERFACTUALS----


# annuals ----
#2020 ----
annuals2020b<-mutate(annuals2020, warmtrt=ifelse(warmtrt==0, "amb", "warm"))

### predict over pm2 (hold seeded_a steady)
annuals2020_ahi <- expand.grid(
  seeded_am2 = 1800 # near the max of seeded_am2
  ,starting_pm2 = as.numeric(seq(1,9, length.out=200))
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2020[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[1, 3]) + starting_pm2*as.numeric(annuals_estimates2020[3, 3])),
                         as.numeric(annuals_estimates2020[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[2, 3]) + starting_pm2*as.numeric(annuals_estimates2020[4, 3]))))

annuals2020_amean <- expand.grid(
  seeded_am2 = mean(annuals$seeded_am2, na.rm=TRUE) #0
  ,starting_pm2 = seq(1,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2020[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[1, 3]) + starting_pm2*as.numeric(annuals_estimates2020[3, 3])),
                         as.numeric(annuals_estimates2020[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[2, 3]) + starting_pm2*as.numeric(annuals_estimates2020[4, 3]))))


annuals2020_alo <- expand.grid(
  seeded_am2 = 75 #above the min
  ,starting_pm2 = seq(1,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2020[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[1, 3]) + starting_pm2*as.numeric(annuals_estimates2020[3, 3])),
                         as.numeric(annuals_estimates2020[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[2, 3]) + starting_pm2*as.numeric(annuals_estimates2020[4, 3]))))



a<-ggplot(data=annuals2020_alo, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2020b, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2020 w/seeded_a=75")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
b<-ggplot(data=annuals2020_amean, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2020b, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2020 w/seeded_a=685")+ scale_color_manual(values=c("blue", "red"))+
theme_classic()+theme(legend.position="none")
c<-ggplot(data=annuals2020_ahi, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2020b, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2020 w/seeded_a=1800")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")


# predict over seeded_a (hold pm2 steady)

annuals2020_phi <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 9 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2020[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[1, 3]) + starting_pm2*as.numeric(annuals_estimates2020[3, 3])),
                         as.numeric(annuals_estimates2020[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[2, 3]) + starting_pm2*as.numeric(annuals_estimates2020[4, 3]))))

annuals2020_pmean <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 5
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2020[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[1, 3]) + starting_pm2*as.numeric(annuals_estimates2020[3, 3])),
                         as.numeric(annuals_estimates2020[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[2, 3]) + starting_pm2*as.numeric(annuals_estimates2020[4, 3]))))

annuals2020_plo <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 1 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2020[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[1, 3]) + starting_pm2*as.numeric(annuals_estimates2020[3, 3])),
                         as.numeric(annuals_estimates2020[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2020[2, 3]) + starting_pm2*as.numeric(annuals_estimates2020[4, 3]))))


d<-ggplot(data=annuals2020_plo, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2020b, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2020 w/starting_p=1")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
e<-ggplot(data=annuals2020_pmean, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2020b, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2020 w/starting_p=5")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
f<-ggplot(data=annuals2020_phi, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2020b, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2020 w/starting_p=9")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()


#visualize conditional plots first with all data
ggarrange(a, b, c, d, e, f, nrow=2, ncol=3, common.legend=T)

#2021 ----
annuals2021b<-mutate(annuals2021, warmtrt=ifelse(warmtrt==0, "amb", "warm"))

### predict over pm2 (hold seeded_a steady)
annuals2021_ahi <- expand.grid(
  seeded_am2 = 1800 # near the max of seeded_am2
  ,starting_pm2 = as.numeric(seq(0,9, length.out=200))
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2021[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[1, 3]) + starting_pm2*as.numeric(annuals_estimates2021[3, 3])),
                         as.numeric(annuals_estimates2021[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[2, 3]) + starting_pm2*as.numeric(annuals_estimates2021[4, 3]))))

annuals2021_amean <- expand.grid(
  seeded_am2 = mean(annuals$seeded_am2, na.rm=TRUE) #0
  ,starting_pm2 = seq(0,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2021[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[1, 3]) + starting_pm2*as.numeric(annuals_estimates2021[3, 3])),
                         as.numeric(annuals_estimates2021[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[2, 3]) + starting_pm2*as.numeric(annuals_estimates2021[4, 3]))))


annuals2021_alo <- expand.grid(
  seeded_am2 = 75 #above the min
  ,starting_pm2 = seq(0,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2021[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[1, 3]) + starting_pm2*as.numeric(annuals_estimates2021[3, 3])),
                         as.numeric(annuals_estimates2021[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[2, 3]) + starting_pm2*as.numeric(annuals_estimates2021[4, 3]))))



a<-ggplot(data=annuals2021_alo, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2021b, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2021 w/seeded_a=75")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
b<-ggplot(data=annuals2021_amean, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2021b, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2021 w/seeded_a=685")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
c<-ggplot(data=annuals2021_ahi, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2021b, aes(x=starting_pm2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2021 w/seeded_a=1800")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")


# predict over seeded_a (hold pm2 steady)

annuals2021_phi <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 9 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2021[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[1, 3]) + starting_pm2*as.numeric(annuals_estimates2021[3, 3])),
                         as.numeric(annuals_estimates2021[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[2, 3]) + starting_pm2*as.numeric(annuals_estimates2021[4, 3]))))

annuals2021_pmean <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 5
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2021[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[1, 3]) + starting_pm2*as.numeric(annuals_estimates2021[3, 3])),
                         as.numeric(annuals_estimates2021[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[2, 3]) + starting_pm2*as.numeric(annuals_estimates2021[4, 3]))))

annuals2021_plo <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 1 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(annuals_estimates2021[5,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[1, 3]) + starting_pm2*as.numeric(annuals_estimates2021[3, 3])),
                         as.numeric(annuals_estimates2021[6,3])/(1+seeded_am2*as.numeric(annuals_estimates2021[2, 3]) + starting_pm2*as.numeric(annuals_estimates2021[4, 3]))))


d<-ggplot(data=annuals2021_plo, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2021b, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2021 w/starting_p=1")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
e<-ggplot(data=annuals2021_pmean, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2021b, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2021 w/starting_p=5")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
f<-ggplot(data=annuals2021_phi, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=annuals2021b, aes(x=seeded_am2, y=percap), shape=1, width=.25)+
  ylab("Annual fecundity 2021 w/starting_p=9")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()


#visualize conditional plots first with all data
ggarrange(a, b, c, d, e, f, nrow=2, ncol=3, common.legend=T)



# perennial adults
# adults ----
#2020 ----
adults2020b<-mutate(adults2020, warmtrt=ifelse(warmtrt==0, "amb", "warm"))

### predict over pm2 (hold seeded_a steady)
adults2020_ahi <- expand.grid(
  seeded_am2 = 1800 # near the max of seeded_am2
  ,starting_pm2 = as.numeric(seq(0,9, length.out=200))
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2020[5,3])/(1+seeded_am2*as.numeric(adults_estimates2020[1, 3]) + starting_pm2*as.numeric(adults_estimates2020[3, 3])),
                         as.numeric(adults_estimates2020[6,3])/(1+seeded_am2*as.numeric(adults_estimates2020[2, 3]) + starting_pm2*as.numeric(adults_estimates2020[4, 3]))))

adults2020_amean <- expand.grid(
  seeded_am2 = mean(annuals$seeded_am2, na.rm=TRUE) #0
  ,starting_pm2 = seq(0,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2020[5,3])/(1+seeded_am2*as.numeric(adults_estimates2020[1, 3]) + starting_pm2*as.numeric(adults_estimates2020[3, 3])),
                         as.numeric(adults_estimates2020[6,3])/(1+seeded_am2*as.numeric(adults_estimates2020[2, 3]) + starting_pm2*as.numeric(adults_estimates2020[4, 3]))))


adults2020_alo <- expand.grid(
  seeded_am2 = 75 #above the min
  ,starting_pm2 = seq(0,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2020[5,3])/(1+seeded_am2*as.numeric(adults_estimates2020[1, 3]) + starting_pm2*as.numeric(adults_estimates2020[3, 3])),
                         as.numeric(adults_estimates2020[6,3])/(1+seeded_am2*as.numeric(adults_estimates2020[2, 3]) + starting_pm2*as.numeric(adults_estimates2020[4, 3]))))



a<-ggplot(data=adults2020_alo, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2020b, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2020 w/seeded_a=75")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
b<-ggplot(data=adults2020_amean, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2020b, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2020 w/seeded_a=685")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
c<-ggplot(data=adults2020_ahi, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2020b, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2020 w/seeded_a=1800")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")


# predict over seeded_a (hold pm2 steady)

adults2020_phi <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 9 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2020[5,3])/(1+seeded_am2*as.numeric(adults_estimates2020[1, 3]) + starting_pm2*as.numeric(adults_estimates2020[3, 3])),
                         as.numeric(adults_estimates2020[6,3])/(1+seeded_am2*as.numeric(adults_estimates2020[2, 3]) + starting_pm2*as.numeric(adults_estimates2020[4, 3]))))

adults2020_pmean <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 5
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2020[5,3])/(1+seeded_am2*as.numeric(adults_estimates2020[1, 3]) + starting_pm2*as.numeric(adults_estimates2020[3, 3])),
                         as.numeric(adults_estimates2020[6,3])/(1+seeded_am2*as.numeric(adults_estimates2020[2, 3]) + starting_pm2*as.numeric(adults_estimates2020[4, 3]))))

adults2020_plo <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 1 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2020[5,3])/(1+seeded_am2*as.numeric(adults_estimates2020[1, 3]) + starting_pm2*as.numeric(adults_estimates2020[3, 3])),
                         as.numeric(adults_estimates2020[6,3])/(1+seeded_am2*as.numeric(adults_estimates2020[2, 3]) + starting_pm2*as.numeric(adults_estimates2020[4, 3]))))


d<-ggplot(data=adults2020_plo, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2020b, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2020 w/starting_p=1")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
e<-ggplot(data=adults2020_pmean, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2020b, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2020 w/starting_p=5")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
f<-ggplot(data=adults2020_phi, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2020b, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2020 w/starting_p=9")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()


#visualize conditional plots first with all data
ggarrange(a, b, c, d, e, f, nrow=2, ncol=3, common.legend=T)

#2021----
adults2021b<-mutate(adults2021, warmtrt=ifelse(warmtrt==0, "amb", "warm"))

### predict over pm2 (hold seeded_a steady)
adults2021_ahi <- expand.grid(
  seeded_am2 = 1800 # near the max of seeded_am2
  ,starting_pm2 = as.numeric(seq(0,9, length.out=200))
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2021[5,3])/(1+seeded_am2*as.numeric(adults_estimates2021[1, 3]) + starting_pm2*as.numeric(adults_estimates2021[3, 3])),
                         as.numeric(adults_estimates2021[6,3])/(1+seeded_am2*as.numeric(adults_estimates2021[2, 3]) + starting_pm2*as.numeric(adults_estimates2021[4, 3]))))

adults2021_amean <- expand.grid(
  seeded_am2 = mean(annuals$seeded_am2, na.rm=TRUE) #0
  ,starting_pm2 = seq(0,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2021[5,3])/(1+seeded_am2*as.numeric(adults_estimates2021[1, 3]) + starting_pm2*as.numeric(adults_estimates2021[3, 3])),
                         as.numeric(adults_estimates2021[6,3])/(1+seeded_am2*as.numeric(adults_estimates2021[2, 3]) + starting_pm2*as.numeric(adults_estimates2021[4, 3]))))


adults2021_alo <- expand.grid(
  seeded_am2 = 75 #above the min
  ,starting_pm2 = seq(0,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2021[5,3])/(1+seeded_am2*as.numeric(adults_estimates2021[1, 3]) + starting_pm2*as.numeric(adults_estimates2021[3, 3])),
                         as.numeric(adults_estimates2021[6,3])/(1+seeded_am2*as.numeric(adults_estimates2021[2, 3]) + starting_pm2*as.numeric(adults_estimates2021[4, 3]))))



a<-ggplot(data=adults2021_alo, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2021b, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2021 w/seeded_a=75")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
b<-ggplot(data=adults2021_amean, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2021b, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2021 w/seeded_a=685")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
c<-ggplot(data=adults2021_ahi, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2021b, aes(x=starting_pm2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2021 w/seeded_a=1800")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")


# predict over seeded_a (hold pm2 steady)

adults2021_phi <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 9 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2021[5,3])/(1+seeded_am2*as.numeric(adults_estimates2021[1, 3]) + starting_pm2*as.numeric(adults_estimates2021[3, 3])),
                         as.numeric(adults_estimates2021[6,3])/(1+seeded_am2*as.numeric(adults_estimates2021[2, 3]) + starting_pm2*as.numeric(adults_estimates2021[4, 3]))))

adults2021_pmean <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 5
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2021[5,3])/(1+seeded_am2*as.numeric(adults_estimates2021[1, 3]) + starting_pm2*as.numeric(adults_estimates2021[3, 3])),
                         as.numeric(adults_estimates2021[6,3])/(1+seeded_am2*as.numeric(adults_estimates2021[2, 3]) + starting_pm2*as.numeric(adults_estimates2021[4, 3]))))

adults2021_plo <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 1 
  ,warmtrt = c("amb","warm")
)%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(adults_estimates2021[5,3])/(1+seeded_am2*as.numeric(adults_estimates2021[1, 3]) + starting_pm2*as.numeric(adults_estimates2021[3, 3])),
                         as.numeric(adults_estimates2021[6,3])/(1+seeded_am2*as.numeric(adults_estimates2021[2, 3]) + starting_pm2*as.numeric(adults_estimates2021[4, 3]))))


d<-ggplot(data=adults2021_plo, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2021b, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2021 w/starting_p=1")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
e<-ggplot(data=adults2021_pmean, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2021b, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2021 w/starting_p=5")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
f<-ggplot(data=adults2021_phi, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=adults2021b, aes(x=seeded_am2, y=fecundity), shape=1, width=.25)+
  ylab("Adult fecundity 2021 w/starting_p=9")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()


#visualize conditional plots first with all data
ggarrange(a, b, c, d, e, f, nrow=2, ncol=3, common.legend=T)

# seedlings ----
#2020 ----
seedlings2020b<-mutate(seedlings2020, warmtrt=ifelse(warmtrt==0, "amb", "warm"))

### predict over pm2 (hold seeded_a&seeded_s steady)
seedlings2020_ashi <- expand.grid(
  seeded_am2 = 1800, # near the max of seeded_am2
  seeded_sm2 = 4200 # near the max of seeded_am2
  ,starting_pm2 = as.numeric(seq(0,9, length.out=200))
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(seedlings_estimates2020[7,3])/ # as.numeric(seedlings_estimates2020[7,3])*seeded_sm2/
                           (1+seeded_am2*as.numeric(seedlings_estimates2020[1, 3]) + 
                              starting_pm2*as.numeric(seedlings_estimates2020[3, 3])+ 
                              seeded_sm2*as.numeric(seedlings_estimates2020[5, 3])),
                         as.numeric(seedlings_estimates2020[8,3])/ # as.numeric(seedlings_estimates2020[8,3])*seeded_sm2/
                            (1+seeded_am2*as.numeric(seedlings_estimates2020[2, 3]) + 
                               starting_pm2*as.numeric(seedlings_estimates2020[4, 3])+ 
                               seeded_sm2*as.numeric(seedlings_estimates2020[6, 3]))))

seedlings2020_asmean <- expand.grid(
  seeded_am2 = 150, # near the mean of seeded_am2
  seeded_sm2 = 2800 # near the mean of seeded_sm2
  ,starting_pm2 = as.numeric(seq(0,9, length.out=200))
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(seedlings_estimates2020[7,3])/ # as.numeric(seedlings_estimates2020[7,3])*seeded_sm2/
                           (1+seeded_am2*as.numeric(seedlings_estimates2020[1, 3]) + 
                              starting_pm2*as.numeric(seedlings_estimates2020[3, 3])+ 
                              seeded_sm2*as.numeric(seedlings_estimates2020[5, 3])),
                         as.numeric(seedlings_estimates2020[8,3])/ # as.numeric(seedlings_estimates2020[8,3])*seeded_sm2/
                           (1+seeded_am2*as.numeric(seedlings_estimates2020[2, 3]) + 
                              starting_pm2*as.numeric(seedlings_estimates2020[4, 3])+ 
                              seeded_sm2*as.numeric(seedlings_estimates2020[6, 3]))))


seedlings2020_aslo <- expand.grid(
  seeded_am2 = 75,  #above the min
  seeded_sm2 = 2200
  ,starting_pm2 = seq(0,9, length.out=200)
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                         as.numeric(seedlings_estimates2020[7,3])/ # as.numeric(seedlings_estimates2020[7,3])*seeded_sm2/
                           (1+seeded_am2*as.numeric(seedlings_estimates2020[1, 3]) + 
                              starting_pm2*as.numeric(seedlings_estimates2020[3, 3])+ 
                              seeded_sm2*as.numeric(seedlings_estimates2020[5, 3])),
                         as.numeric(seedlings_estimates2020[8,3])/ # as.numeric(seedlings_estimates2020[8,3])*seeded_sm2/
                           (1+seeded_am2*as.numeric(seedlings_estimates2020[2, 3]) + 
                              starting_pm2*as.numeric(seedlings_estimates2020[4, 3])+ 
                              seeded_sm2*as.numeric(seedlings_estimates2020[6, 3]))))


a<-ggplot(data=seedlings2020_aslo, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=seedlings2020b, aes(x=starting_pm2, y=survivors/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival 2020 w/seeded_a=75&s=2200")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
b<-ggplot(data=seedlings2020_asmean, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=seedlings2020b, aes(x=starting_pm2, y=survivors/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival 2020 w/seeded_a=685&s=2800")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
c<-ggplot(data=seedlings2020_ashi, aes(x = starting_pm2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=seedlings2020b, aes(x=starting_pm2, y=survivors/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival 2020 w/seeded_a=1800&s=4200")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")


# predict over seeded_a (hold sm2 and pm2 steady)

seedlings2020_sphi <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 9,
  seeded_sm2 = 4200 
  ,warmtrt = c("amb","warm"))%>%
  mutate(estimate=ifelse(warmtrt=="amb", 
                           as.numeric(seedlings_estimates2020[7,3])/ # as.numeric(seedlings_estimates2020[7,3])*seeded_sm2/
                             (1+seeded_am2*as.numeric(seedlings_estimates2020[1, 3]) + 
                                starting_pm2*as.numeric(seedlings_estimates2020[3, 3])+ 
                                seeded_sm2*as.numeric(seedlings_estimates2020[5, 3])),
                           as.numeric(seedlings_estimates2020[8,3])/ # as.numeric(seedlings_estimates2020[8,3])*seeded_sm2/
                             (1+seeded_am2*as.numeric(seedlings_estimates2020[2, 3]) + 
                                starting_pm2*as.numeric(seedlings_estimates2020[4, 3])+ 
                                seeded_sm2*as.numeric(seedlings_estimates2020[6, 3]))))
  
seedlings2020_spmean <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 5,
  seeded_sm2 = 2400 
  ,warmtrt = c("amb","warm"))%>%
    mutate(estimate=ifelse(warmtrt=="amb", 
                           as.numeric(seedlings_estimates2020[7,3])/ # as.numeric(seedlings_estimates2020[7,3])*seeded_sm2/
                             (1+seeded_am2*as.numeric(seedlings_estimates2020[1, 3]) + 
                                starting_pm2*as.numeric(seedlings_estimates2020[3, 3])+ 
                                seeded_sm2*as.numeric(seedlings_estimates2020[5, 3])),
                           as.numeric(seedlings_estimates2020[8,3])/ # as.numeric(seedlings_estimates2020[8,3])*seeded_sm2/
                             (1+seeded_am2*as.numeric(seedlings_estimates2020[2, 3]) + 
                                starting_pm2*as.numeric(seedlings_estimates2020[4, 3])+ 
                                seeded_sm2*as.numeric(seedlings_estimates2020[6, 3]))))
  
seedlings2020_splo <- expand.grid(
  seeded_am2 = seq(75,1800, length.out=100)
  ,starting_pm2 = 1 ,
  seeded_sm2 = 2200 
  ,warmtrt = c("amb","warm"))%>%
    mutate(estimate=ifelse(warmtrt=="amb", 
                           as.numeric(seedlings_estimates2020[7,3])/ # as.numeric(seedlings_estimates2020[7,3])*seeded_sm2/
                             (1+seeded_am2*as.numeric(seedlings_estimates2020[1, 3]) + 
                                starting_pm2*as.numeric(seedlings_estimates2020[3, 3])+ 
                                seeded_sm2*as.numeric(seedlings_estimates2020[5, 3])),
                           as.numeric(seedlings_estimates2020[8,3])/ # as.numeric(seedlings_estimates2020[8,3])*seeded_sm2/
                             (1+seeded_am2*as.numeric(seedlings_estimates2020[2, 3]) + 
                                starting_pm2*as.numeric(seedlings_estimates2020[4, 3])+ 
                                seeded_sm2*as.numeric(seedlings_estimates2020[6, 3]))))

d<-ggplot(data=seedlings2020_splo, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=seedlings2020b, aes(x=seeded_am2, y=survivors/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival 2020 w/p=1&s=2200")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
e<-ggplot(data=seedlings2020_spmean, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=seedlings2020b, aes(x=seeded_am2, y=survivors/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival 2020 w/p=5&s=2800")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")
f<-ggplot(data=seedlings2020_sphi, aes(x = seeded_am2, y = estimate,  color=as.factor(warmtrt))) + 
  geom_point()+
  geom_jitter(data=seedlings2020b, aes(x=seeded_am2, y=survivors/seeded_sm2), shape=1, width=.25)+
  ylab("Seedling survival 2020 w/p=9&s=4200")+ scale_color_manual(values=c("blue", "red"))+
  theme_classic()+theme(legend.position="none")

#still need to do holding perennials and annuals still and varying seedlings. 

#visualize conditional plots first with all data
ggarrange(a, b, c, d, e, f, nrow=2, ncol=3, common.legend=T)
