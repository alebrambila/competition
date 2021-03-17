### Created 8/19/19
### Analysis of ITEX warming chamber trial
### Treatments include standard ITEX (itx), thermal mass ITEX (itm), and ambient (c). 
library(tidyverse)
library(ggpubr)
##FN for Calculating SE
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

theme_set(theme_classic())
#### TEMPERATURES ####
### Read in the data
a2<-read_csv("./temps/summer_20/plot2air.csv", skip=14)%>%
  mutate(plot=2, trt="control",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s2<-read_csv("./temps/summer_20/plot2soil.csv", skip=14)%>%
  mutate(plot=2, trt="control",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a3<-read_csv("./temps/summer_20/plot3air.csv", skip=14)%>%
  mutate(plot=3, trt="control",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s3<-read_csv("./temps/summer_20/plot3soil.csv", skip=14)%>%
  mutate(plot=3, trt="control",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a5<-read_csv("./temps/summer_20/plot5air.csv", skip=14)%>%
  mutate(plot=5, trt="control",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s5<-read_csv("./temps/summer_20/plot5soil.csv", skip=14)%>%
  mutate(plot=5, trt="control",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a7<-read_csv("./temps/summer_20/plot7air.csv", skip=14)%>%
  mutate(plot=7, trt="control",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s7<-read_csv("./temps/summer_20/plot7soil.csv", skip=14)%>%
  mutate(plot=7, trt="control",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a9<-read_csv("./temps/summer_20/plot9air.csv", skip=14)%>%
  mutate(plot=9, trt="warm",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s9<-read_csv("./temps/summer_20/plot9soil.csv", skip=14)%>%
  mutate(plot=9, trt="warm",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a11<-read_csv("./temps/summer_20/plot11air.csv", skip=14)%>%
  mutate(plot=11, trt="warm",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s11<-read_csv("./temps/summer_20/plot11soil.csv", skip=14)%>%
  mutate(plot=11, trt="warm",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a13<-read_csv("./temps/summer_20/plot13air.csv", skip=14)%>%
  mutate(plot=13, trt="warm",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s13<-read_csv("./temps/summer_20/plot13soil.csv", skip=14)%>%
  mutate(plot=13, trt="warm",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a15<-read_csv("./temps/summer_20/plot15air.csv", skip=14)%>%
  mutate(plot=15, trt="warm",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s15<-read_csv("./temps/summer_20/plot15soil.csv", skip=14)%>%
  mutate(plot=15, trt="warm",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a17<-read_csv("./temps/summer_20/plot17air.csv", skip=14)%>%
  mutate(plot=17, trt="control",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s17<-read_csv("./temps/summer_20/plot17soil.csv", skip=14)%>%
  mutate(plot=17, trt="control",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a22<-read_csv("./temps/summer_20/plot22air.csv", skip=14)%>%
  mutate(plot=22, trt="warm",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s22<-read_csv("./temps/summer_20/plot22soil.csv", skip=14)%>%  #filter this until jan 2020 (some weird stuff)
  mutate(plot=22, trt="warm",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  slice(700:n())%>%
  select(-Unit)
a23<-read_csv("./temps/summer_20/plot23air.csv", skip=14)%>%
  mutate(plot=23, trt="control",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s23<-read_csv("./temps/summer_20/plot23soil.csv", skip=14)%>%
  mutate(plot=23, trt="control",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
a24<-read_csv("./temps/summer_20/plot24air.csv", skip=14)%>%
  mutate(plot=24, trt="warm",  medium="air")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)
s24<-read_csv("./temps/summer_20/plot24soil.csv", skip=14)%>%
  mutate(plot=24, trt="warm",  medium="soil")%>%
  rename(time="Date/Time", temp=Value)%>%
  select(-Unit)


fullset<-rbind(a2, s2, a3, s3, a5, s5, a7, s7, a9, s9, a11, s11,
               a13, s13, a15, s15, a17, s17, a22, s22, a23, s23, a24, s24)%>%
  mutate(ztemp0=scale(temp))%>%
  group_by(time)%>%
  mutate(ztemp=scale(temp))%>%
  mutate(daynight=ifelse(ztemp0>0, "day", "night"))%>%
  ungroup()%>%
  mutate(timparse=parse_datetime(time, "%m/%d/%y %I:%M:%S %p"))%>%
  mutate(timparse2=as.Date(timparse))%>%
  mutate(mdy=as.factor(format(timparse, "%m%d%y")))

rm(a2, s2, a3, s3, a5, s5, a7, s7, a9, s9, a11, s11,
   a13, s13, a15, s15, a17, s17, a22, s22, a23, s23, a24, s24)

timelord<-unique(select(fullset, timparse2, mdy))

#jitter of all data
ggplot(subset(fullset), aes(timparse2, temp)) +
  geom_jitter(aes(color=as.factor(trt)))+
  scale_x_date(date_breaks="months", date_labels="%b")+
  facet_grid(medium~plot)# plot 22 has some issues in December - replot without it

ggplot(subset(fullset, plot!=22), aes(timparse2, temp)) +
  geom_jitter(aes(color=as.factor(trt)))+
  scale_x_date(date_breaks="months", date_labels="%b")+facet_wrap(~medium)


fullset0<-fullset%>%
  group_by(mdy, timparse2, trt, plot, medium)%>%
  summarize(maxtemp=as.numeric(max(temp)), mintemp=min(temp), meantemp=mean(temp))%>%
  #filter(temp==maxtemp|temp==mintemp|temp==meantemp)%>%
  gather("maxmin", "temp", maxtemp, mintemp, meantemp)#%>%
  group_by(plot, trt, medium, timparse2, mdy, maxmin, temp2)#%>%
  summarize()
fullset1<-fullset0%>%
  ungroup()%>%
  group_by(temp, trt, plot, medium, mdy, timparse2, maxmin)%>%
  summarize()



#actual max min temps over the time period, but lots of var becase seasonality matters
ggplot(fullset1, aes(x=maxmin, y=temp)) +geom_boxplot(aes(fill=trt)) +facet_wrap(~medium, scales="free")


#do differences within each day (ISSUE - 6w6c, can't subtract for diff. need to do by mean?)
fullset2<-fullset1%>%
  group_by(trt, medium, maxmin, mdy, timparse2)%>%
  summarize(meanmeantemp=mean(temp))%>%
  spread(trt, meanmeantemp)%>%
  mutate(tempdiff=warm-control)%>%
  select(-control, -warm)%>%
  ungroup()#%>%
  group_by(mdy, timparse2, position, maxmin, trt)%>%
  summarize(tempdiff=mean(tempdiff))


## air and soil temp differences in chamber vs control
ggplot(fullset2, aes(x=maxmin, y=tempdiff)) +
  geom_boxplot(aes()) +
  geom_jitter(aes(color=timparse2), size=.5)+
  facet_wrap(~medium, scales="free")+
  ylab("Deviation in Temp > Control (Celsius)") +
  xlab("")+
  geom_hline(yintercept=0)+scale_color_date(date_breaks="months", date_labels="%b")

ggplot(subset(fullset2), aes((timparse2), tempdiff)) +
  geom_point(aes(color=maxmin), size=.7)+geom_smooth(aes(color=maxmin), method="lm")+facet_wrap(~medium)+
scale_x_date(date_breaks="months", date_labels="%b")+geom_hline(aes(yintercept=0))

rm(fullset0, fullset1, fullset, timelord)

#### VEG  Read-in/Cleanup####

### Read data
phytometers<-read_csv("phytometers.csv")%>%
  mutate(sub=toupper(sub))
startingphyt<-read_csv("starting_plugs.csv")%>%
  mutate(sub=toupper(sub))
vegplot<-read_csv("key_cover.csv")%>%
  mutate(sub=toupper(sub))%>%
  mutate(comptrt=ifelse(comptrt=="n", "none", ifelse(comptrt=="a", "annuals", 
                                                     ifelse(comptrt=="p", "adult perennials", 
                                                            ifelse(comptrt=="s", "seedling perennials",
                                                                   ifelse(comptrt=="a_p", "annuals+adults", 
                                                                          ifelse(comptrt=="a_s", "annuals+seedlings", "seedlings+adults")))))))
vegplot$comptrt <- factor(vegplot$comptrt, levels = c("none", "annuals", "seedling perennials", "adult perennials", "annuals+seedlings", "annuals+adults", "seedlings+adults"))

plotkey<-select(vegplot, 1:6)

#clean up phytometers

phytometers1<-left_join(phytometers, select(startingphyt, -date))%>%
  mutate(growth=veg_height-starting_height)%>%
  mutate(widening=circumference-starting_cir)%>%
  mutate(veg_height=ifelse(is.na(veg_height), rep_height, veg_height))%>%
#  select(-1, -8, -date, -starting_height, -starting_cir)%>%
  mutate(type=substr(id, 1, 1), id=substr(id, 2, 2))

phytometers0<-left_join(phytometers1, plotkey)%>%
  mutate(comptrt=ifelse(comptrt=="none", "none", ifelse(comptrt=="a", "annuals", 
                                                     ifelse(comptrt=="p", "adult perennials", 
                                                            ifelse(comptrt=="s", "seedling perennials",
                                                                    ifelse(comptrt=="a_p", "annuals+adults", 
                                                                           ifelse(comptrt=="a_s", "annuals+seedlings", "seedlings+adults")))))))
phytometers0$comptrt <- factor(phytometers0$comptrt, levels = c("none", "annuals", "seedling perennials", "adult perennials", "annuals+seedlings", "annuals+adults", "seedlings+adults"))


#### Part 1: Phytometer growth and reproduction by competition and warming treatments ####

#perennial growth (height)
ggplot(subset(phytometers0, type=="p"), aes(y=growth, x=comptrt))+geom_boxplot(aes(fill=warmtrt))
#perennial growth (circumference)
ggplot(subset(phytometers0, type=="p"), aes(y=widening, x=comptrt))+geom_boxplot(aes(fill=warmtrt))
#perennial reproduction
ggplot(subset(phytometers0, type=="p"), aes(y=tillers, x=warmtrt))+geom_boxplot(aes(fill=warmtrt)) +facet_grid(~comptrt)

#seedling size*diameter
ggplot(subset(phytometers0, type=="s"), aes(y=veg_height*diameter, x=comptrt))+geom_boxplot(aes(fill=warmtrt))

#annual height
ggplot(subset(phytometers0, type=="a"), aes(y=veg_height, x=comptrt))+geom_boxplot(aes(fill=warmtrt))
#annual rep
ggplot(subset(phytometers0, type=="a"), aes(y=tillers, x=comptrt))+geom_boxplot(aes(fill=warmtrt))



#### Part 2: Plot level count, tillers, survival(TODO) ####

#seedling counts
ggplot(subset(vegplot), aes(y=count_s, x=warmtrt))+geom_boxplot(aes(fill=warmtrt))+geom_jitter(width=.1, size=.5)+facet_wrap(~comptrt, scales="free")

#annual counts
ggplot(subset(vegplot), aes(y=count_a, x=warmtrt))+geom_boxplot(aes(fill=warmtrt))+geom_jitter(width=.1, size=.5)+facet_wrap(~comptrt, scales="free")

#annual tillers
ggplot(subset(vegplot), aes(y=til_a, x=warmtrt))+geom_boxplot(aes(fill=warmtrt))+geom_jitter(width=.1, size=.5)+facet_wrap(~comptrt, scales="free")


#annual tillers
ggplot(subset(vegplot), aes(y=til_p, x=warmtrt))+geom_boxplot(aes(fill=warmtrt))+geom_jitter(width=.1, size=.5)+facet_wrap(~comptrt, scales="free")


### Part 3: Cover ####
vegcov<-select(vegplot, 1:9)%>%
  gather("type", "cover", per_s, per_p, per_a)


ggplot(vegcov, aes(x=type, y=cover)) +geom_boxplot(aes(fill=warmtrt))+facet_wrap(~comptrt, scales="free")



#aggregate survival
#surival<-phytometers0%>%
#  group_by(row, column, sub, block, warmtrt, comptrt, type)%>%
#  summarize(survivors=n())%>%
#  mutate(starters=NA)%>%
#  mutate(starters=ifelse(comptrt=="a"&type=="a", 4,
#                        ifelse(comptrt=="a"&type=="s", 4,
#                               ifelse(comptrt=="a"&type=="p", 2,
#                                     ifelse(comptrt=="s"&type=="s", 4,
#                                           ifelse(comptrt=="s"&type=="a", 4,
#                                                 ifelse(comptrt=="s"& type=="p", 2,
#                                                       ifelse(comptrt=="n"& type=="a", 4,
#                                                             ifelse(comptrt=="n"& type=="s" ,4,
#                                                                   ifelse(comptrt=="n"& type=="p", 2,
#                                                                         ifelse(comptrt=="p"& type=="a", 8,
#                                                                               ifelse(comptrt=="p"& type=="s" ,8,
#                                                                                     ifelse(comptrt=="p"& type=="p", 10, starters)))))))))))))%>%
# mutate(survival=survivors/starters)

#survival
#ggplot(surival, aes(x=type, y=survival))+geom_boxplot(aes(fill=warmtrt))+facet_wrap(~comptrt)
