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

                                                                