
library(tidyverse)
library(ggpubr)
##FN for Calculating SE
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

theme_set(theme_classic())
#### VEG  Read-in/Cleanup ####

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
                                                                          ifelse(comptrt=="a_s", "annuals+seedlings", "seedlings+adults")))))))%>%
  mutate(plotid=paste(row, column, sub, sep="."))%>%
  mutate(time="spring2020")
vegplot$comptrt <- factor(vegplot$comptrt, levels = c("none", "annuals", "seedling perennials", "adult perennials", "annuals+seedlings", "annuals+adults", "seedlings+adults"))

# create vegplot2020, cover and counts of each 'species' in each plot at springtime
vegplot2020<-select(vegplot, 18, 19, 1:12)
vegplot2020[is.na(vegplot2020)] <- 0

#create seedling_sumsur2020, calculation of seedling summer survival in 2020
seedling_sumsur2020<-select(vegplot, 18, 19, 10, 15, 16)%>%
  mutate(time="fall2020")%>%
  filter(!is.na(count_s), gopher_ss_correction!=1, count_s!=0)%>%
  mutate(spring20_s=count_s, fall20_s=`20fall_s`)%>%
  mutate(fall20_s.g=fall20_s/(1-gopher_ss_correction))%>% # adjust summer survival for the proportion of the plot that was damaged by gophers
  select(-count_s, `20fall_s`)%>%
  mutate(fall20_s.g=ifelse(fall20_s.g>spring20_s, spring20_s, fall20_s.g))%>% #in some cases my adjustment above put fall>spring, this is to max out at 100% survival
  mutate(sumsur=fall20_s.g/spring20_s)
  

#create plotkey, used to match 
plotkey<-select(vegplot, 18, 3, 5, 6 )


#germination and survival
spr_sur2020<-select(vegplot2020, 1:8, 12, 14)%>%
  filter(comptrt!="none"&comptrt!="adult perennials")%>%
  mutate(seeded_s=ifelse(comptrt=="seedling perennials", 992*6, #calculated from seeds per pound online.  seems pretty high to me...
         ifelse(comptrt=="seedlings+adults"|comptrt=="annuals+seedlings", 992*2, 0)))%>%
  mutate(seeded_a=ifelse(comptrt=="annuals", 478*6,
                         ifelse(comptrt=="annuals+adults"|comptrt=="annuals+seedlings", 478*2, 0)))%>%
  mutate(spring20_s=count_s, spring20_a=count_a)%>%
  select(-count_s, -count_a)
  


#clean up phytometers.  Only to be used for fecundity. 

phytometers1<-left_join(phytometers, select(startingphyt, -date))%>%
  mutate(growth=veg_height-starting_height)%>%
  mutate(widening=circumference-starting_cir)%>%
  mutate(veg_height=ifelse(is.na(veg_height), rep_height, veg_height))%>%
#  select(-1, -8, -date, -starting_height, -starting_cir)%>%
  mutate(type=substr(id, 1, 1), id=substr(id, 2, 2))%>%
  mutate(plotid=paste(row, column, sub, sep="."))

phytometers0<-left_join(phytometers1, plotkey)%>%
  mutate(comptrt=ifelse(comptrt=="none", "none", ifelse(comptrt=="a", "annuals", 
                                                     ifelse(comptrt=="p", "adult perennials", 
                                                            ifelse(comptrt=="s", "seedling perennials",
                                                                    ifelse(comptrt=="a_p", "annuals+adults", 
                                                                           ifelse(comptrt=="a_s", "annuals+seedlings", "seedlings+adults")))))))
phytometers0$comptrt <- factor(phytometers0$comptrt, levels = c("none", "annuals", "seedling perennials", "adult perennials", "annuals+seedlings", "annuals+adults", "seedlings+adults"))

# to do - correct tillers to florets. this year I want to count seeds...
    # to be more accurate. germination is so low that i should. florets can maybe have >1 seed in them
#read_csv("florets.csv")

fecundity<-select(phytometers0, plotid, block, warmtrt, comptrt, type, id, tillers)%>%
  filter(type!="s")%>%
  filter(!is.na(tillers))%>%
  mutate(seeds=ifelse(type=="a", tillers*30, tillers*30))


fecundity_plot<-select(vegplot, plotid, time, count_a, til_a)%>%
  filter(count_a>0)%>%
  mutate(pc_til=til_a/count_a)%>%
  mutate(seeds=pc_til*30)



### Old Visualizations
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
