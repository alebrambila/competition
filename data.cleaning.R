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

### DENSITY ###
### Density is measured as cover and counts/m2.
### I need to correct counts to density per meter squared for some plots.  Full plots are 1m2 and 50:50 plots are .66m2. 
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

vegplot2020<-select(vegplot, 18, 19, 1:12)
vegplot2020[is.na(vegplot2020)] <- 0

density_spring20 <- right_join(plotkey, dplyr::select(vegplot2020, plotid, time, per_a, per_p, per_s, count_a, count_p, count_s))%>%
  mutate(am2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_a, count_a/.66))%>%
  mutate(pm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_p, count_p/.66))%>%
  mutate(sm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_s, count_s/.66))
#density_spring21 <- ### TO BE GATHERED IN JUNE 2021
#am2, pm2 and sm2 are the relevant data here. 


### PLOTKEY ###
plotkey<-select(vegplot, plotid, block, warmtrt, comptrt)

### SUMMER SURVIVAL ###
#create seedling_sumsur2020, calculation of seedling summer survival in 2020
seedling_sumsur2020<-select(vegplot, 18, 19, 10, 15, 16)%>%
  mutate(time="fall2020")%>%
  filter(!is.na(count_s), gopher_ss_correction!=1, count_s!=0)%>%
  mutate(spring20_s=count_s, fall20_s=`20fall_s`)%>%
  mutate(fall20_s.g=fall20_s/(1-gopher_ss_correction))%>% # adjust summer survival for the proportion of the plot that was damaged by gophers
  select(-count_s, `20fall_s`)%>%
  mutate(fall20_s.g=ifelse(fall20_s.g>spring20_s, spring20_s, fall20_s.g))%>% #in some cases my adjustment above put fall>spring, this is to max out at 100% survival
  mutate(sumsur=fall20_s.g/spring20_s)

### Seedling summer survival: 
# phytometers and background within each plot
#spring20_s and fall20_s.g are the relevant data here
sumsur2020s<-right_join(plotkey, seedling_sumsur2020)
#sumsur2021s <- ### TO BE GATHERED JUNE AND SEPTEMBER 2021


#### SPRING SURVIVAL ###
spr_sur2020<-select(vegplot2020, 1:8, 12, 14)%>%
  filter(comptrt!="none"&comptrt!="adult perennials")%>%
  mutate(seeded_s=ifelse(comptrt=="seedling perennials", 750*6, #weighed 1g of seed and counted
         ifelse(comptrt=="seedlings+adults"|comptrt=="annuals+seedlings", 750*2, 0)))%>%
  mutate(seeded_a=ifelse(comptrt=="annuals", 244*6, # from online lit: 3.2-5g/1000 seeds
                         ifelse(comptrt=="annuals+adults"|comptrt=="annuals+seedlings", 244*2, 0)))%>%
  mutate(spring20_s=count_s, spring20_a=count_a)%>%
  select(-count_s, -count_a)

#Annual spring survival
sprsur2020a<-filter(spr_sur2020, seeded_a!=0)%>%
  dplyr::select(-seeded_s, -spring20_s)%>%
  mutate(sprsur_a=spring20_a/seeded_a)

#Seedling spring survival
sprsur2020s <-filter(spr_sur2020, seeded_s!=0)%>%
  dplyr::select(-seeded_a, -spring20_a)%>%
  mutate(sprsur_s=spring20_s/seeded_s)


### FECUNDITY ###
# to do - correct tillers to florets. this year I want to count seeds...
# to be more accurate. germination is so low that i should. florets can maybe have >1 seed in them
#read_csv("florets.csv")
phytometers1<-left_join(phytometers, select(startingphyt, -date))%>%
  mutate(growth=veg_height-starting_height)%>%
  mutate(widening=circumference-starting_cir)%>%
  mutate(veg_height=ifelse(is.na(veg_height), rep_height, veg_height))%>%
  mutate(type=substr(id, 1, 1), id=substr(id, 2, 2))%>%
  mutate(plotid=paste(row, column, sub, sep="."))

fecundity<-select(phytometers1, plotid, type, id, tillers)%>% 
  filter(type!="s")%>%
  filter(!is.na(tillers))%>%
  mutate(seeds=ifelse(type=="a", tillers*30*14, tillers*30*5.5)) #lolium: 30 spikelets w/ up to 8-20 florets, festuca 4-7/spikelet
rm(phytometers1)


fecundity_plot<-select(vegplot, plotid, time, count_a, til_a)%>%
  filter(count_a>0)%>%
  mutate(pc_til=til_a/count_a)%>%
  mutate(seeds=pc_til*30*10)

# annual seed production per capita
fecundity_phyt_a <- filter(fecundity, type=="a") # phytometers
fecundity_plot_a<-fecundity_plot #plot level

# plug seed production per capita
fecundity_phyt_p <-filter(fecundity, type=="p") #phytometers


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
