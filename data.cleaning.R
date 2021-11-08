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
  mutate(sub=toupper(sub))%>%
  dplyr::select(1:7, 10:18)%>%
  mutate(area=ifelse(is.na(area), (as.numeric(diameter)/2)^2*3.14, area))%>%
  mutate(area=ifelse(is.na(area), (as.numeric(circumference)/3.14/2)^2*3.14, area))%>%
  dplyr::select(-diameter, -circumference)

startingphyt<-read_csv("starting_plugs.csv")%>%
  mutate(sub=toupper(sub))

### DENSITY ####
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
vegplot2021<-read_csv("key_cover2021.csv")%>%
  mutate(sub=toupper(sub))%>%
  mutate(comptrt=ifelse(comptrt=="n", "none", ifelse(comptrt=="a", "annuals", 
                                                     ifelse(comptrt=="p", "adult perennials", 
                                                            ifelse(comptrt=="s", "seedling perennials",
                                                                   ifelse(comptrt=="a_p", "annuals+adults", 
                                                                          ifelse(comptrt=="a_s", "annuals+seedlings", "seedlings+adults")))))))%>%
  mutate(plotid=paste(row, column, sub, sep="."))%>%
  mutate(time="spring2021")
vegplot2021$comptrt <- factor(vegplot$comptrt, levels = c("none", "annuals", "seedling perennials", "adult perennials", "annuals+seedlings", "annuals+adults", "seedlings+adults"))
#working on 2021


# CREATE PLOTKEY
plotkey<-select(vegplot, plotid, block, warmtrt, comptrt)

vegplot2020<-select(vegplot, 18, 19, 1:14)
vegplot2020[is.na(vegplot2020)] <- 0

density_spring20 <- right_join(plotkey, dplyr::select(vegplot2020, plotid, time, per_a, per_p, per_s, count_a, count_p, count_s))%>%
  mutate(am2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_a, count_a/.66))%>%
  mutate(pm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_p, count_p/.66))%>%
  mutate(sm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_s, count_s/.66))
density_spring21 <- right_join(plotkey, dplyr::select(vegplot2021, plotid, time, per_a, per_p, per_na, per_s, count_a, count_p, count_na, count_dam_p, count_s, gopher_spring))%>%
  mutate(am2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_a, count_a/.66))%>%
  mutate(pm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_p, count_p/.66))%>%
  mutate(dpm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_dam_p, count_dam_p/.66))%>%
  mutate(nam2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_na, count_na/.66))%>%
  mutate(sm2=ifelse(comptrt%in%c("none", "annuals", "adult perennials", "seedling perennials"), count_s, count_s/.66))
#am2, pm2 and sm2 are the relevant data here. 
#add .g a gopher correction of area for each (NO, it's based on who is actually there)


### SUMMER SURVIVAL ####
#create seedling_sumsur2020, calculation of seedling summer survival in 2020
seedling_sumsur2020<-select(vegplot, 18, 19, 10, 15, 16)%>%
  mutate(time="fall2020")%>%
  filter(!is.na(count_s), gopher_ss_correction!=1, count_s!=0)%>% #get rid of plots w/o seedlings, 100% gopher damage
  mutate(spring=count_s, fall=`20fall_s`)%>% 
  mutate(fall.g=fall/(1-gopher_ss_correction))%>% # adjust summer survival for the proportion of the plot that was damaged by gophers
  select(-count_s, -`20fall_s`)%>%
  mutate(fall.g=ifelse(fall.g>spring, spring, fall.g))%>% #in some cases my adjustment above put fall>spring, this is to max out at 100% survival
  mutate(sumsur=fall.g/spring)%>%
  mutate(gopher=gopher_ss_correction, time=2020)%>%
  select(1:2, 4:8)


seedling_sumsur2021<-select(vegplot2021, plotid, time, count_s, count_s_fall, fall_s_g, gopher_fall)%>%
  mutate(spring=ifelse(!is.na(fall_s_g)&is.na(gopher_fall), `count_s`-fall_s_g, count_s), fall=`count_s_fall`)%>% #if I know a seedling was gone in fall due to gopher, and there is no other accounting for gopher, subtract it out from the beginning
  mutate(gopher_fall=ifelse(is.na(gopher_fall), 0, gopher_fall))%>%
  filter(!is.na(count_s), count_s!=0)%>% #get rid of plots w/o seedlings, 100% gopher damage
  mutate(fall.g=fall/(1-gopher_fall))%>% # adjust summer survival for the proportion of the plot that was damaged by gophers
select(-count_s, -`count_s_fall`)%>%
  mutate(fall.g=ifelse(fall.g>spring, spring, fall.g))%>% #in some cases my adjustment above put fall>spring, this is to max out at 100% survival
  mutate(sumsur=fall.g/spring)%>%
  mutate(gopher=gopher_fall, time=2021)%>%
  select(1, 2, 5:9)

seedling_sumsur<-rbind(seedling_sumsur2020, seedling_sumsur2021)

### Seedling summer survival: 
# phytometers and background within each plot
#spring20_s and fall20_s.g are the relevant data here
seedling_sumsur<-right_join(plotkey, seedling_sumsur)

rm(seedling_sumsur2020, seedling_sumsur2021)


#### SPRING SURVIVAL ####
spr_sur2020<-select(vegplot2020, 1:8, 12, 14)%>%
  filter(comptrt!="none"&comptrt!="adult perennials")%>%
  mutate(seeded_s=ifelse(comptrt=="seedling perennials", 750*6, #weighed 1g of seed and counted
                         ifelse(comptrt=="seedlings+adults"|comptrt=="annuals+seedlings", 750*2, 0)))%>%
  mutate(seeded_a=ifelse(comptrt=="annuals", 310*6, # from online lit: 3.2-5g/1000 seeds
                         ifelse(comptrt=="annuals+adults"|comptrt=="annuals+seedlings", 310*2, 0)))%>%
  mutate(time=2020)
 
spr_sur2021<-select(vegplot2021, plotid, time, 1:6, count_s, count_a, gopher_spring)%>%
  filter(comptrt!="none"&comptrt!="adult perennials")%>%
  mutate(seeded_s=ifelse(comptrt=="seedling perennials", 750*6, #weighed 1g of seed and counted
                         ifelse(comptrt=="seedlings+adults"|comptrt=="annuals+seedlings", 750*2, 0)))%>%
  mutate(seeded_a=ifelse(comptrt=="annuals", 310*6,
                         ifelse(comptrt=="annuals+adults"|comptrt=="annuals+seedlings", 310*2, 0)))%>%
  mutate(time=2021)%>%
  mutate(seeded_s=seeded_s*gopher_spring, seeded_a=seeded_a*gopher_spring)%>%
  select(-gopher_spring)
  
sprsur<-rbind(spr_sur2020, spr_sur2021)
rm(spr_sur2020, spr_sur2021)

#Annual spring survival
annual_sprsur<-filter(sprsur, seeded_a!=0)%>%
  dplyr::select(-seeded_s, -count_s)%>%
  mutate(survival=count_a/seeded_a)

#Seedling spring survival
seedling_sprsur <-filter(sprsur, seeded_s!=0)%>%
  dplyr::select(-seeded_a, -count_a)%>%
  mutate(survival=count_s/seeded_s)


### FECUNDITY ####
#read_csv("florets.csv")
#phytometers1<-left_join(phytometers, select(startingphyt, -date))%>%
#  mutate(growth=veg_height-starting_height)%>%
##  mutate(widening=circumference-starting_cir)%>%
#  mutate(veg_height=ifelse(is.na(veg_height), rep_height, veg_height))%>%
##  mutate(type=substr(id, 1, 1), id=substr(id, 2, 2))%>%
#  mutate(plotid=paste(row, column, sub, sep="."))




phytometers1<-phytometers%>%
  group_by(year_data, year_planted, row, column, sub)%>%
  mutate(plotid=paste(row, column, sub, sep="."))%>%
  mutate(tillers=as.numeric(ifelse(is.na(tillers), 0, tillers)))%>%
  mutate(type=(ifelse(is.na(type), "plug", type)))%>%
  mutate(type2=ifelse(type=="adult"&year_data==2021, `2021_spring`, ifelse(type=="adult"&year_data==2020, "plug", type)))%>%
  mutate(type2=paste(type2, year_data, sep="_"))%>%
  mutate(type2=ifelse(type2=="alive_2021", "plug_2021", type2))%>%
  mutate(type2=ifelse(type2=="partial_2021", "gopherplug_2021", type2))
phytometers1$type2 <- factor(phytometers1$type2, levels = c("plug_2019", "plug_2020", "plug_2021", "gopherplug_2021", "newadult_2021"))

#phytometers2, to compare differences between different types of adults (over time, damaged, new)
phytometers2<-left_join(phytometers1, plotkey)#%>%
  mutate(growth=veg_height-starting_height)%>%
  mutate(widening=circumference-starting_cir)%>%
  mutate(veg_height=ifelse(is.na(veg_height), rep_height, veg_height))%>%
  #  select(-1, -8, -date, -starting_height, -starting_cir)%>%
  mutate(type=substr(id, 1, 1), id=substr(id, 2, 2))%>%
  
#visualize plugs, new adults, and existing adults (damaged, vs good)
ggplot(subset(phytometers2, type=="adult"|type=="plug"|type=="newadult"), aes(x=type2, y=tillers)) + geom_boxplot(aes(fill=comptrt))
ggplot(subset(phytometers2, type=="adult"|type=="plug"|type=="newadult"), aes(x=type2)) + geom_histogram(stat="count")  


fecundity<-select(phytometers1, plotid, type, id, tillers, type2)%>% 
  filter(type!="seedling")%>%
  mutate(type=ifelse(type2=="gopherplug_2021", "adult_g", type))%>%
  mutate(type=ifelse(is.na(type), "annual", type))%>%
  filter(!is.na(tillers))%>%
  select(-type2)%>%
  mutate(seeds=ifelse(type=="annual", tillers*18*10.6, tillers*25.5*6)) #lolium: 30 spikelets w/ up to 8-20 florets, festuca 25 spikelets with 6 florets
rm(phytometers1)

#ANNUALS DATA
# annual seed production per capita
fecundity_plot_a2020<-select(vegplot2020, plotid, time, count_a, til_a)%>%
  filter(count_a>0)%>%
  mutate(pc_til=til_a/count_a)%>%
  mutate(seeds=pc_til*18*10.6) #annuals
fecundity_plot_a2021<-select(vegplot2021, plotid, time, count_a, til_a)%>%
  filter(count_a>0)%>%
  mutate(pc_til=til_a/count_a)%>%
  mutate(seeds=pc_til*18*10.6) #annuals
#combine both years
fecundity_plot_a<-rbind(fecundity_plot_a2020, fecundity_plot_a2021)
#annual phytometers
fecundity_phyt_a <- filter(fecundity, type=="annual") # phytometers

#PERENNIAL FECUNDITY
# plug seed production per capita
fecundity_phyt_p <-filter(fecundity, type!="annual") #phytometers

rm(phytometers, vegplot, vegplot2020, vegplot2021, fecundity_plot_a2020, fecundity_plot_a2021, phytometers2, sprsur2020s, sumsur2020s)


