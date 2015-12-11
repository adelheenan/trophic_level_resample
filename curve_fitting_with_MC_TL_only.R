### fitting the curves - once
### once the curve fitting step is complete - need to put it in a loop that draws from the distribution of trophic levels

### step 1 updating the curve fitting to run for one year

# Matteo Zucchetta
# Marco Anelli Monti
# R script to fitting cumulative curves of biomasses with trophic levels.
# Estimate parameters of curves: lower asymptote, inflection point, maximum slope of the line passing through the infl. point.
# Mar. 23, 2015
rm(list=ls())

library(msm)
library(reshape)
library(drc)
library(ggplot2)
library(ggthemes)
library(Rmisc)

getwd()
setwd("~/Documents/pRojects/exploratory_ideas/sigmoidal_cumulative_biomass_curves/fish-paste/")
# Load your dataset
#biom<-read.csv("ch_bio.csv",sep=";",dec=".",na.string="NA") # carico le biomasse
load('~/Documents/pRojects/exploratory_ideas/sigmoidal_cumulative_biomass_curves/fish-paste/data/analysis_ready_curve_fitting_data.Rdata')

oy4<-oy4[-which(oy4$ISLAND == "South Bank"),]
oy4<-oy4[-which(oy4$ISLAND == "Necker"),]
oy4<-oy4[-which(oy4$ISLAND == "Nihoa"),]
oy4<-oy4[-which(oy4$ISLAND == "Gardner"),]
#oy4<-oy4[-which(oy4$ISLAND == "Johnston"),]
#oy4<-oy4[-which(oy4$ISLAND == "Rota"),]



oy4<-droplevels(oy4)
oy4<-oy4[complete.cases(oy4),]
biom<-oy4
summary(biom)

### creating the data from which the simulated data will be drawn
n<-nrow(biom)
biom$TL_CV<-0.1
biom$Upper<-biom$TROPHICLEVEL + biom$TROPHICLEVEL * biom$TL_CV
biom$Lower<-biom$TROPHICLEVEL - biom$TROPHICLEVEL * biom$TL_CV

biom<-biom[-which(biom[,4] == 0),]
#biom$r<-biom$BIOMASS^2/biom$S_E^2 
#biom$l<-biom$S_E/biom$BIOMASS

BIOM<-biom

mc_results<-NULL
nreps<-1000
start<-proc.time()

for(jj in 1:nreps){
### creating a simulated dataset
biom<-BIOM
simdata<-biom[,c("SPECIES", "REGION", "ISLAND", "BIOMASS", "TROPHICLEVEL")]
n<-nrow(biom)

#simdata$BIOMASS_SIM<-rgamma(1:n, shape = biom$l, r=biom$r)
simdata$TROPHICLEVEL<- rtnorm(n, mean = biom$TROPHICLEVEL, sd = biom$TL_CV, lower = biom$Lower, upper = biom$Upper)


## match by col name and order
simdata<-simdata[,c("SPECIES", "TROPHICLEVEL", "ISLAND", "BIOMASS")]
names(simdata)<-c("SP", "TL", "Basin", "1981")

dataset<-simdata

#######################################################
#######################################################
#######################################################
# define some variables
#######################################################
#######################################################
#######################################################

# define the width of the step
step<-0.1

# starting and ending years
y0<-1981  # definisco l'anno di inizio
y1<-1981  #definisco l'anno finale

# Columns with biomasses
c0<-4   #numero della prima colonna contenente la biomassa
c1<-4  #numero dell'ultima colonna contenente la biomassa

# folder of output (optional, outputs go in the WD) load('~/Documents/pRojects/exploratory_ideas/sigmoidal_cumulative_biomass_curves/fish-paste/data/analysis_ready_curve_fitting_data.Rdata')
myout<-"out2"

# The script can works at the same time in different basins, in column 3 you can define the name of each basins
mybas<-3

#######################################################
#######################################################

# List of basins

filelist<-sub(",", "_",levels(biom[,mybas]),fixed=TRUE)
filelist<-sub(" ", "_",filelist,fixed=TRUE)
filelist<-sub(" ", "_",filelist,fixed=TRUE)
filelist<-sub("__", "_",filelist,fixed=TRUE)

# preparo un vettore con gli anni
years<-c(y0:y1)

# tabella che raccoglie i risultati
results<-expand.grid(year=years,bas=filelist,infe=NA,steep=NA,infl=NA,bio=NA)

# inizializzo un contatore
counter<-1

# The script loads the functions (from the other file)
source(file="scr/function_curve.r")

# Start a cycle
for(j in 1:length(filelist)){
  
  #j<-11
  
  temp<-subset(dataset,dataset[,mybas]==levels(dataset[,mybas])[j])
  temp<-subset(temp,!is.na(TL))  
  names(temp)[c(1,2)]<-c("SP","TL")
  biom<-temp
  rm(temp)
  
  
  mc<-mycast(biom,0.1,c0,c1,y0,y1)
   #png(filename=paste(myout,filelist[j],"_",step,".png",sep=""),width=1600,pointsize = 22,height=1200)
  
  nf<-layout(matrix(c(1:8),4,2,byrow=TRUE)) # Graphic parameters
  par(mar=c(0,0,0,0))
  
  x<-as.numeric(levels(mc$fTLr)[mc$fTLr])
  for (i in 2:ncol(mc)){
  i<-2	
    y<-cumsum(mc[,i])
    year<-years[i-1]
    mypar<-my.plot(x,y)
    results$year[counter]<-years[i-1]
    results$infe[counter]<-mypar[1]
    results$steep[counter]<-mypar[2]
    results$infl[counter]<-mypar[3]
    results$bio[counter]<-mypar[4]
    counter<-counter+1
  }
  #dev.off()
  
}

mc_results<-rbind(mc_results, results)
print("whoop whoop")
}
end<-proc.time()

end-start

write.csv(mc_results,file="mc_TL_only_bio_out.csv") # Output files
############################ end of curve fitting
### create summary data
tmp<-read.csv("mc_TL_only_bio_out.csv")

###### make up some summary plots - without CI for now
### add in the island level data from Julia
MC_RESULTS_SAVED<-tmp
mc_results<-tmp
#mc_results<-MC_RESULTS_SAVED
#########################

std <- function(x) sd(x)/sqrt(length(x))


theme_bw <- function(base_size = 12,  base_family = "") {
         # Starts with theme_grey and then modify some parts
          theme_grey(base_size = base_size, base_family = base_family) %+replace%
            theme(
              axis.text         = element_text(size = rel(0.8)),
              axis.ticks        = element_line(colour = "black"),
              axis.ticks.x      = element_blank(),
              axis.text.x       = element_text(angle = 90, hjust =1, vjust =0.4),
              legend.key        = element_rect(colour = "grey80"),
              legend.position   = "none",
              panel.background  = element_rect(fill = "white"),
              panel.border      = element_rect(fill = NA, colour = "grey50"),
              panel.grid.major  = element_line(colour = "grey90", size = 0.2),
              panel.grid.minor  = element_line(colour = "grey98", size = 0.5),
              strip.background  = element_rect(fill = "grey80", colour = "grey50"),
              strip.background  = element_rect(fill = "grey80", colour = "grey50")
            )
        }
       
######################       


names(mc_results)<-c("DUMMY1","DUMMY","ISLAND", "infe", "steep", "infl", "bio")

mc_results_sd<-aggregate(cbind(infl, bio, infe, steep) ~ ISLAND, data = mc_results,  std)

mc_results_mean<-aggregate(cbind(infl, bio, infe, steep) ~ ISLAND, data = mc_results,  mean)

is<-read.csv("data/Table of Island Information_v1_slope.csv")

tmp<-is[,c(1,4,7,8,9,10)]

names(tmp)[1]<-"ISLAND"

results<-merge(mc_results_mean, tmp, by = "ISLAND")

names(results)[6:10]<-c("ISLANDTYPE", "LANDAREA", "REEFAREA30","SLOPE_MEAN", "SLOPE_STDEV")
results<-merge(results, mc_results_sd, by = "ISLAND")

names(results)[c(2:5, 11:14)]<-c("infl", "bio", "infe", "steep", "infl_SE", "bio_SE", "infe_SE", "steep_SE") 


results$REGION<-results$ISLAND
levels(results$REGION)<-c(levels(results$REGION), levels(oy4$REGION))

results[results$ISLAND %in% c("Guam", "Rota", "Aguijan", "Tinian", "Saipan"),]$REGION<-"S.MARIAN"
results[results$ISLAND %in% c("Alamagan","Guguan","Sarigan","Pagan", "Agrihan", "Asuncion", "Maug", "Farallon_de_Pajaros", "AGS"),]$REGION<-"N.MARIAN"
results[results$ISLAND %in% c("Hawaii", "Kauai", "Lanai", "Maui", "Molokai", "Niihau", "Oahu"),]$REGION<-"MHI"
results[results$ISLAND %in% c("French_Frigate", "Kure", "Laysan", "Lisianski", "Maro", "Midway", "Pearl_&_Hermes"),]$REGION<-"NWHI"
results[results$ISLAND %in% c("Baker", "Howland", "Jarvis", "Johnston", "Kingman", "Palmyra", "Wake"),]$REGION<-"PRIAs"
results[results$ISLAND %in% c("Ofu_&_Olosega", "Rose", "Swains", "Tau","Tutuila"),]$REGION<-"SAMOA"
results<-droplevels(results)

#results<-results[,c(1, 2:11)]

levels(results$ISLANDTYPE)<-c(levels(results$ISLANDTYPE), "Island", "Atoll")
results[results$ISLANDTYPE %in% c("Guam", "Rota", "Aguijan", "Tinian", "Saipan"),]$REGION<-"S.MARIAN"

results[results$ISLANDTYPE %in% c("Closed atoll", "Open atoll"),]$ISLANDTYPE<-"Atoll"

############

p <- ggplot(results, aes(infl, bio, colour = REGION))
p + geom_point()+ scale_color_gdocs() + facet_grid(.~REGION) + geom_errorbarh(aes(xmax = infl + infl_SE, xmin = infl - infl_SE))+ geom_errorbar(aes(ymax = bio + bio_SE, ymin = bio - bio_SE)) +
ggtitle("Values of % cumulative biomass and trophic level inflection point")


############ steepness by island
cc<-c("#3d78d6","#a51d02","#efc033","#6aa74f","#674ea6","#3e84c4")

p <- ggplot(results, aes(ISLAND, steep))
p + geom_bar(stat="identity", aes(fill = REGION)) + scale_color_gdocs() + facet_grid(.~REGION, scales="free") + geom_errorbar(aes(ymax = steep + steep_SE, ymin = steep - steep_SE)) + scale_fill_manual(values = cc) +
ggtitle("Values of % cumulative biomass and trophic level inflection point") + theme_bw()

############ infl by island
cc<-c("#3d78d6","#a51d02","#efc033","#6aa74f","#674ea6","#3e84c4")

p <- ggplot(results, aes(ISLAND, infl))
p + geom_bar(stat="identity", aes(fill = REGION)) + scale_color_gdocs() + facet_grid(.~REGION, scales="free") + geom_errorbar(aes(ymax = infl + infl_SE, ymin = infl - infl_SE)) + scale_fill_manual(values = cc) +
ggtitle("Trophic level inflection point") + theme_bw()

############ cumbio by island
cc<-c("#3d78d6","#a51d02","#efc033","#6aa74f","#674ea6","#3e84c4")

p <- ggplot(results, aes(ISLAND,bio))
p + geom_bar(stat="identity", aes(fill = REGION)) + scale_color_gdocs() + facet_grid(.~REGION, scales="free") + geom_errorbar(aes(ymax = bio + bio_SE, ymin = bio - bio_SE)) + scale_fill_manual(values = cc) +
ggtitle("Proportion of cumulative biomass where inflection occurs") + theme_bw()

############ infe by island
cc<-c("#3d78d6","#a51d02","#efc033","#6aa74f","#674ea6","#3e84c4")

p <- ggplot(results, aes(ISLAND,infe))
p + geom_bar(stat="identity", aes(fill = REGION)) + scale_color_gdocs() + facet_grid(.~REGION, scales="free") + geom_errorbar(aes(ymax = infe + infe_SE, ymin = infe - infe_SE)) + scale_fill_manual(values = cc) +
ggtitle("Proportion of cumulative biomass where inflection occurs") + theme_bw()


############

p <- ggplot(results, aes(infl, bio, colour = REGION))
p + geom_point()+ scale_color_gdocs() + geom_errorbarh(aes(xmax = infl + infl_SE, xmin = infl - infl_SE))+ geom_errorbar(aes(ymax = bio + bio_SE, ymin = bio - bio_SE))+ 
ggtitle("Values of % cumulative biomass and trophic level inflection point")


###################

p <- ggplot(results, aes(infl, bio), colour = REGION)
p + geom_point(aes(colour = factor(REGION), shape = factor(ISLANDTYPE)))+ scale_color_gdocs() + geom_errorbarh(aes(xmax = infl + infl_SE, xmin = infl - infl_SE))+ geom_errorbar(aes(ymax = bio + bio_SE, ymin = bio - bio_SE))+ 
ggtitle("Values of % cumulative biomass and trophic level inflection point")




, scale_shape_manual(values=c(15,16,17,18)))) + scale_colour_manual(values=c("red", "blue", "black", "yellow", "red", "red"))

theme_calc()
 + ggtitle("Diamonds")
 + scale_color_calc())




