## eDNA mtDNA amplicons
## Final R Script
## January 25th 2021
pal9 <- c("#920000","#db6d00","#24ff24","#009292","#004949","#b6dbff","#006ddb","#b66dff","#490092")

#setwd("~/Documents/PhD/2nd-Project-emtDNA/data/28-12-2020")
#setwd("/Users/clareadams/Documents/mtDNA-eDNA-project/eDNA-part/data/eDNA-mtDNA-amplicon-results")

setwd("/Volumes/CLARE_PHD/Documents/mtDNA-eDNA-project/eDNA-part/data/eDNA-mtDNA-amplicon-results")

eDNA.amp<-read.csv("Paua.eDNA.ee10.csv")  # ee=0.10


library(tidyverse)
library(reshape2)
library(ggplot2)


## Longform the data

Long.eDNA<-melt(eDNA.amp, id.vars=c("OTU.ID"))

## Get rid of the Negatives and Nothings
Long.eDNA<-Long.eDNA[!(Long.eDNA$variable=="Nothing5"),]
Long.eDNA<-Long.eDNA[!(Long.eDNA$variable=="Nothing6"),]

## Add in something for the samples

Long.eDNA$Sample <- ifelse(grepl("Sample_10_", Long.eDNA$variable, ignore.case = T),
                          "Sample_10",ifelse(grepl("Sample_1_", Long.eDNA$variable, ignore.case = T),
                          "Sample_1",ifelse(grepl("Sample_2_", Long.eDNA$variable, ignore.case = T),
                          "Sample_2",ifelse(grepl("Sample_3_", Long.eDNA$variable, ignore.case = T),
                          "Sample_3",ifelse(grepl("Sample_4_", Long.eDNA$variable, ignore.case = T),
                          "Sample_4",ifelse(grepl("Sample_5_", Long.eDNA$variable, ignore.case = T),
                          "Sample_5",ifelse(grepl("Sample_6_", Long.eDNA$variable, ignore.case = T),
                          "Sample_6",ifelse(grepl("Sample_7_", Long.eDNA$variable, ignore.case = T),
                          "Sample_7",ifelse(grepl("Sample_8_", Long.eDNA$variable, ignore.case = T),
                          "Sample_8",ifelse(grepl("Sample_9_", Long.eDNA$variable, ignore.case = T),
                          "Sample_9","Ohau"))))))))))

Long.eDNA$Location <- ifelse(grepl("Sample_", Long.eDNA$variable, ignore.case = T),
                             "Warrington","Ohau Point")

head(Long.eDNA)

library(dplyr)
total<-aggregate(value~OTU.ID+Location, Long.eDNA, sum)
totalW<-filter(total, Location == "Warrington")
totalO<-filter(total, Location == "Ohau")


library(dplyr)
TissueW<-read.csv("WarHaps.csv")
Combined=Long.eDNA %>% bind_rows(TissueW)

# plot all three -- warrington tissue, warrington eDNA and Ohau eDNA
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pal9 <- c("#920000","#db6d00","#009292","#004949","#b6dbff","#006ddb","#b66dff","#490092")

##"#24ff24",    +ggtitle("Tissue and environmental DNA Haplotype ratios")

withr::with_options(
  list(ggplot2.discrete.fill = pal9),
  print(ggplot(Combined, aes(fill=OTU.ID, y=value, x=Location)) +
          geom_bar(position="fill", stat="identity")+theme_bw()+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1","Zotu2", "Zotu3", "Zotu4", "Zotu5", "Zotu6"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+scale_x_discrete(limits = rev(levels(as.factor(Combined$Location))))+ylab("Percentage of total sampled")+
          theme(axis.text=element_text(size=20),
                axis.title=element_text(size=20,face="bold"),
                legend.text=element_text(size=16),
                legend.title = element_text(size=20))))

# plot by sample

totalW<-filter(Long.eDNA, Location == "Warrington")

# Re-level
totalW$Sample<- factor(totalW$Sample, levels=c("Sample_1","Sample_2","Sample_3","Sample_4","Sample_5","Sample_6","Sample_8","Sample_9","Sample_10"))

#+ggtitle("eDNA haplotypes per sample"

# each eDNA haplotypes per sample for Warrington
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
withr::with_options(
  list(ggplot2.discrete.fill = pal9),
  print(ggplot(totalW, aes(fill=OTU.ID, y=value, x=Sample)) +
          geom_bar(position="stack", stat="identity")+theme_bw())+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1", "Zotu2", "Zotu3", "Zotu4", "Zotu5", "Zotu6"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+ylab("Number of Reads")+
          theme(axis.text=element_text(size=15),
                axis.title=element_text(size=25,face="bold"),
                legend.text=element_text(size=16),
                legend.title = element_text(size=20,face="bold")))

Long.eDNA$Sample = factor(Long.eDNA$Sample, levels = c('Sample_1', 'Sample_2','Sample_3','Sample_4','Sample_5','Sample_6','Sample_8','Sample_9','Sample_10', 'Ohau'))
withr::with_options(
  list(ggplot2.discrete.fill = pal9),
  print(ggplot(Long.eDNA, aes(fill=OTU.ID, y=value, x=Sample)) +
          geom_bar(position="stack", stat="identity")+theme_bw())+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1", "Zotu2", "Zotu3", "Zotu4", "Zotu5", "Zotu6"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+ylab("Number of Reads")+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=25,face="bold"),
          legend.text=element_text(size=16),
          legend.title = element_text(size=20,face="bold")))



# Warrington eDNA vs Warrington Tissue

totalWC<-filter(Combined, Combined$Location == "Warrington" | Combined$Location == "Warrington - Tissue")
totalWC$Location[totalWC$Location == "Warrington"] <- "Warrington - eDNA"
totalWC$Location <- factor(totalWC$Location, levels = c("Warrington - Tissue", "Warrington - eDNA"))


okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
withr::with_options(
  list(ggplot2.discrete.fill = pal9),
  print(ggplot(totalWC, aes(fill=OTU.ID, y=value, x=Location)) +
          geom_bar(position="fill", stat="identity")+theme_bw()+scale_fill_discrete(name = "Haplotype", labels = c("Haplotype 1", "Haplotype 2", "Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+ylab("Proportion of total sampled")+
          theme(axis.text=element_text(size=20),
                axis.title=element_text(size=25,face="bold"),
                legend.text=element_text(size=16),
                legend.title = element_text(size=20,face="bold"))))

## eDNA - Ohau vs Warrington

str(Long.eDNA)

Long.eDNA$Location <- factor(Long.eDNA$Location, levels = c("Warrington", "Ohau Point"))

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
withr::with_options(
  list(ggplot2.discrete.fill = pal9),
  print(ggplot(Long.eDNA, aes(fill=OTU.ID, y=value, x=Location)) +
          geom_bar(position="fill", stat="identity")+theme_bw()+ scale_fill_discrete(name = "Haplotype", labels = c("Haplotype 1", "Haplotype 2", "Haplotype 3","Haplotype 4"))+ylab("Proportion of total sampled")+
          theme(axis.text=element_text(size=20),
                axis.title=element_text(size=25,face="bold"),
                legend.text=element_text(size=16),
                legend.title = element_text(size=20,face="bold"))))


## Haplotype Accumulation plot


library(pegas)
library(ape)
library(spider)

setwd("/Users/clareadams/Documents/mtDNA-eDNA-project/eDNA-part/data/SpeciesAccumulationCurve")

PauaTissueHaps<-read.FASTA("Align71Seqs-realAmplicon.fasta") #read in FASTA files
haplotype(PauaTissueHaps) #counts number of haplotypes

x <- haploAccum(PauaTissueHaps, method = "random", permutations = 1000) #makes haploaccumulation data taking your file

plot(x) #plots your data

# three or four is low but plausible.




################### graph of map #########################

## doing a lil R project to make a map of my sampling site!

library(maps)
library(mapdata)
library(mapproj)
library(maptools)

map('nzHires', xlim=c(165,179), ylim=c(-50,-35))
nz<-map('nzHires')
visual <- ggplotly(nz, height = 1.2 * 600, width = 600, tooltip=c("text"),
                   hoverinfo='hide', dynamicTicks = FALSE)

samps <- read.csv("SamplingSite.csv", header=TRUE, stringsAsFactors=T)

#plot all city names of NZ onto
#map.cities(country="New Zealand", label=TRUE, cex=1,xlim=c(165,179), ylim=c(-50,-35), pch=20)

map('nzHires', xlim=c(164,179), ylim=c(-50,-35))
points(samps$Long, samps$Lat, pch=19, col="red",cex=1)
pointLabel(samps$Long, samps$Lat, samps$Sampling_site)



###
##alternatively, ggmap

library(ggmap)

map.nz <- map_data(map = "nz")
p1 <- ggplot(map.nz, aes(x = long, y = lat, group=group))
p1 <- p1 + geom_polygon()
p1 <- p1 + labs(title = "New Zealand")
p1
setwd("~/Documents/mtDNA-eDNA-project/eDNA-part/data/eDNA-mtDNA-amplicon-results")
samps<- read.csv("MapLatLong.csv", header=TRUE, stringsAsFactors = T)

realsamps<-samps[1:2,]

p2=p1+geom_point(data = realsamps, aes(x = Longitude, y = Latitude,shape = Sampling_Site, colour = Sampling_Site), size = 7,inherit.aes = FALSE)+theme_bw()

p2

pal9<-c("#920000", "#009292","#000000")

p2+ geom_text(data = realsamps, aes(x = Longitude, y = Latitude, label = Sampling_Site), hjust = -0.2, colour = "black", size = 5,inherit.aes = FALSE)+scale_color_manual(values=pal9)+xlab("Longitude")+ylab("Latitude")+theme(axis.title=element_text(size=25),axis.text = element_text(size=15),legend.title=element_text(size=15),legend.text = element_text(size=10))+scale_fill_discrete("Sampling Site")

p2+scale_color_manual(values=pal9)+xlab("Longitude")+ylab("Latitude")+theme(axis.title=element_text(size=15),legend.title=element_text(size=15))
