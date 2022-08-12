## eDNA mtDNA amplicons

# first R analysis

#setwd("~/Documents/PhD/2nd-Project-emtDNA/data/28-12-2020")
setwd("/Users/clareadams/Documents/mtDNA-eDNA-project/eDNA-part/data/eDNA-mtDNA-amplicon-results")

#eDNA.amp<-read.csv("Paua.eDNA.1.csv")   # ee=1.00
#or
eDNA.amp<-read.csv("Paua.eDNA.ee10.csv")  # ee=0.10

# ee-0.1 is more restricitve (you get ~75%) of ee1.0 amplicons
# and it makes almost no difference in the final ratio outcome
# For this paper, we used ee=0.1 because we wanted strict filtering
# However, it is what we found to be effective enough for haplotype stuff

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
     "Warrington","Ohau")


head(Long.eDNA)


## a generic plot of

p10<-ggplot(Long.eDNA, aes(x=OTU.ID, y=value, color=OTU.ID))+facet_wrap(Sample ~ .) + geom_boxplot()+theme_bw()+ggtitle("ZOTUs-ee0.10")
p10

library(gridExtra)
grid.arrange(p50,p10, nrow=1, ncol=2)

# Basically the same graph, but you're throwing out about 1/4th of the sequences with ee0.1

### Making a pie chart
library(dplyr)

total<-aggregate(value~OTU.ID+Location, Long.eDNA, sum)
totalW<-filter(total, Location == "Warrington")
totalO<-filter(total, Location == "Ohau")

library(ggplot2)
#create pie chart
PieW<-ggplot(totalW, aes(x="", y=value, fill=OTU.ID)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void()+ geom_text(aes(label = paste0(value)), position = position_stack(vjust=0.5))+ ggtitle("eDNA Haplotype read counts\nfor Warrington, Otago")+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1", "Zotu2", "Zotu3", "Zotu4"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4"))

PieO<-ggplot(totalO, aes(x="", y=value, fill=OTU.ID)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void()+ geom_text(aes(label = paste0(value)), position = position_stack(vjust=0.5))+ ggtitle("eDNA Haplotype read counts\nfor Ohau Point, Kaikoura")+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1", "Zotu2", "Zotu3", "Zotu4"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4"))


#library(gridExtra)
grid.arrange(PieW,PieO, nrow=1, ncol=2)


### Add in Expected Tissue Haplotype counts
Tissue<-read.csv("TissueHaps.csv")
head(Tissue)

PieT<-ggplot(Tissue, aes(x="", y=Count, fill=Haplotype)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void()+ geom_text(aes(label = paste0(Count)), position = position_stack(vjust=0.5))+ ggtitle("Individuals' tissue haplotypes\nfor Warrington, Otago\nNumber of individuals")+scale_fill_discrete(name="Haplotypes")
PieT

grid.arrange(PieT,PieW,PieO, nrow=1, ncol=3)

#############################
#### Stacked Bar Charts #####
#############################

## Just eDNA haplotype ratios ONLY

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
withr::with_options(
  list(ggplot2.discrete.fill = okabe),
  print(ggplot(Long.eDNA, aes(fill=OTU.ID, y=value, x=Location)) +
          geom_bar(position="fill", stat="identity")+theme_bw()+ggtitle("eDNA haplotype ratios")+ scale_fill_discrete(name = "Haplotype", labels = c("Haplotype 1", "Haplotype 2", "Haplotype 3","Haplotype 4"))+ylab("Percentage of total sampled")+
          theme(axis.text=element_text(size=14),
                axis.title=element_text(size=14,face="bold"),
                legend.text=element_text(size=12),
                legend.title = element_text(size=14),
                plot.title = element_text(size=40,face = "bold"))))


ggplot(Long.eDNA, aes(fill=OTU.ID, y=value, x=Location)) +
  geom_bar(position="fill", stat="identity")+theme_bw()+ggtitle("eDNA Haplotype ratios")+ scale_fill_discrete(name = "Haplotype", labels = c("Haplotype 1", "Haplotype 2", "Haplotype 3","Haplotype 4"))+ylab("Percentage of total sampled")


fit <- oneway.test(sqrt(value) ~ OTU.ID+Location, data=Long.eDNA)
summary(fit)
anova(fit)
hist(residuals(fit),col="darkgray")
plot(fitted(fit),residuals(fit))
leveneTest(Long.eDNA$value,Long.eDNA$Location, center=median)
median()kruskal.test(value ~ OTU.ID, data = Long.eDNA)
kruskal.test(value ~ Location, data = Long.eDNA)

install.packages("ggpubr")
library("ggpubr")
ggline(Long.eDNA, x = "Location", y = "value", color = "OTU.ID",
       add = c("mean_se", "dotplot"))

TukeyHSD(fit, which = "OTU.ID")


## With Tissue!

library(dplyr)

TissueW<-read.csv("WarHaps.csv")

Combined=Long.eDNA %>% bind_rows(TissueW)


# plot all three -- warrington tissue, warrington eDNA and Ohau eDNA
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
withr::with_options(
  list(ggplot2.discrete.fill = okabe),
  print(ggplot(Combined, aes(fill=OTU.ID, y=value, x=Location)) +
          geom_bar(position="fill", stat="identity")+theme_bw()+ggtitle("Tissue and environmental DNA Haplotype ratios")+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1","Zotu2", "Zotu3", "Zotu4", "Zotu5", "Zotu6"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+scale_x_discrete(limits = rev(levels(as.factor(Combined$Location))))+ylab("Percentage of total sampled")+
          theme(axis.text=element_text(size=14),
                axis.title=element_text(size=14,face="bold"),
                legend.text=element_text(size=12),
                legend.title = element_text(size=14),
                plot.title = element_text(size=40,face = "bold")))
)


ggplot(Combined, aes(fill=OTU.ID, y=value, x=Location)) +
  geom_bar(position="fill", stat="identity")+theme_bw()+ggtitle("Tissue and environmental DNA Haplotype ratios")+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1","Zotu2", "Zotu3", "Zotu4", "Zotu5", "Zotu6"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+scale_x_discrete(limits = rev(levels(as.factor(Combined$Location))))+ylab("Percentage of total sampled")

# stack instead of fill does not work well because there are few invidivuals compared to the amount of sequences we have


# plot by sample
library(dplyr)
# Filter for only Warrington samples
totalW<-filter(Long.eDNA, Location == "Warrington")

# Re-level
totalW$Sample<- factor(totalW$Sample, levels=c("Sample_1","Sample_2","Sample_3","Sample_4","Sample_5","Sample_6","Sample_8","Sample_9","Sample_10"))

ggplot(totalW, aes(fill=OTU.ID, y=value, x=Sample)) +
  geom_bar(position="stack", stat="identity")+theme_bw()+ggtitle("Environmental DNA haplotypes per sample")+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1", "Zotu2", "Zotu3", "Zotu4", "Zotu5", "Zotu6"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+ylab("Number of Reads")+scale_fill_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#000000"))

# each eDNA haplotypes per sample for Warrington
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
withr::with_options(
  list(ggplot2.discrete.fill = okabe),
  print(ggplot(totalW, aes(fill=OTU.ID, y=value, x=Sample)) +
          geom_bar(position="stack", stat="identity")+theme_bw()+ggtitle("eDNA haplotypes per sample")+scale_fill_discrete(name="Haplotypes", breaks=c("Zotu1", "Zotu2", "Zotu3", "Zotu4", "Zotu5", "Zotu6"),labels=c("Haplotype 1","Haplotype 2","Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+ylab("Number of Reads")+
          theme(axis.text=element_text(size=14),
                axis.title=element_text(size=14,face="bold"),
                legend.text=element_text(size=12),
                legend.title = element_text(size=14),
                plot.title = element_text(size=40,face = "bold"))))




# Warrington eDNA vs Warrington Tissue

totalWC<-filter(Combined, Combined$Location == "Warrington" | Combined$Location == "Warrington - Tissue")

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
withr::with_options(
  list(ggplot2.discrete.fill = okabe),
  print(ggplot(totalWC, aes(fill=OTU.ID, y=value, x=Location)) +
          geom_bar(position="fill", stat="identity")+theme_bw()+ggtitle("Warrington eDNA vs Tissue")+ scale_fill_discrete(name = "Haplotype", labels = c("Haplotype 1", "Haplotype 2", "Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+ylab("Percentage of total sampled")+
          theme(axis.text=element_text(size=14),
                axis.title=element_text(size=14,face="bold"),
                legend.text=element_text(size=12),
                legend.title = element_text(size=14),
                plot.title = element_text(size=40,face = "bold")) ))

ggplot(totalWC, aes(fill=OTU.ID, y=value, x=Location)) +
  geom_bar(position="fill", stat="identity")+theme_bw()+ggtitle("Warrington eDNA vs Tissue")+ scale_fill_discrete(name = "Haplotype", labels = c("Haplotype 1", "Haplotype 2", "Haplotype 3","Haplotype 4","Haplotype 5","Haplotype 6"))+ylab("Percentage of total sampled")


#### Looking at rank order abundances

library(RADanalysis)

SummedValues<-aggregate(value ~ OTU.ID + Location, data = Combined, FUN = sum)

#find the totals
hyp_grouped <- aggregate(SummedValues$value, by=list(Category=SummedValues$Location), FUN=sum)
# for some reason these need to be highlighted and run in chunks to work
# add in the total values onto the dataframe
dataframe1 <- SummedValues %>% rowwise %>%
  do({
    result = as_tibble(.)
    result$Total = hyp_grouped[hyp_grouped$Category == result$Location, 2]
    result
  })
# find the percentages
dataframe2 <- dataframe1 %>%
  rowwise %>%
  do({
    result = as_tibble(.)
    result$Percentage = (result$value/(result$Total))*100
    result
  })
# now you have a percentage'd dataframe
head(dataframe2)

### Okay, now you can create a rank-order abundance

# sort the data

Ordered <- dataframe2[order(-dataframe2$Percentage),]
head(Ordered)

# remove white space
Ordered$Location<-gsub(" ", "", Ordered$Location, fixed = TRUE)
Ordered$Location<-gsub("-", "", Ordered$Location, fixed = TRUE)
# split data by location
SplitData<-split(Ordered, Ordered$Location)

#plot out each location

plot(SplitData$WarringtonTissue$Percentage,type="o",xlim=c(1,6),lwd=2,col="green", pch=1,cex=2,xlab="Haplotype Rank",ylab="Proportion of Abundance", main="Whittaker rank-abundance plot of haplotypes at Warrington and Ohau")
lines(SplitData$Warrington$Percentage,type="o",xlim=c(1,6),lwd=2,col="blue", cex=2,pch=2)
lines(SplitData$Ohau$Percentage,type="o",xlim=c(1,6),lwd=2,col="red", cex=2,pch=3)
legend(5,45,legend=c("Warrington Tissue", "Warrington eDNA", "Ohau eDNA"),col=c("green","blue","red"), lty=1:1,pch=1:3)
