
library(pegas)
library(ape)
library(spider)

setwd("/Users/clareadams/Documents/mtDNA-eDNA-project/eDNA-part/data/SpeciesAccumulationCurve")

PauaTissueHaps<-read.FASTA("Align71Seqs-realAmplicon.fasta") #read in FASTA files
haplotype(PauaTissueHaps) #counts number of haplotypes

x <- haploAccum(PauaTissueHaps, method = "random", permutations = 1000) #makes haploaccumulation data taking your file

par(mar=c(5,6,4,1)+.1)
plot(x,xlab="Sequences", main=NULL, cex.lab=3, cex.axis=2)   #plots your data

# three or four is low but plausible.