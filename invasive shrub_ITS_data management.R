library(phytools)
library(reshape2)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')

#mr bayes tree: 2 million runs, including native, invasive, and worldwide invasive species
ITStreeNexus<-read.nexus(file='ITS data\\MrBayes tree\\2 mil runs\\La Pierre_invasive shrub_ITS_consensus sequences_97similarity.nexus.con.tre')

#makes tree into a .phylo file instead of .newick or .tre
ITStree<-as.phylo(ITStreeNexus) 

#root tree with Mesorhizobium outgroup
ITStreeR<-root(ITStree, "USDA_3622_liaoningense", r=T)

#plot trees
plot.phylo(ITStree, use.edge.length=T)
plot.phylo(ITStree, use.edge.length=F)

#making a tree with only Bay Area strains (i.e., not worldwide strains)
BayAreaITStree<-drop.tip(ITStreeR, c('U115', 'UU22sfb'))
plot.phylo(BayAreaITStree, use.edge.length=F)

#making a tree with only Brady japonicum strains
BjITStree<-drop.tip(ITStreeR, c('ITS K01a', 'ITS V01s', 'ITS W01a', 'ITS W01n', 'ITS X01z', 'ITS E01h')) 
plot.phylo(BjITStree, use.edge.length=F)

#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreaITStree<-drop.tip(ITStreeR, c('U115', 'UU22sfb', 'ITS K01a', 'ITS V01s', 'ITS W01a', 'ITS W01n', 'ITS X01z', 'ITS E01h'))
plot.phylo(BjBayAreaITStree, use.edge.length=F)