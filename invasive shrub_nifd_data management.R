library(phytools)
library(reshape2)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')

#mr bayes tree: 2 million runs, including native, invasive, and worldwide invasive species
nifdtreeNexus<-read.nexus(file='nifd data\\mr bayes tree\\2 mil runs\\La Pierre_invasive shrub_nifd_consensus_95sim.nexus.con.tre')

#makes tree into a .phylo file instead of .newick or .tre
nifdtree<-as.phylo(nifdtreeNexus) 

#root tree with Mesorhizobium outgroup
nifdtreeR<-root(nifdtree, "USDA_3622_liaoningense", r=T)

#plot trees
plot.phylo(nifdtreeR, use.edge.length=T)
plot.phylo(nifdtreeR, use.edge.length=F)

#making a tree with only Bay Area strains (i.e., not worldwide strains)
BayAreanifdtree<-drop.tip(nifdtreeR, c('nifd_03'))
plot.phylo(BayAreanifdtree, use.edge.length=F)

#not relevent because all strains had brady japonicum
# #making a tree with only Brady japonicum strains
# Bjnifdtree<-drop.tip(nifdtreeR, c('ITS K01a', 'ITS V01s', 'ITS W01a', 'ITS W01n', 'ITS X01z', 'ITS E01h')) # plot.phylo(Bjnifdtree, use.edge.length=F)

#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreanifdtree<-drop.tip(nifdtreeR, c('nifd_03'))
plot.phylo(BjBayAreanifdtree, use.edge.length=F)