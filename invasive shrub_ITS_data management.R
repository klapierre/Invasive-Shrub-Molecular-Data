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
ITStreeR<-root(ITStree, "ITS_Mesorhizobium", r=T)

#plot trees
plot.phylo(ITStree, use.edge.length=T)
plot.phylo(ITStree, use.edge.length=F)

#making a tree with only Bay Area strains (i.e., not worldwide strains)
BayAreaITStree<-drop.tip(ITStreeR, c("GV101", "GV104", "GV115", "GV102", "GV144", "GV159", "U118", "U115", "UU22sfb", "U29", "U21", "U23", "U13", "U214", "U12", "GV135", "GV137"))
plot.phylo(BayAreaITStree, use.edge.length=F)

#making a tree with only Brady japonicum strains
BjITStree<-drop.tip(ITStreeR, c("Bsp01", "Rl01", "Rl02", "Rl03", "Rl04", "Rl05", "Rl06", "Rl07", "GV137")) 
plot.phylo(BjITStree, use.edge.length=F)

#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreaITStree<-drop.tip(ITStreeR, c("Bsp01", "Rl01", "Rl02", "Rl03", "Rl04", "Rl05", "Rl06", "Rl07", "GV101", "GV104", "GV115", "GV102", "GV144", "GV159", "U118", "U115", "UU22sfb", "U29", "U21", "U23", "U13", "U214", "U12", "GV135", "GV137"))
plot.phylo(BjBayAreaITStree, use.edge.length=F)