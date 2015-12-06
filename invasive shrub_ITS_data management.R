library(phytools)
library(reshape2)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')


#####################
# tree management
#####################

#mr bayes tree: 2 million runs, including native, invasive, and worldwide invasive species
ITStreeNexus<-read.nexus(file='ITS data\\MrBayes tree\\2 mil runs\\La Pierre_invasive shrub_ITS_consensus sequences_97similarity.nexus.con.tre')

#makes tree into a .phylo file instead of .newick or .tre
ITStree<-as.phylo(ITStreeNexus) 

#root tree with Mesorhizobium outgroup
ITStreeR<-root(ITStree, "USDA_3622_liaoningense", r=T)

#plot trees
plot.phylo(ITStreeR, use.edge.length=T)
plot.phylo(ITStreeR, use.edge.length=F)

# #making a tree with only Bay Area strains (i.e., not worldwide strains)
# BayAreaITStree<-drop.tip(ITStreeR, c('U115', 'UU22sfb'))
# plot.phylo(BayAreaITStree, use.edge.length=F)
# 
# #making a tree with only Brady japonicum strains
# BjITStree<-drop.tip(ITStreeR, c('ITS_K01a', 'IITS_V01s', 'ITS_W01a', 'ITS_W01n', 'ITS_X01z', 'ITS_E01h')) 
# plot.phylo(BjITStree, use.edge.length=F)

#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreaITStree<-drop.tip(ITStreeR, c('U115', 'UU22sfb', 'ITS_K01a', 'ITS_V01s', 'ITS_W01a', 'ITS_W01n', 'ITS_X01z', 'ITS_E01h'))
plot.phylo(BjBayAreaITStree, use.edge.length=F)


#####################
# strain data management
#####################
#read in nodule data
ITSnodules <- read.csv('strain data\\La Pierre_invasion molecular manuscript_strain information_092515.csv')%>%
  select(plant_species, plant_status, nodule_ID, ITS_contig_97sim)%>%
  #remove strains not in tree
  filter(ITS_contig_97sim!='ITS_23_5N', ITS_contig_97sim!='ITS_29_2N')
  
#create an interaction matrix of strains for each plant species
ITSinteractionMatrix <- ITSnodules%>%
  select(plant_species, plant_status, nodule_ID, ITS_contig_97sim)%>%
  filter(ITS_contig_97sim!='')%>%
  mutate(interaction=1)%>%
  spread(key=ITS_contig_97sim, value=interaction, fill=0)
  
# #subset out only Bay Area strains
# ITSBayAreaInteractionMatrix <- ITSinteractionMatrix%>%
#   filter(plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
#   select(-U115, -UU22sfb)%>%
#   #get summary interaction matrix (sum of interactions by species)
#   gather(key=ITS_contig_97sim, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
#   group_by(plant_status, plant_species, ITS_contig_97sim)%>%
#   summarise(interaction=sum(interaction))%>%
#   spread(key=ITS_contig_97sim, value=interaction)
#   
# #subset out only brady japonicum strains
# ITSbjInteractionMatrix <- ITSinteractionMatrix%>%
#   select(-ITS_K01a, -ITS_V01s, -ITS_W01a, -ITS_W01n, -ITS_X01z, -ITS_E01h)%>%
#   filter(nodule_ID!='K01a', nodule_ID!='V01s', nodule_ID!='W01a', nodule_ID!='W01n', nodule_ID!='ITS_X01z', nodule_ID!='E01h')%>%
#   #get summary interaction matrix (sum of interactions by species)
#   gather(key=ITS_contig_97sim, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
#   group_by(plant_status, plant_species, ITS_contig_97sim)%>%
#   summarise(interaction=sum(interaction))%>%
#   spread(key=ITS_contig_97sim, value=interaction)

#subset out only Bay Area and brady japonicum strains
ITSbjBayAreaInteractionMatrix <- ITSinteractionMatrix%>%
  select(-ITS_K01a, -ITS_V01s, -ITS_W01a, -ITS_W01n, -ITS_X01z, -ITS_E01h, -U115, -UU22sfb)%>%
  filter(nodule_ID!='K01a', nodule_ID!='V01s', nodule_ID!='W01a', nodule_ID!='W01n', nodule_ID!='ITS_X01z', nodule_ID!='E01h', plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
  #get summary interaction matrix (sum of interactions by species)
  gather(key=ITS_contig_97sim, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
  group_by(plant_status, plant_species, ITS_contig_97sim)%>%
  summarise(interaction=sum(interaction))%>%
  spread(key=ITS_contig_97sim, value=interaction)














