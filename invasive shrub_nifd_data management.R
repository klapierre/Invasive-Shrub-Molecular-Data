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

# #making a tree with only Bay Area strains (i.e., not worldwide strains)
# BayAreanifdtree<-drop.tip(nifdtreeR, c('nifd_03'))
# plot.phylo(BayAreanifdtree, use.edge.length=F)

#not relevent because all strains had brady japonicum
# #making a tree with only Brady japonicum strains
# Bjnifdtree<-drop.tip(nifdtreeR, c('ITS K01a', 'ITS V01s', 'ITS W01a', 'ITS W01n', 'ITS X01z', 'ITS E01h')) # plot.phylo(Bjnifdtree, use.edge.length=F)

#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreanifdtree<-drop.tip(nifdtreeR, c('nifd_03'))
plot.phylo(BjBayAreanifdtree, use.edge.length=F)


#####################
# strain data management
#####################
#read in nodule data
nifdNodules <- read.csv('strain data\\La Pierre_invasion molecular manuscript_strain information_092515.csv')%>%
  select(plant_species, plant_status, nodule_ID, nifd_contig_95sim)

#create an interaction matrix of strains for each plant species
nifdInteractionMatrix <- nifdNodules%>%
  select(plant_species, plant_status, nodule_ID, nifd_contig_95sim)%>%
  filter(nifd_contig_95sim!='')%>%
  mutate(interaction=1)%>%
  spread(key=nifd_contig_95sim, value=interaction, fill=0)

# #subset out only Bay Area strains
# nifdBayAreaInteractionMatrix <- nifdInteractionMatrix%>%
#   filter(plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
#   #get summary interaction matrix (sum of interactions by species)
#   gather(key=nifd_contig_95sim, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
#   group_by(plant_status, plant_species, nifd_contig_95sim)%>%
#   summarise(interaction=sum(interaction))%>%
#   spread(key=nifd_contig_95sim, value=interaction)

# #not relevent, because all non-brady japonicum strains grouped with a strain that was a bj
# #subset out only brady japonicum strains
# ITSbjInteractionMatrix <- nifdInteractionMatrix%>%
#   #get summary interaction matrix (sum of interactions by species)
#   gather(key=nifd_contig_95sim, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
#   group_by(plant_status, plant_species, nifd_contig_95sim)%>%
#   summarise(interaction=sum(interaction))%>%
#   spread(key=nifd_contig_95sim, value=interaction)

#subset out only Bay Area and brady japonicum strains
nifdBjBayAreaInteractionMatrix <- nifdInteractionMatrix%>%
  filter(plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
  #get summary interaction matrix (sum of interactions by species)
  gather(key=ITS_contig_97sim, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
  group_by(plant_status, plant_species, ITS_contig_97sim)%>%
  summarise(interaction=sum(interaction))%>%
  spread(key=ITS_contig_97sim, value=interaction)
