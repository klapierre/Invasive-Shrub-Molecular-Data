library(phytools)
library(reshape2)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')

#mr bayes tree: 2 million runs, including native, invasive, and worldwide invasive species
nifdtreeNexus<-read.nexus(file='nifd data\\mr bayes tree_01272016\\nifd_consensus_outgroups_01252016.nexus.con.tre')

#makes tree into a .phylo file instead of .newick or .tre
nifdtree<-as.phylo(nifdtreeNexus) 

#root tree with Mesorhizobium outgroup
nifdtreeR<-root(nifdtree, "Mesorhizobium_ciceri_USDA3383_nifd_GQ167280", r=T)

#plot trees
plot.phylo(nifdtreeR, use.edge.length=T)
plot.phylo(nifdtreeR, use.edge.length=F)

# #making a tree with only Bay Area strains (i.e., not worldwide strains)
# BayAreanifdtree<-drop.tip(nifdtreeR, c('nifd_03'))
# plot.phylo(BayAreanifdtree, use.edge.length=F)

#remove strains that could not be concatenated with ITS sequences
BjBayAreanifdtree<-drop.tip(nifdtreeR, c('nifd_028')) # plot.phylo(Bjnifdtree, use.edge.length=F)
plot.phylo(BjBayAreanifdtree, use.edge.length=F)


#making a tree with only Brady japonicum strains, excluding Australian Ulex
Bjnifdtree<-drop.tip(nifdtreeR, c('nifd_028', 'nifd_019', 'nifd_021', 'nifd_029', 'Rhizobium_leguminosarum_X01z_nifd'))
plot.phylo(Bjnifdtree, show.node.label=T)
nodelabels(frame='none', cex=0.8)


#####################
# strain data management
#####################
#read in nodule data
nifdNodules <- read.csv('strain data\\La Pierre_invasion molecular manuscript_strain information_01252016.csv')%>%
  select(plant_species, plant_status, nodule_ID, nifd_OTU_99, concatenated_OTU_98)

#create an interaction matrix of strains for each plant species
nifdInteractionMatrix <- nifdNodules%>%
  select(plant_species, plant_status, nodule_ID, nifd_OTU_99, concatenated_OTU_98)%>%
  filter(nifd_OTU_99!='', nifd_OTU_99!='nifd_028', concatenated_OTU_98!='')%>%
  select(-concatenated_OTU_98)%>%
  mutate(interaction=1)%>%
  spread(key=nifd_OTU_99, value=interaction, fill=0)

#subset out only Bay Area and brady japonicum strains
nifdBjBayAreaInteractionMatrix <- nifdInteractionMatrix%>%
  select(-nifd_019, -nifd_021, -nifd_029)%>%
  filter(plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
  #get summary interaction matrix (sum of interactions by species)
  gather(key=nifd_OTU_99, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
  group_by(plant_status, plant_species, nifd_OTU_99)%>%
  summarise(interaction=sum(interaction))%>%
  spread(key=nifd_OTU_99, value=interaction)
