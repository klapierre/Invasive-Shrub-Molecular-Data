library(phytools)
library(reshape2)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')


#####################
# tree management
#####################

#mr bayes tree: 2 million runs, including native, invasive, and worldwide invasive species
conctreeNexus<-read.nexus(file='concatenated data\\MrBayes tree_01272016\\concatenated_consensus_outgroups_01272016.nexus.con.tre')

#makes tree into a .phylo file instead of .newick or .tre
contree<-as.phylo(conctreeNexus) 

#root tree with Mesorhizobium outgroup
contreeR<-root(contree, "Mesorhizobium_ciceri_USDA3383_ITS_AF345262_1", r=T)

#plot trees
plot.phylo(contreeR, use.edge.length=T)
plot.phylo(contreeR, use.edge.length=F)


#making a tree with only Brady japonicum strains and not including Australian Ulex
Bjconctree<-drop.tip(contreeR, c('conc_019', 'conc_011', 'Rhizobium_leguminosarum_X01z'))
plot.phylo(Bjconctree, use.edge.length=F)
nodelabels(frame='none', cex=0.8)


#making a tree with only Brady japonicum strains and not including Australian Ulex and dropping conc_020 because such a long branch
Bjconctree<-drop.tip(contreeR, c('conc_019', 'conc_011', 'Rhizobium_leguminosarum_X01z', 'conc_020'))
plot.phylo(Bjconctree, use.edge.length=T)
nodelabels(frame='none', cex=0.8)


#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreaconctree<-drop.tip(Bjconctree, c('conc_011', 'conc_021', 'conc_022')) 
plot.phylo(BjBayAreaconctree, use.edge.length=F)



#####################
# strain data management
#####################
#read in nodule data
concnodules <- read.csv('strain data\\La Pierre_invasion molecular manuscript_strain information_01252016.csv')%>%
  select(plant_species, plant_status, nodule_ID, concatenated_OTU_98)

#create an interaction matrix of strains for each plant species
concinteractionMatrix <- concnodules%>%
  select(plant_species, plant_status, nodule_ID, concatenated_OTU_98)%>%
  filter(concatenated_OTU_98!='')%>%
  mutate(interaction=1)%>%
  spread(key=concatenated_OTU_98, value=interaction, fill=0)

#subset out only Bay Area and brady japonicum strains
concbjBayAreaInteractionMatrix <- concinteractionMatrix%>%
  select(-conc_019, -conc_011, -conc_021, -conc_022)%>%
  filter(nodule_ID!='X01z', plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
  #get summary interaction matrix (sum of interactions by species)
  gather(key=concatenated_OTU_98, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
  group_by(plant_status, plant_species, concatenated_OTU_98)%>%
  summarise(interaction=sum(interaction))%>%
  spread(key=concatenated_OTU_98, value=interaction)