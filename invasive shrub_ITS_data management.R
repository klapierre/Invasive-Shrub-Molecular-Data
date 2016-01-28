library(phytools)
library(reshape2)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')


#####################
# tree management
#####################

#mr bayes tree: 2 million runs, including native, invasive, and worldwide invasive species
ITStreeNexus<-read.nexus(file='ITS data\\MrBayes tree_01272016\\ITS_consensus_outgroups_01272016.nexus.con.tre')

#makes tree into a .phylo file instead of .newick or .tre
ITStree<-as.phylo(ITStreeNexus) 

#root tree with Mesorhizobium outgroup
ITStreeR<-root(ITStree, "Mesorhizobium_ciceri_USDA3383_ITS_AF345262_1", r=T)

#plot trees
plot.phylo(ITStreeR, use.edge.length=T)
plot.phylo(ITStreeR, use.edge.length=F)

#remove strains that could not be concatenated with nifd
concITStree<-drop.tip(ITStreeR, c('ITS_013', 'ITS_018', 'ITS_020', 'ITS_034', 'ITS_035', 'ITS_036'))
plot.phylo(concITStree, use.edge.length=F)


#making a tree with only Brady japonicum strains
BjITStree<-drop.tip(concITStree, c('ITS_026', 'ITS_028', 'ITS_029', 'ITS_031', 'ITS_032', 'ITS_033'))
plot.phylo(BjITStree, use.edge.length=F)


#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreaITStree<-drop.tip(BjITStree, c('ITS_011', 'ITS_021', 'ITS_027')) 
plot.phylo(BjBayAreaITStree, use.edge.length=F)



#####################
# strain data management
#####################
#read in nodule data
ITSnodules <- read.csv('strain data\\La Pierre_invasion molecular manuscript_strain information_01252016.csv')%>%
  select(plant_species, plant_status, nodule_ID, ITS_OTU_97, concatenated_OTU_98)
  
#create an interaction matrix of strains for each plant species
ITSinteractionMatrix <- ITSnodules%>%
  select(plant_species, plant_status, nodule_ID, ITS_OTU_97, concatenated_OTU_98)%>%
  filter(ITS_OTU_97!='', concatenated_OTU_98!='')%>%
  select(-concatenated_OTU_98)%>%
  mutate(interaction=1)%>%
  spread(key=ITS_OTU_97, value=interaction, fill=0)

#subset out only Bay Area and brady japonicum strains
ITSbjBayAreaInteractionMatrix <- ITSinteractionMatrix%>%
  select(-ITS_032, -ITS_011, -ITS_021, -ITS_027)%>%
  filter(nodule_ID!='ITS_X01z', plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
  #get summary interaction matrix (sum of interactions by species)
  gather(key=ITS_OTU_98, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
  group_by(plant_status, plant_species, ITS_OTU_98)%>%
  summarise(interaction=sum(interaction))%>%
  spread(key=ITS_OTU_98, value=interaction)














