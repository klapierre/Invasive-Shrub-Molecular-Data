library(picante)
library(plyr)
library(ggplot2)
library(reshape2)
library(grid)
library(vegan)
library(tidyr)
library(dplyr)
library(bipartite)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=32, vjust=-0.35), axis.text.x=element_text(size=28),
             axis.title.y=element_text(size=32, angle=90, vjust=0.5), axis.text.y=element_text(size=28),
             plot.title = element_text(size=24, vjust=2),
             axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

#source data management code
source('Invasive Shrub Molecular Data Analysis\\Invasive-Shrub-Molecular-Data\\invasive shrub_nifd_data management.R')

#drop reference strains from tree
BjBayAreanifdtree<-drop.tip(BjBayAreanifdtree, c('USDA_38_japonicum', 'USDA_3622_liaoningense', 'USDA_94_elkanii', 'M76_sinorhizobium', 'NZP514_rhizobium'))
plot.phylo(BjBayAreanifdtree, use.edge.length=F)

nifdBjBayAreaInteractionMatrix <- nifdBjBayAreaInteractionMatrix%>%
  #create a column combining status and species
  mutate(plant_code=ifelse(plant_species=='Acmispon angustissimus', 'ACAN', ifelse(plant_species=='Acmispon glaber', 'ACGL', ifelse(plant_species=='Acmispon heermannii', 'ACHE', ifelse(plant_species=='Acmispon micranthus', 'ACMI', ifelse(plant_species=='Acmispon strigosus', 'ACST', ifelse(plant_species=='Genista monspessulana', 'GEMO', ifelse(plant_species=='Lupinus bicolor', 'LUBI', ifelse(plant_species=='Spartium junceum', 'SPJU', ifelse(plant_species=='Lupinus arboreous', 'LUAR', 'ULEU'))))))))), type=as.character(paste(plant_code, plant_status, sep='_')))


#calculate diversity metrics
PD <- pd(nifdBjBayAreaInteractionMatrix[,c(-1:-2, -28:-29)], BjBayAreanifdtree, include.root=F) #phylogenetic diversity

MPD <- mpd(nifdBjBayAreaInteractionMatrix[,c(-1:-2, -28:-29)], cophenetic(BjBayAreanifdtree)) #mean pairwise distance

phydist <- cophenetic(BjBayAreanifdtree) #matrix of phylogenetic distances among all pairs

ses.mpd <- ses.mpd(nifdBjBayAreaInteractionMatrix[,c(-1:-2, -28:-29)], phydist, null.model="taxa.labels", abundance.weighted=T, runs=100)%>%
  select(-runs, -ntaxa) #net relatedness index

ses.mntd <- ses.mntd(nifdBjBayAreaInteractionMatrix[,c(-1:-2, -28:-29)], phydist, null.model="taxa.labels", abundance.weighted=T, runs=100) #nearest taxon index

#combine diversity metrics into one dataset
diversity <- cbind(PD, MPD, ses.mpd, ses.mntd)%>%
  mutate(NRI=mpd.obs.z*(-1), NTI=mntd.obs.z*(-1))%>%
  select(PD, SR, MPD, NRI, NTI)

#get back species type
diversity <- cbind(nifdBjBayAreaInteractionMatrix$type, nifdBjBayAreaInteractionMatrix$plant_status, nifdBjBayAreaInteractionMatrix$plant_species, diversity)
diversity[is.na(diversity)] <- 0
names(diversity)[names(diversity)=='nifdBjBayAreaInteractionMatrix$type'] <- 'type'
names(diversity)[names(diversity)=='nifdBjBayAreaInteractionMatrix$plant_status'] <- 'plant_status'
names(diversity)[names(diversity)=='nifdBjBayAreaInteractionMatrix$plant_species'] <- 'plant_species'


###student's ttest (assumes equal variances)
t.test(SR~plant_status, diversity, var.equal=T) #SR different, t=5.5979, p=0.001383, df=6
t.test(PD~plant_status, diversity, var.equal=T) #PD different, t=3.1932, p=0.01876, df=6
t.test(MPD~plant_status, diversity, var.equal=T) #MPD not different, t=0.76807, p=0.4716, df=6
t.test(NRI~plant_status, diversity, var.equal=T) #NRI different, t=-2.5108, p=0.04585, df=6
t.test(NTI~plant_status, diversity, var.equal=T) #NTI not different, t=-1.7147, p=0.1372, df=6

#PD and MPD
PDplot<-ggplot(data=barGraphStats(data=diversity, variable="PD", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.25, 0.05), name="PD") +
  scale_x_discrete(limits=c('native', 'invasive')) +
  coord_cartesian(ylim=c(0, 0.25)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none") +
  annotate('text', x=1, y=0.11, label='a', size=10) +
  annotate('text', x=2, y=0.23, label='b', size=10)

MPDplot<-ggplot(data=barGraphStats(data=diversity, variable="MPD", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.05, 0.01), name="MPD") +
  scale_x_discrete(limits=c('native', 'invasive')) +
  coord_cartesian(ylim=c(0, 0.05)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none")

pushViewport(viewport(layout=grid.layout(1,2))) 
print(PDplot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(MPDplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))

#NRI and NTI
NRIplot<-ggplot(data=barGraphStats(data=diversity, variable="NRI", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_x_discrete(limits=c('native', 'invasive')) +
  scale_y_continuous(breaks=seq(0, 1.6, 0.2), name="NRI") +
  coord_cartesian(ylim=c(0, 1.6)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none") +
  annotate('text', x=1, y=1.55, label='a', size=10) +
  annotate('text', x=2, y=0.75, label='b', size=10)

NTIplot<-ggplot(data=barGraphStats(data=diversity, variable="NTI", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_x_discrete(limits=c('native', 'invasive')) +
  scale_y_continuous(breaks=seq(0, 1.2, 0.2), name="NTI") +
  coord_cartesian(ylim=c(0, 1.2)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none")

pushViewport(viewport(layout=grid.layout(1,2)))
print(NRIplot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(NTIplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))


###Chao richness estimates
plantStrainRichness <- specnumber(nifdBjBayAreaInteractionMatrix[,3:27]) #gives strain richness for each plant species based on abundance data

totalStrainPool <- specpool(nifdBjBayAreaInteractionMatrix[,3:27]) #estimate total strain pool across all plant species

speciesStrainRichness <- data.frame(estimateR(nifdBjBayAreaInteractionMatrix[,3:27])) #estimates Chao strain richness for each plant species 

speciesStrainRichness <- cbind(row_names=rownames(speciesStrainRichness), speciesStrainRichness)%>%
  gather(key=type, value=estimate, -row_names)%>%
  filter(type!='row_names')%>%
  mutate(plant_species=ifelse(type=='X1', 'Genista monspessulana', ifelse(type=='X2', 'Spartium junceum', ifelse(type=='X3', 'Ulex europaeus', ifelse(type=='X4', 'Acmispon glaber', ifelse(type=='X5', 'Acmispon heermannii', ifelse(type=='X6', 'Acmispon micranthus', ifelse(type=='X7', 'Acmispon strigosus', ifelse(type=='X8', 'Lupinus arboreous', 'Lupinus bicolor')))))))))%>%
  spread(key=row_names, value=estimate)%>%
  mutate(plant_status=ifelse(plant_species=='Genista monspessulana', 'invasive', ifelse(plant_species=='Spartium junceum', 'invasive', ifelse(plant_species=='Ulex europaeus', 'invasive', 'native'))))
  
###ttest for Chao richness
t.test(S.chao1~plant_status, speciesStrainRichness, var.equal=T) #Chao richness estimate not different, t=1.3327, p=0.2244, df=7

chaoPlot <- ggplot(data=barGraphStats(data=speciesStrainRichness, variable="S.chao1", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_x_discrete(limits=c('native', 'invasive')) +
  scale_y_continuous(breaks=seq(0, 12, 2), name="Chao Richness Estimate") +
  coord_cartesian(ylim=c(0, 12)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none")


#figure of Chao richness, PD, and NRI
pushViewport(viewport(layout=grid.layout(1,3)))
print(chaoPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(PDplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(NRIplot, vp=viewport(layout.pos.row=1, layout.pos.col=3))





# ###network analysis - determine specialization and generalization
# matrix <- (nifdBjBayAreaInteractionMatrix[,2:47])
# rownames(matrix) <- matrix$plant_species
# matrix <-  matrix%>%
#   select(-plant_species)
# 
# visweb(t(matrix))
# PDI_NifD <- as.data.frame(PDI(t(matrix)))%>%
#   add_rownames('plant')%>%
#   mutate(plant_status=ifelse(plant=='Genista monspessulana', 'invasive', ifelse(plant=='Spartium junceum', 'invasive', ifelse(plant=='Ulex europaeus', 'invasive', ifelse(plant=='Acmispon glaber', 'invasive', 'native')))))
# 
# t.test(PDI(t(matrix))~plant_status, PDI_NifD, var.equal=T) #SR different, t=-1.228, p=0.2654, df=6
# 
# 
# #network test from both directions
# networkIndices <- specieslevel(t(matrix))
# 
# #plot the interaction web
# 
# plantcolors<-c("#00990099", "#00990099", "#00990099","#00990099", "#FF990099", "#FF990099", "#FF990099",  "#FF990099", "#FF990099")
# straincolors<-c(rep("black"))
# 
# 
# plotweb(t(matrix), bor.col.interaction=plantcolors,  col.interaction=plantcolors, arrow="no", method="normal", labsize=1.5,
#         text.rot=90, low.lab.dis=NULL,
#         y.width.high=0.005,y.width.low=0.005, col.low="black", col.high="black", 
#         bor.col.high="black", bor.col.low="black")















