library(picante)
library(plyr)
library(ggplot2)
library(reshape2)
library(grid)
library(vegan)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=32, vjust=-0.35), axis.text.x=element_text(size=28),
             axis.title.y=element_text(size=32, angle=90, vjust=0.5), axis.text.y=element_text(size=28),
             plot.title = element_text(size=24, vjust=2),
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
source('Invasive Shrub Molecular Data Analysis\\Invasive-Shrub-Molecular-Data\\invasive shrub_concatenated_data management.R')

#drop reference strains from tree
BjBayAreaconctree<-drop.tip(BjBayAreaconctree, c('Mesorhizobium_ciceri_USDA3383_ITS_AF345262_1', 'Rhizobium_leguminosarum_X01z', 'B_elkanii_USDA76_ITS_AF345254_1', 'B_yuanmingense_LMG21827_ITS_AY386734_1', 'B_canariense_BTA1_ITS_AY386708_1', 'B_liaoningense_USDA3622_ITS_AF345256_1'))
plot.phylo(BjBayAreaconctree, use.edge.length=F)

concbjBayAreaInteractionMatrix <- concbjBayAreaInteractionMatrix%>%
  #create a column combining status and species
  mutate(plant_code=ifelse(plant_species=='Acmispon angustissimus', 'ACAN', ifelse(plant_species=='Acmispon glaber', 'ACGL', ifelse(plant_species=='Acmispon heermannii', 'ACHE', ifelse(plant_species=='Acmispon micranthus', 'ACMI', ifelse(plant_species=='Acmispon strigosus', 'ACST', ifelse(plant_species=='Genista monspessulana', 'GEMO', ifelse(plant_species=='Lupinus bicolor', 'LUBI', ifelse(plant_species=='Spartium junceum', 'SPJU', ifelse(plant_species=='Lupinus arboreous', 'LUAR', 'ULEU'))))))))), type=as.character(paste(plant_code, plant_status, sep='_')))

#this needs the fixed tree with new OTU names in it
#calculate diversity metrics
PD <- pd(concbjBayAreaInteractionMatrix[,c(-1:-2, -22:-23)], BjBayAreaconctree, include.root=F) #phylogenetic diversity

MPD <- mpd(concbjBayAreaInteractionMatrix[,c(-1:-2, -22:-23)], cophenetic(BjBayAreaconctree)) #mean pairwise distance

phydist <- cophenetic(BjBayAreaconctree) #matrix of phylogenetic distances among all pairs

ses.mpd <- ses.mpd(concbjBayAreaInteractionMatrix[,c(-1:-2, -22:-23)], phydist, null.model="taxa.labels", abundance.weighted=T, runs=100)%>%
  select(-runs, -ntaxa) #net relatedness index

ses.mntd <- ses.mntd(concbjBayAreaInteractionMatrix[,c(-1:-2, -22:-23)], phydist, null.model="taxa.labels", abundance.weighted=T, runs=100) #nearest taxon index

#combine diversity metrics into one dataset
diversity <- cbind(PD, MPD, ses.mpd, ses.mntd)%>%
  mutate(NRI=mpd.obs.z*(-1), NTI=mntd.obs.z*(-1))%>%
  select(PD, SR, MPD, NRI, NTI)

#get back species type
diversity <- cbind(concbjBayAreaInteractionMatrix$type, concbjBayAreaInteractionMatrix$plant_status, concbjBayAreaInteractionMatrix$plant_species, diversity)
diversity[is.na(diversity)] <- 0
names(diversity)[names(diversity)=='concbjBayAreaInteractionMatrix$type'] <- 'type'
names(diversity)[names(diversity)=='concbjBayAreaInteractionMatrix$plant_status'] <- 'plant_status'
names(diversity)[names(diversity)=='concbjBayAreaInteractionMatrix$plant_species'] <- 'plant_species'


###student's ttest (assumes equal variances)
t.test(SR~plant_status, diversity, var.equal=T) #SR not different, t=2.0837, p=0.07567, df=7
t.test(PD~plant_status, diversity, var.equal=T) #PD not different, t=0.56525, p=0.5896, df=7
t.test(MPD~plant_status, diversity, var.equal=T) #MPD not different, t=-0.13967, p=0.8929, df=7
t.test(NRI~plant_status, diversity, var.equal=T) #NRI is different, t=3.9128, p=0.005802, df=7
t.test(NTI~plant_status, diversity, var.equal=T) #NTI not different, t=1.8548, p=0.106, df=7

#PD and MPD
PDplot<-ggplot(data=barGraphStats(data=diversity, variable="PD", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.2, 0.05), name="Phylogenetic Diversity") +
  scale_x_discrete(limits=c('native', 'invasive')) +
  coord_cartesian(ylim=c(0, 0.2)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none")

MPDplot<-ggplot(data=barGraphStats(data=diversity, variable="MPD", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02), name="Mean Pairwise Distance") +
  scale_x_discrete(limits=c('native', 'invasive')) +
  coord_cartesian(ylim=c(0, 0.08)) +
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
  scale_y_continuous(breaks=seq(0, 2, 0.2), name="Net Relatedness Index") +
  coord_cartesian(ylim=c(0, 2)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none") +
  annotate('text', x=1, y=0.8, label='a', size=10) +
  annotate('text', x=2, y=1.9, label='b', size=10)

NTIplot<-ggplot(data=barGraphStats(data=diversity, variable="NTI", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_x_discrete(limits=c('native', 'invasive')) +
  scale_y_continuous(breaks=seq(0, 1.6, 0.2), name="Nearest Taxon Index") +
  coord_cartesian(ylim=c(0, 1.6)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none")

pushViewport(viewport(layout=grid.layout(1,2)))
print(NRIplot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(NTIplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))


###Chao richness estimates
plantStrainRichness <- specnumber(concbjBayAreaInteractionMatrix[,3:21]) #gives strain richness for each plant species based on abundance data

totalStrainPool <- specpool(concbjBayAreaInteractionMatrix[,3:21]) #estimate total strain pool across all plant species

speciesStrainRichness <- data.frame(estimateR(concbjBayAreaInteractionMatrix[,3:21])) #estimates Chao strain richness for each plant species

speciesStrainRichness <- cbind(row_names=rownames(speciesStrainRichness), speciesStrainRichness)%>%
  gather(key=type, value=estimate, -row_names)%>%
  filter(type!='row_names')%>%
  mutate(plant_species=ifelse(type=='X1', 'Acmispon angustissimus', ifelse(type=='X2', 'Genista monspessulana', ifelse(type=='X3', 'Spartium junceum', ifelse(type=='X4', 'Ulex europaeus', ifelse(type=='X5', 'Acmispon glaber', ifelse(type=='X6', 'Acmispon heermannii', ifelse(type=='X7', 'Acmispon micranthus', ifelse(type=='X8', 'Acmispon strigosus', ifelse(type=='X9', 'Lupinus arboreous', 'Lupinus bicolor'))))))))))%>%
  spread(key=row_names, value=estimate)%>%
  mutate(plant_status=ifelse(plant_species=='Acmispon angustissimus', 'invasive', ifelse(plant_species=='Genista monspessulana', 'invasive', ifelse(plant_species=='Spartium junceum', 'invasive', ifelse(plant_species=='Ulex europaeus', 'invasive', 'native')))))
  
###ttest for Chao richness
t.test(S.chao1~plant_status, speciesStrainRichness, var.equal=T) #Chao richness estimate not different, t=0.89117, p=0.4024, df=7

chaoPlot <- ggplot(data=barGraphStats(data=speciesStrainRichness, variable="S.chao1", byFactorNames=c("plant_status")), aes(x=plant_status, y=mean, fill=plant_status)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_x_discrete(limits=c('native', 'invasive')) +
  scale_y_continuous(breaks=seq(0, 10, 2), name="Chao Richness Estimate") +
  coord_cartesian(ylim=c(0, 10)) +
  xlab("Plant Status") +
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  theme(legend.position="none")


#figure of Chao richness, PD, and MPD
pushViewport(viewport(layout=grid.layout(1,3)))
print(chaoPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(PDplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(MPDplot, vp=viewport(layout.pos.row=1, layout.pos.col=3))

