library(plyr)
library(ggplot2)
library(reshape2)
library(grid)
library(tidyr)
library(dplyr)

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

#import MDS data
nifdMDS <- read.csv('PERMANOVA\\La Pierre_invasion molecular manuscript_nifd_MDS results.csv')

#figure
nifdPlot <- ggplot(data=nifdMDS, aes(x=MDS_1, y=MDS_2, colour=plant_status)) +
  geom_point(shape=19, size=15) +
  geom_text(aes(label=plant_code), hjust=0, vjust=0, colour='black', size=10) +
  scale_color_manual(values=c("#FF9900", "#009900")) +
  xlab('MDS 1') + ylab('MDS 2') +
  theme(legend.position="none") +
  annotate('text', x=-1.25, y=1.2, label='(b) nifD', size=12, hjust=0)
  
#putting ITS and nifD MDS plots together
pushViewport(viewport(layout=grid.layout(1,2)))
print(itsMDS, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(nifdMDS, vp=viewport(layout.pos.row=1, layout.pos.col=2))















