library(maps)
library(maptools)
library(mapdata)
library(ggplot2)
library(ggmap)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\DNA work\\molecular data analysis\\strain data')

dat <- read.csv('La Pierre_invasion molecular manuscript_site coordinates.csv')%>%
  filter(!is.na(lat))

lat<-as.numeric(dat$lat)
lon<-as.numeric(dat$lon)
site<-dat$Site_ID

#get map data
caData <- map_data("state")%>%
  filter(region=="california")

#set bounds as min and max of the data (left, bottom, right, top of map)
bounds <- (c(min(lon)-0.1,min(lat)-0.1,max(lon)+0.1,max(lat)+0.1))

#make map
map <-  get_map(location=bounds, maptype="watercolor", source="stamen")

ggmap(map) +
  geom_point(data=dat, aes(x=lon,y=lat), color='black', size=3) +
  geom_text(data=dat, aes(label=as.character(Site_ID)), hjust=-0.3, angle=40)
#export at 1000x1000