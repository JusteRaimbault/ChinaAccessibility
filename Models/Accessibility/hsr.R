
setwd(paste0(Sys.getenv('CS_HOME'),'/ChinaAccessibility/Models/Accessibility'))

library(raster)
library(ggplot2)
library(dplyr)

source(paste0(Sys.getenv('CN_HOME'),'/Models/TransportationNetwork/NetworkAnalysis/network.R'))
source(paste0(Sys.getenv('CN_HOME'),'/Models/SpatioTempCausality/functions.R'))
source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

datadir = paste0(Sys.getenv('CS_HOME'),'/ChinaAccessibility/Data/')

resdir=paste0(Sys.getenv('CS_HOME'),'/ChinaAccessibility/Results/HSR/');dir.create(resdir)

trainlines = readOGR(paste0(datadir,'hsr/'),layer = 'lines')

cities=read.csv(paste0(datadir,'population/ChinaCities_Swerts.csv'),sep=';')
# remove rows with na (after 2000)
cities=cities[!is.na(cities$X2000)&!is.na(cities$X2010)&!is.na(cities$Lat)&!is.na(cities$Long),]
citiespoints = SpatialPointsDataFrame(cities[,c("Long","Lat")],cities,proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
citiespoints$id = as.character(1:length(citiespoints))
# strange outlyer far in the north -> filter
citiespoints = citiespoints[citiespoints$Lat<100,]
# reproject in same crs as train nw
citiespoints=spTransform(citiespoints,crs(trainlines))
# filter on city populations (to have decent nw)
citiespoints=citiespoints[citiespoints$X2010>500000,]


# train layer
network = addTransportationLayer(stations_layer = paste0(datadir,'hsr/stations.shp'),
                                link_layer = paste0(datadir,'hsr/lines.shp'),snap = 10000 # high snap as projection is in m
                                ,e_attr_names = c("speed","opening")
)
# ! speed is in km/h here

#plot(network,vertex.size=2,vertex.label=NA,vertex.color=components(network)$membership)

fullnetwork = addAdministrativeLayer(g=network,
                                     admin_layer=citiespoints,
                                     connect_speed = 100,
                                     attributes=list("pop2000"="X2000","pop2010"="X2010","id"="id"),
                                     empty_graph_heuristic="full"
                                     )
E(fullnetwork)$opening[is.na(E(fullnetwork)$opening)]=0

#sizes=log(as.numeric(V(fullnetwork)$pop2000));sizes[is.na(sizes)]=0
#plot(fullnetwork,vertex.size=sizes,vertex.label=NA,vertex.color=components(network)$membership)

# compute travel times between cities
years = unique(E(fullnetwork)$opening);years=as.character(years[years>0])
#timemats=list()
citiesind=which(!is.na(V(fullnetwork)$pop2010))
totalpop=sum(as.numeric(V(currentnetwork)$pop2010[citiesind]))

accesses=c();cyears=c();decays=c();cityids=c()
for(year in sort(years)){
  show(year)
  currentnetwork=subgraph.edges(fullnetwork,which(E(fullnetwork)$opening<=as.numeric(year)))
  citiesind=which(!is.na(V(currentnetwork)$pop2010))
  currenttmat=distances(currentnetwork,v = citiesind,to= citiesind,weights = E(currentnetwork)$length/ (1000*E(currentnetwork)$speed))
  #timemats[[year]]=currenttmat
  show(summary(c(currenttmat)))
  
  for(decay in c(6,12,24)){
    access = computeAccess(accessorigdata = data.frame(id=rownames(currenttmat),var=rep(1,nrow(currenttmat)),year=rep(year,nrow(currenttmat))),
                           accessdestdata = data.frame(id=rownames(currenttmat),var=as.numeric(V(currentnetwork)$pop2010[citiesind]),year=rep(year,nrow(currenttmat))) ,
                           matfun=exp(-currenttmat/decay)
    )
    accesses=append(accesses,access$var/totalpop)
    cyears=append(cyears,access$year);decays=append(decays,rep(decay,length(access$var)));
    cityids=append(cityids,citiesind)
  }
  show(summary(access$var/totalpop))
}

# TODO : should add city names to be able to follow specific trajectories

res = data.frame(access=accesses,year=cyears,decay=decays,city=cityids)

g=ggplot(res,aes(x=access,group=interaction(year,decay),color=year,linetype=as.character(decay)))
g+geom_density()+scale_linetype_discrete(name="Decay (h)")+xlab("Relative accessibility")+ylab("Density")+stdtheme
ggsave(file=paste0(resdir,'access_hist_all.png'),width=20,height=15,units='cm')


g=ggplot(res,aes(x=access,group=year,color=year))
g+geom_density()+facet_wrap(~decay)+scale_linetype_discrete(name="Decay (h)")+xlab("Relative accessibility")+ylab("Density")+stdtheme
ggsave(file=paste0(resdir,'access_hist_facet.png'),width=20,height=15,units='cm')


# select years and decay as we do not see anything
g=ggplot(res[as.character(res$year)%in%c("2010","2015","2020")&res$decay%in%c(6,12),],aes(x=access,group=year,color=year))
g+geom_density()+facet_wrap(~decay)+xlab("Relative accessibility")+ylab("Density")+stdtheme
ggsave(file=paste0(resdir,'access_hist_facet_selected.png'),width=20,height=15,units='cm')

#g=ggplot(res,aes(x=year,y=access,group=decay,color=decay))
#g+geom_point(pch='.')+stat_smooth(span = 0.1,method = 'loess')

hierarchy<-function(x){reg=lm(data=data.frame(y=log(sort(x,decreasing = T)),x=log(1:length(x))),y~x);return(-reg$coefficients[2])}

## summary
sres=as.tbl(res)%>%group_by(year,decay)%>%summarise(sdaccess=sd(access),hierarchy=hierarchy(access),access=mean(access))

g=ggplot(sres,aes(x=year,y=access,ymin=access-sdaccess,ymax=access+sdaccess,color=as.character(decay),group=as.character(decay)))
g+geom_point()+geom_line()+geom_errorbar()+xlab("Year")+ylab("Relative accessibility")+scale_color_discrete(name="Decay (h)")+stdtheme
ggsave(file=paste0(resdir,'access_summary.png'),width=30,height=20,units='cm')

# hierarchy
g=ggplot(sres,aes(x=year,y=hierarchy,color=as.character(decay),group=as.character(decay)))
g+geom_point()+geom_line()+xlab("Year")+ylab("Hierarchy of accessibility")+scale_color_discrete(name="Decay (h)")+stdtheme
ggsave(file=paste0(resdir,'access_hierarchy.png'),width=30,height=20,units='cm')






