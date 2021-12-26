
setwd(paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/ChinaAccessibility/Models/AccessProfiles'))

library(raster)
library(ggplot2)
library(sf)

source(paste0(Sys.getenv('CN_HOME'),'/Models/TransportationNetwork/NetworkAnalysis/network.R'))
source(paste0(Sys.getenv('CN_HOME'),'/Models/SpatioTempCausality/functions.R'))
source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

cities = c("Beijing","Chengdu","Chongqing","Guangzhou","Nanjing","Shanghai","Shenzhen","Wuhan","Xian")

ptlayers = list("Beijing"="beijing","Chengdu"="chengdu","Chongqing"="chongqing","Guangzhou"="guangzhou","Nanjing"="nanjing","Shanghai"="shanghai","Shenzhen"="shenzhen","Wuhan"="wuhan","Xian"="xian")


fuas <- st_read(paste0(Sys.getenv('CS_HOME'),'/Data/JRC_EC/GHS/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.gpkg'))

# specific selection for guangzhou/shenzhen -> use both name and id
fuaids = list("Beijing"=c(5536),"Chengdu"=c(4034),"Chongqing"=c(2267),"Guangzhou"=c(4413,2151,681),"Nanjing"=c(2419),"Shanghai"=c(5878),"Shenzhen"=c(583,247,1573),"Wuhan"=c(2682),"Xian"=c(2083))
fuanames = list("Beijing"="Beijing","Chengdu"="Chengdu","Chongqing"="Chongqing","Guangzhou"="Guangzhou","Nanjing"="Nanjing","Shanghai"="Shanghai","Shenzhen"="Guangzhou","Wuhan"="Wuhan","Xian"="Xi'an")
#for(city in cities){show(city);show(st_bbox(fuas[fuas$FUA_area%in%fuaids[[city]]&fuas$eFUA_name==fuanames[[city]],]))}

datadir = paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/ChinaAccessibility/Data/metro/')

resdir=paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/ChinaAccessibility/Results/AccessProfiles/');dir.create(resdir)

pop90 = raster(x = paste0(Sys.getenv('CS_HOME'),'/Data/JRC_EC/GHS/GHS_POP_GPW41990_GLOBE_R2015A_54009_1k_v1_0/GHS_POP_GPW41990_GLOBE_R2015A_54009_1k_v1_0.tif'))
pop00 = raster(x = paste0(Sys.getenv('CS_HOME'),'/Data/JRC_EC/GHS/GHS_POP_GPW42000_GLOBE_R2015A_54009_1k_v1_0/GHS_POP_GPW42000_GLOBE_R2015A_54009_1k_v1_0.tif'))
pop15 = raster(x = paste0(Sys.getenv('CS_HOME'),'/Data/JRC_EC/GHS/GHS_POP_GPW42015_GLOBE_R2015A_54009_1k_v1_0/GHS_POP_GPW42015_GLOBE_R2015A_54009_1k_v1_0.tif'))
populations = list("1990" = pop90, "2000" = pop00, "2015" = pop15)


# construct graphs
# city = 'Xian'
for(city in cities){
  
  show(city)
  
  currentextent = fuas[fuas$FUA_area%in%fuaids[[city]]&fuas$eFUA_name==fuanames[[city]],]
  currentbbox = st_bbox(currentextent)
  
  trgraph=addTransportationLayer(link_layer = paste0(datadir,ptlayers[[city]]),speed=0.0012,snap=200,
                                 e_attr_names=c("year"),reprojection=crs(pop90))
  if(city=='Guangzhou'){
    # specific case of guangzhou-foshan
    trgraph=addTransportationLayer(g=trgraph,link_layer = paste0(datadir,'foshan'),speed=0.0012,snap=200,
                                   e_attr_names=c("year"),reprojection=crs(pop90))
  }
 
  
  # add population points (for each year)
  fullgraph=trgraph
  for(year in c("1990","2000","2015")){
    show(year)
    popraster = populations[[year]]
    maxrow = rowFromY(popraster,currentbbox$ymin);minrow = rowFromY(popraster,currentbbox$ymax)
    mincol = colFromX(popraster,currentbbox$xmin);maxcol = colFromX(popraster,currentbbox$xmax)
  
    popdf = data.frame(cbind(
      x = rep(xFromCol(popraster,mincol:maxcol),maxrow-minrow+1),
      y = c(t(matrix(rep(yFromRow(popraster,minrow:maxrow),maxcol-mincol+1),ncol=maxcol-mincol+1,byrow = F))),
      pop=getValuesBlock(popraster,row=minrow,nrows=maxrow-minrow+1,col=mincol,ncols=maxcol-mincol+1))
    )
  
    popdf$pop[is.na(popdf$pop)]=0
    popdf$year=rep(as.numeric(year),nrow(popdf))
    poppoints = SpatialPointsDataFrame(popdf[popdf$pop>0,c("x","y")],popdf[popdf$pop>0,])
    poppoints$id = as.character(1:length(poppoints))
    # issue with calling addPoints with null empty graph heuristic
    fullgraph = addAdministrativeLayer(fullgraph,poppoints,connect_speed = 0.006,attributes=list("pop"="pop","id"="id","year"="year"),empty_graph_heuristic="NA")
  }
  
  show(max(E(fullgraph)$year,na.rm=T))
  
  save(fullgraph,file=paste0(datadir,'processed/',city,'_multiyearGHSL.RData'))
   
}





avgrelaccesses=c();
#decays=c();
years=c();ccities=c()









