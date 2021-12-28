
setwd(paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/ChinaAccessibility/Models/AccessProfiles'))

library(raster)
library(ggplot2)
library(sf)
library(dplyr)

source('network.R')
source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

cities = c("Beijing","Chengdu","Chongqing","Guangzhou","Nanjing","Shanghai","Shenzhen","Wuhan","Xian")
# cities = c("Shanghai","Shenzhen","Wuhan","Xian")

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


getPopPoints <- function(popraster,currentbbox,year){
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
  return(poppoints)
}




#access=c();
#decays=c(); # compute access with fixed decay here -> 60min
#years=c();ccities=c();radius=c()


decay = 60

for(city in cities){
  
  # compute at different network stages -> 1990 (! no nw?) -> 2000 (or min date); 2015; 2030?
  #load(paste0(datadir,'processed/',city,'_multiyearGHSL.RData'))
  #show(summary(E(fullgraph)$year))
  
  currentextent = fuas[fuas$FUA_area%in%fuaids[[city]]&fuas$eFUA_name==fuanames[[city]],]
  currentbbox = st_bbox(currentextent)
  
  xcenter = st_coordinates(st_centroid(currentextent))[1];ycenter=st_coordinates(st_centroid(currentextent))[2]
  
  trgraph=addTransportationLayer(link_layer = paste0(datadir,ptlayers[[city]]),speed=0.0012,snap=500,
                                 e_attr_names=c("year"),reprojection=crs(pop90))
  if(city=='Guangzhou'){
    # specific case of guangzhou-foshan
    trgraph=addTransportationLayer(g=trgraph,link_layer = paste0(datadir,'foshan'),speed=0.0012,snap=500,
                                   e_attr_names=c("year"),reprojection=crs(pop90))
  }
  
  fullgraph=trgraph
  #for(year in c("1990","2000","2015")){
  # 
  #  popraster = populations[[year]]
  #  poppoints = getPopPoints(popraster,currentbbox)
  #  # issue with calling addPoints with null empty graph heuristic
  #  fullgraph = addAdministrativeLayer(fullgraph,poppoints,connect_speed = 0.006,attributes=list("pop"="pop","id"="id","year"="year"),empty_graph_heuristic="NA")
  #}
  #save(fullgraph,file=paste0(datadir,'processed/',city,'_multiyearGHSL.RData'))
  
  # year=2030
  
  # data for profile plots
  pops=c();accesses=c();dists=c();years=c()
  
  for(year in c(2000,2015,2030)){
    currentyear=year
    if(year==2000){currentyear=max(c(min(E(fullgraph)$year,na.rm = T),year))}
    show(currentyear)
    currentgraph = subgraph.edges(fullgraph,which(E(fullgraph)$year<=currentyear),delete.vertices = T)
    # need to re-add pop points - full compute was not necessary
    popyear = ifelse(currentyear<2015,"2000","2015")
    popraster = populations[[popyear]]
    poppoints = getPopPoints(popraster,currentbbox,year)
    currentgraph = addAdministrativeLayer(currentgraph,poppoints,connect_speed = 0.002,attributes=list("pop"="pop","id"="id","year"="year"),empty_graph_heuristic="NA")
    dmat = distances(graph = currentgraph,v=V(currentgraph)[!is.na(V(currentgraph)$id)],to=V(currentgraph)[!is.na(V(currentgraph)$id)],weights = E(currentgraph)$speed*E(currentgraph)$length)
    rownames(dmat)<-V(currentgraph)$id[!is.na(V(currentgraph)$id)];colnames(dmat)<-V(currentgraph)$id[!is.na(V(currentgraph)$id)]
    
    access = computeAccess(accessorigdata = data.frame(id=rownames(dmat),var=rep(1,nrow(dmat)),year=rep(currentyear,nrow(dmat))),
                           #accessdestdata = data.frame(id=rownames(dmat),var=as.numeric(V(currentgraph)$pop[!is.na(V(currentgraph)$pop)]),year=rep(currentyear,nrow(dmat))) ,
                           accessdestdata = data.frame(id=rownames(dmat),var=rep(1,nrow(dmat)),year=rep(currentyear,nrow(dmat))),
                          matfun=exp(-dmat/decay)
    )
    
    
    ### Access map
    xr=xres(popraster);yr=yres(popraster)
    sppol = SpatialPolygons(apply(poppoints@data,1,function(r){r = as.numeric(r)
        Polygons(list(Polygon(matrix(c(r[1]-xr/2,r[2]+yr/2,r[1]+xr/2,r[2]+yr/2,r[1]+xr/2,r[2]-yr/2,r[1]-xr/2,r[2]-yr/2,r[1]-xr/2,r[2]+yr/2),ncol=2,byrow = T))),ID = as.character(r[5]))
      })
      ,proj4string = crs(popraster))
    spdat=poppoints@data;rownames(spdat)<-sapply(sppol@polygons,function(p){p@ID})
    spdf=SpatialPolygonsDataFrame(sppol,data = spdat)
    
    nwlayer = spTransform(readOGR(datadir,ptlayers[[city]]),crs(spdf))
    
    # left_join(spdf@data,access,by=c('id'='id'))
    
    trnw=nwlayer[as.numeric(nwlayer$year)<=currentyear,]
    
    # ! there may be some connectivity issues with Guangzhou network (connection GUangzhou - Foshan?)
    map(data=access,layer=spdf,spdfid="id",dfid="id",variable="var",
        filename=paste0(resdir,'access_',city,'_',year,'.png'),title=paste0(city," (",currentyear,")"),
        legendtitle = "Accessibility",extent=spdf,
        nclass=8,
        width=15,height=15,palette='div',lwd=0.01,
        pngResolution = 300,
        additionalLinelayers=list(
          list(trnw,'blue')
        ),
        withScale=NULL,
        legendPosition = "bottomright"
    )
    
    
    ######
    # profiles
    allpoints = left_join(spdf@data,access,by=c('id'='id'))
    centerdist= sqrt((allpoints$x - xcenter)^2 + (allpoints$y - ycenter)^2)
    distquantiles = c(0,quantile(centerdist,seq(0.05,1.0,0.05)))
    
    for(i in 2:length(distquantiles)){
      inds = which(centerdist>distquantiles[i-1]&centerdist<=distquantiles[i]);
      pops = append(pops,allpoints$pop[inds]);accesses=append(accesses,allpoints$var[inds]);
      dists=append(dists,rep((distquantiles[i]+distquantiles[i-1])/2,length(inds)))
      years=append(years,rep(year,length(inds)))
    }
    
  }
  
  plotdata=data.frame(population = pops,accessibility = accesses,distance=dists,year=as.character(years))
  plotdata=plotdata[plotdata$distance<100000,]
  # pop profile
  g=ggplot(data = plotdata,aes(x=distance,y=population,group=year,color=year))
  #g+geom_point(pch='.')+geom_smooth() # points too broad
  ggsave(plot = g+geom_smooth()+stdtheme, filename = paste0(resdir,'popprofile-smooth_',city,'.png'),width=25,height=20,units='cm')
  
  #sdata = as_tibble(data.frame(population = pops,distance=dists,year=as.character(years))) %>% group_by(distance,year) %>% summarise(popsd=sd(population), population = mean(population))
  #g=ggplot(data = sdata,aes(x=distance,y=population,group=year,color=year))
  ##g+geom_point(pch='.')+geom_smooth()
  #ggsave(plot = g+geom_point()+geom_errorbar(aes(ymin=population-popsd,ymax=population+popsd))+stdtheme, filename = paste0(resdir,'popprofile-sd_',city,'.png'),width=25,height=20,units='cm')
  
  # access profile
  g=ggplot(data = plotdata,aes(x=distance,y=accessibility,group=year,color=year))
  ggsave(plot = g+geom_smooth()+ggtitle(city)+stdtheme, filename = paste0(resdir,'accessprofile-smooth_',city,'.png'),width=25,height=20,units='cm')
  
  
}









