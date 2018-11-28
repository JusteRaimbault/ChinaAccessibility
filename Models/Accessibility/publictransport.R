
setwd(paste0(Sys.getenv('CS_HOME'),'/ChinaAccessibility/Models/Accessibility'))

library(raster)
library(ggplot2)

source(paste0(Sys.getenv('CN_HOME'),'/Models/TransportationNetwork/NetworkAnalysis/network.R'))
source(paste0(Sys.getenv('CN_HOME'),'/Models/SpatioTempCausality/functions.R'))
source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

# city names correspond to layers names 
cities = c("beijing","chengdu","chongqing","guangzhou","nanjing","shanghai","shenzhen","wuhan","xian")
#cities=c("chengdu","chongqing")

datadir = paste0(Sys.getenv('CS_HOME'),'/ChinaAccessibility/Data/metro/')

resdir=paste0(Sys.getenv('CS_HOME'),'/ChinaAccessibility/Results/TCAccess/');dir.create(resdir)

popraster = raster(x = paste0(Sys.getenv('CS_HOME'),'/ChinaAccessibility/Data/population/PopulationGrid_China2010.tif'))

avgrelaccesses=c();accesshierarchies=c();decays=c();years=c();ccities=c()
# construct graphs
for(city in cities){
  # city="beijing"
  show(city)
  trgraph=addTransportationLayer(link_layer = paste0(datadir,city),speed=0.0012,snap=200,
                                 e_attr_names=c("year"),reprojection=crs(popraster))
  if(city=='guangzhou'){
    # specific case of guangzhou-foshan
    trgraph=addTransportationLayer(trgraph,link_layer = paste0(datadir,'foshan'),speed=0.0012,snap=200,
                                 e_attr_names=c("year"),reprojection=crs(popraster))
  }
  #components(trgraph)
  #plot(trgraph,vertex.size=0,vertex.label=NA,vertex.color=components(trgraph)$membership)
  # get population points
  maxrow = rowFromY(popraster,min(V(trgraph)$y));minrow = rowFromY(popraster,max(V(trgraph)$y))
  mincol = colFromX(popraster,min(V(trgraph)$x));maxcol = colFromX(popraster,max(V(trgraph)$x))
  popdf = data.frame(cbind(
                      x = rep(xFromCol(popraster,mincol:maxcol),maxrow-minrow+1),
                      y = c(t(matrix(rep(yFromRow(popraster,minrow:maxrow),maxcol-mincol+1),ncol=maxcol-mincol+1,byrow = F))),
                      pop=getValuesBlock(popraster,row=minrow,nrows=maxrow-minrow+1,col=mincol,ncols=maxcol-mincol+1)))
  #g=ggplot(popdf,aes(x=x,y=y,fill=pop))
  #g+geom_raster()
  poppoints = SpatialPointsDataFrame(popdf[popdf$pop>0,c("x","y")],popdf[popdf$pop>0,])
  poppoints$id = as.character(1:length(poppoints))
  #fullgraph = addAdministrativeLayer(trgraph,poppoints,connect_speed = 0.006,attributes=list("pop"="pop","id"="id"))
  #plot(fullgraph,vertex.size=0,vertex.label=NA,edge.width=1/(50*E(fullgraph)$speed))
  # -> for accurate computation of accessibility each year, should redo the connection each year
  
  # year by year
  for(year in unique(E(trgraph)$year)){
    #year=2008
    show(year)
    currentgraph = subgraph.edges(trgraph,which(E(trgraph)$year<=year),delete.vertices = T)
    currentgraph = addAdministrativeLayer(currentgraph,poppoints,connect_speed = 0.006,attributes=list("pop"="pop","id"="id"))
    #plot(currentgraph,vertex.size=0,vertex.label=NA,edge.width=1/(50*E(fullgraph)$speed))
    
    dmat = distances(graph = currentgraph,weights = E(currentgraph)$speed*E(currentgraph)$length)
    # summary(c(dmat))
    #decay = 30
    for(decay in c(30,60,120)){
    
      access = computeAccess(accessorigdata = data.frame(id=rownames(dmat),var=rep(1,nrow(dmat)),year=rep(year,nrow(dmat))),
                          accessdestdata = data.frame(id=rownames(dmat),var=as.numeric(V(currentgraph)$pop),year=rep(year,nrow(dmat))) ,
                          matfun=exp(-dmat/decay)
                           )
      avgrelaccess = mean(access$var/sum(as.numeric(V(currentgraph)$pop),na.rm = T))
      #hist(access$var/sum(as.numeric(V(currentgraph)$pop),na.rm = T),breaks=100)
      accesshierarchy = -lm(data=data.frame(access=log(sort(access$var,decreasing = T)),rank=log(1:nrow(access))))$coefficients[2]
    
      #res = rbind(res,c(avgrelaccess=avgrelaccess,accesshierarchy=accesshierarchy,decay=decay,city=city,year=year))
      avgrelaccesses=append(avgrelaccesses,avgrelaccess);accesshierarchies=append(accesshierarchies,accesshierarchy)
      decays=append(decays,decay);ccities=append(ccities,city);years=append(years,year)
      
      # TODO add some maps -> access gains ? avgrelaccess ? (comparable scale !)
      
      
      }
  }
}

res=data.frame(accesshierarchy=accesshierarchies,avgrelaccess=avgrelaccesses,decay=decays,year=years,city=ccities)

g=ggplot(res,aes(x=year,y=accesshierarchy,colour=city,group=interaction(city,decay),linetype=as.character(decay),shape=as.character(decay)))
g+geom_point()+geom_line()+xlab("year")+ylab("Hierarchy of accessibility")+
  scale_shape_discrete(name="Decay")+scale_linetype_discrete(name="Decay")+
  scale_color_discrete(name="City")+
  stdtheme
ggsave(file=paste0(resdir,'test_accesshierarchy.png'),width=15,height=10,units='cm')
# TODO capitalize cities

g=ggplot(res,aes(x=year,y=avgrelaccess,colour=city,group=interaction(city,decay),linetype=as.character(decay),shape=as.character(decay)))
g+geom_point()+geom_line()+xlab("year")+ylab("Average relative accessibility")+stdtheme
ggsave(file=paste0(resdir,'test_avgaccess.png'),width=15,height=10,units='cm')


