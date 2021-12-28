
## TC Network transportation time estimation

library(igraph)
library(rgdal)
library(rgeos)
library(cartography)
library(classInt)
library(RColorBrewer)



map<-function(data,layer,spdfid,dfid,variable,filename,title,legendtitle="",legendPosition="topleft",
              extent=NULL,withLayout=T,legendRnd=2,width=15,height=10,nclass=10,breaks=NULL,
              palette='Spectral',lwd=0.01,
              pngResolution = 1200,
              additionalLinelayers=NULL,additionalPointlayers=NULL,additionalLabelLayers=NULL,
              withScale=NULL){
  
  graphicsext = strsplit(filename,split='.',fixed=T)[[1]][2]
  if(graphicsext=='png'){png(file=filename,width=width,height=height,units='cm',res=pngResolution)}
  if(graphicsext=='pdf'){pdf(file=filename,width=width,height=height)}
  
  if(withLayout==T){par(mar = c(0.4,0.4,2,0.4))}else{par(mar = c(0.4,0.4,0.4,0.4))}
  
  #plot.new()
  
  layoutLayer(title = ifelse(withLayout,title,""), sources = "",
              author = "", col = ifelse(withLayout,"grey","white"), coltitle = "black", theme = NULL,
              bg = NULL, scale=withScale , frame = withLayout, north = F, south = FALSE,extent=extent)
  
  plot(layer, border = NA, col = "white",add=T)
  
  if(is.numeric(data[,variable])){
    cols <- carto.pal(pal1 = "green.pal",n1 = floor(nclass/2), pal2 = "red.pal",n2 = floor(nclass/2))
    if(palette=='Spectral'){cols = rev(brewer.pal(nclass,'Spectral'))}
    
    #breaks=classIntervals(data[,variable],nclass,style = 'equal')
    if(is.null(breaks)){breaks=classIntervals(data[,variable],nclass);breaks=breaks$brks}
    show(breaks)
    
    choroLayer(spdf = layer,spdfid = spdfid,
               df = data,dfid = dfid,
               var=variable,
               col=cols,colNA='lightgrey',breaks=breaks,
               add=TRUE,lwd = lwd,
               legend.pos = "n"
    )
  }
  
  if(is.factor(data[,variable])){
    typoLayer(spdf = layer,spdfid = spdfid,
              df = data,dfid = dfid,
              var=variable,
              #col=cols,colNA='lightgrey',breaks=breaks$brks,
              add=TRUE,lwd = lwd,
              legend.pos = "topleft"
    )
  }
  
  # additional layers
  if(!is.null(additionalLinelayers)){
    for(l in additionalLinelayers){
      if(class(l[[1]])=="SpatialPolygonsDataFrame"){
        plot(l[[1]], border = l[[2]],add=T)
      }else{
        plot(l[[1]], col = l[[2]],add=T)
      }
    }
  }
  
  if(!is.null(additionalPointlayers)){
    for(l in additionalPointlayers){
      plot(l[[1]], col = l[[2]],pch=20,size=1,add=T)
    }
  }
  
  if(!is.null(additionalLabelLayers)){
    for(l in additionalLabelLayers){
      labelLayer(l[[1]],col = l[[2]],txt=l[[3]])
    }
  }
  
  if(is.numeric(data[,variable])){
    legendChoro(pos =  legendPosition,title.txt = legendtitle,
                title.cex = 0.8, values.cex = 0.6, breaks, cols, cex = 0.7,
                values.rnd = legendRnd, nodata = TRUE, nodata.txt = "No data",
                nodata.col = 'lightgrey', frame = FALSE, symbol = "box"
    )
  } 
  #if(is.factor(data[,variable])){
  #  legendTypo(pos =  "bottomleft",title.txt = legendtitle,
  #              title.cex = 0.8, values.cex = 0.6, breaks$brks, cols, cex = 0.7,
  #              values.rnd = legendRnd, nodata = TRUE, nodata.txt = "No data",
  #              nodata.col = 'lightgrey', frame = FALSE, symbol = "box"
  #  )
  #}
  
  
  dev.off()
  
}




#'
#' @description Compute accessibilities
#' @param allmats : acc mat , or function of year giving it
computeAccess<-function(accessorigdata,accessdestdata=NULL,matfun,d0=1,mode="time"){
  if(is.null(dim(accessorigdata))){
    # simple data vector
    weightmat = exp(-matfun / d0)
    data = accessorigdata
    potod = which(names(data)%in%rownames(weightmat)&!is.na(data))
    if(mode=="time"){Pi=rep(1,length(potod));Pj=matrix(rep(1,length(potod)),nrow=length(potod))}
    if(mode=="weighteddest"){Pi=rep(1,length(potod));Pj=matrix(data[potod],nrow=length(potod))}
    if(mode=="weightedboth"){Pi=data[potod];Pj=matrix(data[potod],nrow=length(potod))}
    access = Pi*((weightmat[,names(data)[potod]]%*%Pj)[names(data)[potod],])
    return(access)
  }else{
    # data frame of yearly data
    accessyears = unique(accessorigdata$year)
    allaccess = data.frame()
    for(year in accessyears){
      if(is.matrix(matfun)){weightmat=matfun}else{weightmat=matfun(year)}
      yearlyaccessorig = accessorigdata[accessorigdata$year==year,]
      potorig = which(yearlyaccessorig$id%in%rownames(weightmat)&!is.na(yearlyaccessorig$var))
      potdest = which(accessdestdata$id%in%colnames(weightmat)&!is.na(accessdestdata$var))
      Pi = yearlyaccessorig$var[potorig];
      Ej = matrix(accessdestdata$var[potdest],nrow=length(potdest))
      access = Pi*(weightmat[,accessdestdata$id[potdest]]%*%Ej)[yearlyaccessorig$id[potorig],]/(sum(Pi))#*sum(Ej))*sum(weightmat[yearlyaccessorig$id[potorig],accessdestdata$id[potdest]]))
      allaccess=rbind(allaccess,data.frame(var = access,id=names(access),year=rep(year,length(access))))
      rm(weightmat);gc()
    }
    allaccess$id=as.character(allaccess$id);allaccess$year=as.character(allaccess$year)
    return(allaccess)
  }
}



#' 
#' 
addAdministrativeLayer<-function(g=empty_graph(0)$fun(0),
                                 admin_layer,
                                 connect_speed=1,
                                 attributes=list("CP"="INSEE_COMM"),
                                 empty_graph_heuristic=NULL
                                 ){
  if(is.character(admin_layer)){
    spath = strsplit(strsplit(admin_layer,'.shp')[[1]][1],'/')[[1]]
    admin <- readOGR(paste(spath[1:(length(spath)-1)],collapse="/"),spath[length(spath)])
  }else{
    admin = admin_layer
  }
  
  centroids = gCentroid(admin,byid = TRUE)
  attrvals = list()
  for(attr in names(attributes)){attrvals[[attr]]=as.character(admin@data[,attributes[[attr]]])}
  return(addPoints(g,centroids@coords,attrvals,list(speed=rep(connect_speed,length(admin))),empty_graph_heuristic))
}

#testadmin = addAdministrativeLayer(trgraph,"data/gis/communes.shp")


addPointsLayer<-function(g,points_layer,connect_speed=1){
  spath = strsplit(strsplit(points_layer,'.shp')[[1]][1],'/')[[1]]
  points <- readOGR(paste(spath[1:(length(spath)-1)],collapse="/"),spath[length(spath)])
  return(addPoints(g,points@coords,list(pointname=as.character(points$name)),list(speed=rep(connect_speed,length(points)))))
}


#' 
#' add vertices, connecting to closest stations
addPoints<-function(g,
                    coords,
                    v_attr_list,
                    e_attr_list,
                    empty_graph_heuristic=NULL
                    ){
  currentvid = ifelse(length(V(g))>0,max(as.numeric(V(g)$name)) + 1,1)
  attrs = v_attr_list
  attrs[["name"]] = as.character(currentvid:(currentvid+nrow(coords)-1))
  attrs[["station"]] = rep(FALSE,nrow(coords))
  attrs[["x"]] = coords[,1];attrs[["y"]] = coords[,2]
  currentg = add_vertices(graph = g,nv=nrow(coords),attr = attrs)
  
  # add connectors
  etoadd = c();elengths=c()
  for(k in 1:nrow(coords)){
    if(k%%1000==0){show(k)}
    if(length(which(V(currentg)$station==TRUE))>0){
      stationscoords = data.frame(id=V(currentg)$name[V(currentg)$station==TRUE],x=V(currentg)$x[V(currentg)$station==TRUE],y=V(currentg)$y[V(currentg)$station==TRUE])
      stationscoords$id=as.character(stationscoords$id)
      dists = sqrt((stationscoords$x - coords[k,1])^2 + (stationscoords$y - coords[k,2])^2)
      closest_station_id = stationscoords$id[which(dists==min(dists))[1]]
      etoadd = append(etoadd,c(attrs[["name"]][k],closest_station_id))
      elengths=append(elengths,min(dists))
    }
      
    if(empty_graph_heuristic=="nearest-neighbor"){
      xcoords=coords$x[-k];ycoords=coords$y[-k];ids=attrs[["name"]][-k]
      dists = sqrt((xcoords - coords[k,1])^2 + (ycoords - coords[k,2])^2)
      #neigh = c(xcoords[which(dists==min(dists))[1]],ycoords[which(dists==min(dists))[1]])
      etoadd = append(etoadd,c(attrs[["name"]][k],ids[which(dists==min(dists))[1]]))
      elengths=append(elengths,min(dists))
    }
    if(empty_graph_heuristic=="full"){
      for(l in 1:nrow(coords)){
         if(k!=l){
           etoadd = append(etoadd,c(attrs[["name"]][k],attrs[["name"]][l]))
           elengths = append(elengths,sqrt((coords[l,1] - coords[k,1])^2 + (coords[l,2] - coords[k,2])^2))
         }
      }
    }
  }
  attrs = e_attr_list
  attrs[["length"]] = elengths
  currentg = add_edges(graph=currentg,edges = etoadd,attr = attrs)
  return(currentg)
}

#testg = addPoints(trgraph,matrix(c(661796,6881580),nrow=1),list(),list(speed=c(0.5)))


#'
#' 
addFootLinks<-function(g,walking_speed=1,snap=100){
  distmat = dist(data.frame(x=V(g)$x[V(g)$station==TRUE],y=V(g)$y[V(g)$station==TRUE]))
  # TODO

}


#  sqrt(sum((stations@coords[stations$OBJECTID==210,]-stations@coords[stations$OBJECTID==467,])^2))
# -> 415.3246 between gare du nord et de l'est : (checked geoportail) : Lambert93 meters

#'
#' @name addTransportationLayer 
#' @description Construct tarnsportation graph by adding layers successively
#' 
#' @param 
#' 
addTransportationLayer<-function(stations_layer=NULL,
                                 link_layer,
                                 g=empty_graph(0)$fun(0),
                                 speed=1,
                                 snap=100,
                                 e_attr_names=NULL,
                                 reprojection=NULL
                                 ){
  show(paste0('Adding transportation network : stations = ',stations_layer,' ; links = ',link_layer))
  
  # construct vertex set
  vertexes = data.frame()
  if(length(V(g))>0){
    vertexes = data.frame(id=V(g)$name,x=V(g)$x,y=V(g)$y,station=V(g)$station)
    vertexes$id=as.numeric(as.character(vertexes$id))
    currentvid = vertexes$id[nrow(vertexes)] + 1 
  }
  
  if(!is.null(stations_layer)){
    if(is.character(stations_layer)){
      spath = strsplit(strsplit(stations_layer,'.shp')[[1]][1],'/')[[1]]
      stations <- readOGR(paste(spath[1:(length(spath)-1)],collapse="/"),spath[length(spath)])
    }else{stations <- stations_layer}
    
    
    if(!is.null(reprojection)){
      stations <- spTransform(stations, reprojection)
    }
    
    if(length(V(g))>0){
      coords=stations@coords
      for(i in 1:length(stations)){
        statdist = apply(vertexes[,c("x","y")] - matrix(rep(coords[i,],nrow(vertexes)),ncol=2,byrow=TRUE),1,function(r){sqrt(r[1]^2+r[2]^2)})
        # create only if does not exist
        #show(min(statdist))
        if(statdist[statdist==min(statdist)]>snap){
          vertexes=rbind(vertexes,c(id=currentvid,x=coords[i,1],y=coords[i,2],station=TRUE))
          currentvid=currentvid+1
        }
      }
    }else{
      vertexes=rbind(vertexes,data.frame(id=(nrow(vertexes)+1):(nrow(vertexes)+length(stations)),x=stations@coords[,1],y=stations@coords[,2],station=rep(TRUE,length(stations))))
      vertexes$id=as.numeric(as.character(vertexes$id))
    }
  }
  
  
  # links
  if(is.character(link_layer)){
    lpath = strsplit(strsplit(link_layer,'.shp')[[1]][1],'/')[[1]]
    links <- readOGR(paste(lpath[1:(length(lpath)-1)],collapse="/"),lpath[length(lpath)])
  }else{links <- link_layer}
  
  if(!is.null(reprojection)){
    links <- spTransform(links, reprojection)
  }
  
  edges = data.frame()
  #show(g)
  if(length(E(g))>0){
    edges = data.frame(from=tail_of(g,E(g))$name,to=head_of(g,E(g))$name,speed=E(g)$speed,length=E(g)$length)
    for(attrname in e_attr_names){edges=cbind(edges,get.edge.attribute(g,attrname));colnames(edges)[ncol(edges)]=attrname}
  }
  
  currentvid = ifelse(nrow(vertexes)>0,as.numeric(as.character(vertexes$id))[nrow(vertexes)] + 1,1)
  
  edges$from=as.character(edges$from);edges$to=as.character(edges$to)
  
  # convert shitty factor types
  for(j in 1:ncol(links@data)){links@data[,j]=as.numeric(as.character(links@data[,j]))}
  
  for(l in 1:length(links)){
    #show(l)
    #currentAdditionalAttrs=as.numeric(as.character(links@data[l,e_attr_names]))
    currentAdditionalAttrs=links@data[l,e_attr_names]
    #show(currentAdditionalAttrs)
    for(i in 1:length(links@lines[[l]]@Lines)){
      coords = links@lines[[l]]@Lines[[i]]@coords
      vids = c()
      #mincoords=apply(stations@coords,1,function(r){l=links@lines[[l]]@Lines[[i]]@coords;return(min(apply(abs(l-matrix(data=rep(r,nrow(l)),ncol=2,byrow = TRUE)),1,function(r){sqrt(r[1]^2+r[2]^2)})))})
      for(k in 1:nrow(coords)){
        if(nrow(vertexes)>0){
          statdist = apply(vertexes[,c("x","y")] - matrix(rep(coords[k,],nrow(vertexes)),ncol=2,byrow=TRUE),1,function(r){sqrt(r[1]^2+r[2]^2)})
        }else{statdist=c(2*snap)}
        if(statdist[statdist==min(statdist)]<snap){
          vids=append(vids,vertexes$id[statdist==min(statdist)])
          #show(paste0('existing : ',vids))
        }else{
          # else create new vertex
          vids=append(vids,currentvid)
          #show(paste0('new : ',vids))
          vertexes=rbind(vertexes,c(id=currentvid,x=coords[k,1],y=coords[k,2],station=FALSE))
          names(vertexes)<-c("id","x","y","station")
          currentvid=currentvid+1
        }
        #show(vertexes)
      }
      # add edges
      for(k in 2:nrow(coords)){
        addedge=c(from=vids[k-1],
                    to=vids[k],
                    speed=speed,
                    length=sqrt((coords[k-1,1]-coords[k,1])^2+(coords[k-1,2]-coords[k,2])^2),
                    currentAdditionalAttrs
                  )
        #show(addedge)
        edges=rbind(edges,addedge)
      }
    }
  }
  
  #show(edges)
  names(edges)<-c("from","to","speed","length",e_attr_names)
  
  
  comblist=list(speed="mean",length="sum")
  for(attrname in e_attr_names){comblist[[attrname]]="min"}
  res = simplify(graph_from_data_frame(edges,directed=FALSE,vertices = vertexes),edge.attr.comb = comblist)
  
  g = induced_subgraph(res,which(degree(res)>0))
  
  if(is.null(stations_layer)){
     V(g)$station = rep(TRUE,vcount(g))
  }
  
  return(g)
  
}

#g=addTransportationLayer('data/gis/gares.shp','data/gis/rer_lignes.shp')


# tests
#plot(stations)
#for(i in 1:length(links@lines[[1]]@Lines)){points(links@lines[[1]]@Lines[[i]]@coords,col='yellow')}
#for(i in 1:length(links@lines[[2]]@Lines)){points(links@lines[[2]]@Lines[[i]]@coords,col='blue')}
#for(i in 1:length(links@lines[[3]]@Lines)){points(links@lines[[3]]@Lines[[i]]@coords,col='green')}
#for(i in 1:length(links@lines[[4]]@Lines)){points(links@lines[[4]]@Lines[[i]]@coords,col='purple')}
#for(i in 1:length(links@lines[[4]]@Lines)){points(links@lines[[4]]@Lines[[i]]@coords,col='purple')}


#mincoords=apply(stations@coords,1,function(r){l=links@lines[[1]]@Lines[[2]]@coords;return(min(rowSums(abs(l-matrix(data=rep(r,nrow(l)),ncol=2,byrow = TRUE)))))})
#mincoords[mincoords<200]


#getDistMat<-function(g,fromids,toids){
#  fromnames=c();for(cp in iris$DCOMIRIS){fromids=append(fromids,which(V(g)$IRIS==cp));if(cp%in%V(g)$IRIS){fromnames=append(fromnames,cp)}}
#  tonames=c();for(cp in communes$INSEE_COMM){toids=append(toids,which(V(g)$CP==cp));if(cp%in%V(g)$CP){tonames=append(tonames,cp)}}
#  res = distances(graph = g,v = fromids,to = toids,weights = E(g)$speed*E(g)$length)
#  rownames(res)<-fromnames;colnames(res)<-tonames
#  return(res)
#}







