require(raster)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)

##############
#Read in population locations & climate data
##############
pops <- read.csv("./data/popInfo.csv")

#choose predictors
predNames <- c("bio2","bio3","bio10","bio11","bio18","bio19")
presClim <- stack("C:/Users/Andy/Documents/GitHub/Ch2-novelDisappearing/currentClimate/current_northAmerica.grd")
presClim <- presClim[[predNames]]

#Creates pred data for gf (cols = population name, long, lat, climate data)
pred <- data.frame(pop=pops$code,long=pops$long, lat=pops$lat, extractNA(rast=presClim, x=pops$long, y=pops$lat),
                   stringsAsFactors=FALSE)

######################
#GF model
######################
#read in maf data
snps <- read.csv("./data/mafTab.csv", stringsAsFactors=FALSE)

#get snp names
snpsAll <- colnames(snps)[-c(1:3)]

#merge climate data and maf
snps <- merge(pred, snps, by.x="pop", by.y="code", all.x=TRUE)

#run gradient forest model
gfMod <- gradientForest(data=snps, predictor.vars=predNames,response.vars=snpsAll,
                        ntree=500, 
                        maxLevel=log2(0.368*nrow(snps)/2), trace=T, 
                        corr.threshold=0.70)

#read in shape file
shp <- shapefile("./popubals.shp")

#load future climate data
futClims <- stack(...) #stack future climate layers
futClims <- futClims[[predNames]]
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)

#transform future climate data with gradient forest model
futClimDatGF <- data.frame(futClimDat[,c("x","y")], predict(gfMod,futClimDat[,predNames])) 

#transform current climate data with gradient forest model
popDatGF <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDatGF <- data.frame(popDatGF[,c("x","y")], predict(gfMod, popDatGF[,predNames]))
popDatGF <- split(popDatGF, seq(nrow(popDatGF)))

###############
#Forward offset calculation
##############
cl <- makeCluster(15)
registerDoParallel(cl)
forwardOffsetGF <- foreach(i = 1:length(popDat), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population
  onePopGF <- popDatGF[[i]]
  
  #get destination populations and add gf distance
  combinedDatGF <- futClimDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], futClimDatGF[,predNames]))
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  
  #choose the pixels with the minimum fst
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dist == min(minCoordsGF$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  
  #get local offset
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y ==coordGF$y),"gfOffset"]
  
  #get the minimum predicted fst - forward offset in this case
  minValGF <- minCoordsGF$gfOffset
  
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF, forwardOffset=minValGF, predDist=toGoGF, bearing=bearGF,x2=minPtGF[[1]],y2=minPtGF[[2]])
  
}

stopCluster(cl)


#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#forwardOffset: forward offset
#predDist: distance to site of forward offset
#bearing: bearing to site of forward offset
#x2/y2: coordinate of site of forward offset
forwardOffsetGF <- do.call(rbind, forwardOffsetGF)

write.csv(forwardOffsetGF,paste0("./forwardOffsetGF.csv"), row.names=FALSE)


###############
#Reverse offset calculation
##############
#Getting all coordinates in the range in current climate and transform using gf
popDatGF <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDatGF <- data.frame(popDatGF[,c("x","y")], predict(gfMod, popDatGF[,predNames]))

#Gets climate data from the range in future climate and transform using gf
futClimMask <- mask(futClims, mask=shp)
futClimDat <- as.data.frame(futClimMask, xy=TRUE, na.rm=TRUE)
futClimDatGF <- data.frame(futClimDat[,c("x","y")],predict(gfMod,futClimDat[,predNames]))


###############
#Reverse offset calculation
##############
cl <- makeCluster(15)
registerDoParallel(cl)
reverseOffsetGF <- foreach(i = 1:nrow(futClimDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population in future climate
  onePopGF <- futClimDatGF[i,]
  
  #make prediction between focal population and current climate
  combinedDatGF <- popDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], popDatGF[,predNames]))
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  
  #choose the pixels with the minimum offset
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  
  #get local offset
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y == coordGF$y),"gfOffset"]
  
  #get the minimum predicted offset - reverse offset in this case
  minValGF <- minCoordsGF$gfOffset
  
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]],local=offsetGF, reverseOffset=minValGF, predDist=toGoGF, bearing=bearGF, x2=minPtGF[[1]],y2=minPtGF[[2]])
  
}

stopCluster(cl)


#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#reverseOffset: reverse offset
#predDist: distance to site of reverse offset
#bearing: bearing to site of reverse offset
#x2/y2: coordinate of site of reverse offset
reverseOffsetGF <- do.call(rbind, reverseOffsetGF)

write.csv(reverseOffsetGF,paste0("./reverseOffsetGF.csv"), row.names=FALSE)