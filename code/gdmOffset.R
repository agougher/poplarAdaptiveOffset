require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)

#function to remove intercept from gdm predictions
removeIntercept <- function(mod,pred){
  adjust <- 0 - log(1-pred) - mod$intercept
  adjustDissim <- 1-exp(0-adjust)
  return(adjustDissim)
}

##############
#Load and prep FST data
##############
adptMat <- read.csv("./data/fstTab.csv")

#Read in population locations
pops <- read.csv("./data/popInfo.csv")

#make sure the two are in the same order
all(adptMat$pop == pops$code)

##############
#Load and prep shapefile and climate data
##############
#load shapefile
shp <- shapefile("./popubals.shp")

#choose predictors
predNames <- c("bio2","bio3","bio10","bio11","bio18","bio19")

presClim <- stack(...) #stack current climate layers
presClim <- presClim[[predNames]]

#Creates pred data for gdm (cols = population name, long, lat, climate data)
pred <- data.frame(pop=pops$code,long=pops$long, lat=pops$lat, extract(presClim, y=pops[,c("long","lat")]),
                   stringsAsFactors=FALSE)

######################
#GDM model
######################
#Create site pair table
sitePair <- formatsitepair(bioDat = adptMat, bioFormat=3, siteColumn="pop", XColumn="long",YColumn="lat",predData=pred)

#Create and plot gdm
mod <- gdm(na.omit(sitePair), geo=FALSE)

#load future climate data
futClims <- stack(...) #stack future climate layers
futClims <- futClims[[predNames]]
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)

#Getting all coordinates in range map
popDat <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDat <- data.frame(distance=1, weight=1, popDat)
popDat <- split(popDat, seq(nrow(popDat)))


###############
#Forward offset calculation
##############
cl <- makeCluster(15) #ideally should be run in parallel to reduce computing time
registerDoParallel(cl)
forwardOffsetGDM <- foreach(i = 1:length(popDat), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population
  onePop <- popDat[[i]]
  
  #set up a dataframe where the first site is the focal population, and the second population
  #are sites across North America
  setUp <- cbind(onePop,futClimDat)
  colnames(setUp) <- c("distance","weights",
                       "s1.xCoord", "s1.yCoord",paste("s1.", predNames, sep=""), 
                       "s2.xCoord", "s2.yCoord",paste("s2.", predNames, sep=""))
  
  #rearrange the colums for the gdm prediction
  dat <- setUp[,c("distance","weights","s1.xCoord", "s1.yCoord","s2.xCoord", "s2.yCoord",
                  paste("s1.", predNames, sep=""), 
                  paste("s2.", predNames, sep=""))]
  
  #do the prediction and set up a dataframe with second sites x/y and predicted Fst
  combinedDat <- predict(object=mod, dat, time=FALSE)
  combinedDat <- data.frame(dat[,c("s2.xCoord","s2.yCoord")], predFst=removeIntercept(mod, combinedDat))
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coord <- onePop[,c("x","y")]
  
  #choose the pixels with the minimum fst
  minCoords <- combinedDat[which(combinedDat$predFst == min(combinedDat$predFst)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoords["dists"] <- distGeo(p1=coord, p2=minCoords[,1:2])
  minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoords <- minCoords[sample(1:nrow(minCoords),1),]
  
  #get local offset
  offset <- combinedDat[which(combinedDat$s2.xCoord == coord$x & combinedDat$s2.yCoord == coord$y),"predFst"]
  
  #get the minimum predicted fst - forward offset in this case
  minVal <- minCoords$predFst
  
  #get distance and coordinates of site that minimizes fst
  toGo <- minCoords$dists
  minPt <- minCoords[,c("s2.xCoord", "s2.yCoord")]
  
  #get bearing to the site that minimizes fst
  bear <- bearing(coord, minPt)
  
  #write out
  out <- c(x1=coord[[1]], y1=coord[[2]],local=offset,forwardFst=minVal, predDist=toGo, bearing=bear,x2=minPt[[1]],y2=minPt[[2]])
  
}

stopCluster(cl)

#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#forwardFst: forward offset
#predDist: distance to site of forward offset
#bearing: bearing to site of forward offset
#x2/y2: coordinate of site of forward offset
forwardOffsetGDM <- do.call(rbind, forwardOffsetGDM)

write.csv(forwardOffsetGDM,paste0("./forwardOffsetGDM.csv"), row.names=FALSE)


###############
#Reverse offset calculation
##############
#Getting all coordinates in the range in current climate
popDat <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDat <- data.frame(popDat)

#Gets climate data from the range in future climate
futClimMask <- mask(x=futClims, mask=shp)
futClimDat <- as.data.frame(futClimMask, xy=TRUE, na.rm=TRUE)

#set up for prediction
futClimDat <- data.frame(distance=1, weight=1, futClimDat)

###############
#Reverse offset calculation
##############
cl <- makeCluster(15)
registerDoParallel(cl)
reverseOffsetGDM <- foreach(i = 1:nrow(futClimDat), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population
  onePop <- futClimDat[i,]
  
  #set up a dataframe where the first site is the focal population, and the second population
  #are sites across the range
  setUp <- cbind(onePop,popDat)
  colnames(setUp) <- c("distance","weights",
                       "s1.xCoord", "s1.yCoord",paste("s1.", predNames, sep=""), 
                       "s2.xCoord", "s2.yCoord",paste("s2.", predNames, sep=""))
  
  #rearrange the colums for the gdm prediction
  dat <- setUp[,c("distance","weights","s1.xCoord", "s1.yCoord","s2.xCoord", "s2.yCoord",
                  paste("s1.", predNames, sep=""), 
                  paste("s2.", predNames, sep=""))]
  
  #do the prediction and set up a dataframe with second sites x/y and predicted Fst
  combinedDat <- predict(object=mod, dat, time=FALSE)
  combinedDat <- data.frame(dat[,c("s2.xCoord","s2.yCoord")], predFst=removeIntercept(mod, combinedDat))
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coord <- onePop[,c("x","y")]
  
  #choose the pixels with the minimum fst
  minCoords <- combinedDat[which(combinedDat$predFst == min(combinedDat$predFst)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoords["dists"] <- distGeo(p1=coord, p2=minCoords[,1:2])
  minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoords <- minCoords[sample(1:nrow(minCoords),1),]
  
  #get local offset
  offset <- combinedDat[which(combinedDat$s2.xCoord == coord$x & combinedDat$s2.yCoord == coord$y),"predFst"]
  
  #get the minimum predicted fst - reverse offset in this case
  minVal <- minCoords$predFst
  
  #get distance and coordinates of site that minimizes fst
  toGo <- minCoords$dists
  minPt <- minCoords[,c("s2.xCoord", "s2.yCoord")]
  
  #get bearing to the site that minimizes fst
  bear <- bearing(coord, minPt)
  
  #write out
  out <- c(x1=coord[[1]], y1=coord[[2]],local=offset,reverseFst=minVal, predDist=toGo, bearing=bear,x2=minPt[[1]],y2=minPt[[2]])
  
}

stopCluster(cl)

#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset - included as sanity check - should be identical to 'offset' from the calculation of forward offset above
#reverseFst: reverse offset
#predDist: distance to site of reverse offset
#bearing: bearing to site of reverse offset
#x2/y2: coordinate of site of reverse offset
reverseOffsetGDM <- do.call(rbind, reverseOffsetGDM)

write.csv(reverseOffsetGDM,paste0("./reverseOffsetGDM.csv"), row.names=FALSE)
