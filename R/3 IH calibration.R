library(GeoLight)
library(TwGeos)
library(GeoLocTools)
library(SGAT)
library(readxl)
library(stringr)
library(zoo)
library(MASS)
library(raster)
source("R/FUNCTIONS.R")


load("output/meta data of individuals with logger files.RData")
load("output/all twilight transitions of all logger run through twilightEdit and midnight sun cleaned.RData")
load("output/all activity data of all logger.RData")

# run inhabitat calibration on all locations in beginnning of May 
# as all birds are assumed to be at the colony then
# -> also remove all twilights were bird is in salt water
ids <- unique(twl4.all$id)
twl5ih <- NULL
for(ii in 1:length(ids)){
  
  twl4 <- twl4.all[twl4.all$id == ids[ii],]
  twl4$doy <- as.numeric(strftime(twl4$Twilight,"%j"))
  act  <- t2.all  [t2.all$id == ids[ii],]
  actsum <- data.frame(sum=tapply(act$act, act$doy,sum)*3/86400)
  actsum$doy <- as.numeric(rownames(actsum))
  
  # inhabitat calibration
  twl4ih <- twl4[twl4$Twilight >= as.POSIXct(paste(max(twl4$year),"-04-25",sep="")) & 
                 twl4$Twilight <= as.POSIXct(paste(max(twl4$year),"-05-12",sep="")),]
  # twl4ih <- twl4[twl4$Twilight >= as.POSIXct(paste(max(twl4$year),"-04-29",sep="")) & 
  #                  twl4$Twilight <= as.POSIXct(paste(max(twl4$year),"-05-09",sep="")),]
  twl4ih <- twl4ih[twl4ih$doy %in% actsum$doy[actsum$sum==0],]
  if(nrow(twl4ih)>0){
    twl4ih$colony <- meta2$Colony[meta2$animal.id==twl4$animal.id[1]][1]
    #twl4ih <- twl4[twl4$Twilight >= as.POSIXct(paste(min(twl4$year),"-07-01",sep="")) & twl4$Twilight <= as.POSIXct(paste(min(twl4$year),"-08-20",sep="")),]
    if(is.null(twl5ih)) twl5ih <- twl4ih else twl5ih <- rbind(twl5ih,twl4ih)
  }
}
# only consider birds from Vikingavatn and disregard one year track as that bird came back late
twl6ih<- twl5ih[twl5ih$colony=="Vikingavatn" & twl5ih$id != "339162 3479001 2009",] 
twl6ih$nestrise.diff <- as.numeric(twl6ih$Twilight- twl6ih$nest.rise)/60
twl6ih$nestset.diff <- as.numeric(twl6ih$Twilight- twl6ih$nest.set)/60
twl6ih <- twl6ih[abs(twl6ih$nestrise.diff) < 30 | abs(twl6ih$nestset.diff) < 30,]
twl6ih <- twl6ih[!is.na(twl6ih$Twilight),]

png("figures/IH calibration.png",units = "cm",width=15,height=15,res = 500)
calib <- MythresholdCalibration(twilight = twl6ih$Twilight, 
                                rise = twl6ih$Rise, 
                                lon = -mean(meta2$nest.lon[meta2$Colony=="Vikingavatn"]), 
                                lat = mean(meta2$nest.lat[meta2$Colony=="Vikingavatn"]), 
                                method = "gamma")
dev.off()

save(calib,file="output/IH calibration.RData")




# individual calibration
for(id in unique(twl6ih$animal.id)){
  png(paste("figures/IH calibration ",id,".png",sep=""),units = "cm",width=15,height=15,res = 500)
  calib <- MythresholdCalibration(twilight = twl6ih$Twilight[twl6ih$animal.id==id], 
                                  rise = twl6ih$Rise[twl6ih$animal.id==id], 
                                  lon = -meta2$nest.lon[meta2$animal.id==id][1], 
                                  lat = meta2$nest.lat[meta2$animal.id==id][1], 
                                  method = "gamma")
  dev.off()
}
