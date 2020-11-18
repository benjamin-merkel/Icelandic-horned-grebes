library(GeoLight)
library(TwGeos)
library(GeoLocTools)
library(SGAT)
library(readxl)
library(stringr)
library(zoo)
library(geosphere)
library(FLightR)
library(MASS)
source("R/FUNCTIONS.R")


load("output/all twilight transitions of all logger run through twilightEdit and midnight sun cleaned.RData")
load("output/meta data of individuals with logger files.RData")
load("output/IH calibration.RData")

Seasonal_palette <- grDevices::colorRampPalette(grDevices::hsv(1 - ((1:365) + (365/4))%%365/365, s = 0.8, v = 0.8), 
                                                space = "Lab")


# transform back to old GeoLight format
ids <- unique(twl4.all$id)
for(ii in 1:length(ids)){
  
  cat("\r",ii,"  ")
  
  twl <- twl4.all[twl4.all$id == ids[ii],]
  twlb <- export2GeoLight(twl)
  
  twlb$year  <- as.numeric(strftime(twlb$tFirst,'%Y'))
  twlb$month <- as.numeric(strftime(twlb$tFirst,'%m'))
  twlb$year2 <- twlb$year
  twlb$year2[twlb$month<7] <- twlb$year2[twlb$month<7] - 1
  
  twlb$logger.id <- twl$logger.id[1]
  twlb$animal.id <- twl$animal.id[1]
  
  if(ii==1) twl4 <- twlb else twl4 <- rbind(twl4,twlb)
}


tot <- 0.18#0.23#0.18
ids <- unique(twl4.all$id)
for(ii in 1:length(ids)){
  
  cat("\r",ii,"  ")
  
  twl <- twl4.all[twl4.all$id == ids[ii],]
  
  
  path <- thresholdLocation (twl$Twilight, twl$Rise, zenith = calib[1], tol=tot)
  path <- data.frame(lon=path[[2]][,1], lat=path[[2]][,2],datetime=path[[1]])
  
  path$month <- as.numeric(strftime(path$datetime,"%m"))
  path$year  <- as.numeric(strftime(path$datetime,'%Y'))
  path$year2 <- path$year
  path$year2[path$month<7] <- path$year2[path$month<7] - 1
  path$doy <- as.numeric(strftime(path$datetime,"%j"))
  path <- path[path$lon < 15 & path$lon > -26,]
  
  path$s5lon <- rollapply(path$lon, width = 5, FUN = mean, na.rm = TRUE, align = "center", partial = T)
  path$s5lat <- rollapply(path$lat, width = 5, FUN = mean, na.rm = TRUE, align = "center", partial = T)
  
  path$logger.id <- twl$logger.id[1]
  path$animal.id <- twl$animal.id[1]
  
  if(ii==1) path2 <- path else path2 <- rbind(path2,path)
}
path2$id <- paste(path2$animal.id,path2$logger.id,path2$year2)
path2$week <- as.numeric(strftime(path2$datetime,'%W'))
path2$week2 <- path2$week - 26
path2$week2[path2$week2 < 0] <- path2$week2[path2$week2 < 0] + 53
path2$doy <- path2$doy2 <- as.numeric(strftime(path2$datetime,"%j"))
path2$doy2 <- path2$doy2 - 182
path2$doy2[path2$doy2<1] <- path2$doy2[path2$doy2<1] + 366
path2$animal.id <- as.factor(str_split_fixed(path2$id," ",2)[,1])

#save(path2, file = "output/threshold locations all logger.RData")
load("output/threshold locations all logger.RData")


path2b <- path2[!is.na(path2$s5lat),]


# plot 5 location smoothed positions
cols <- colorRampPalette(c("orange",2,4,3))(366)
cols <- Seasonal_palette(366)
opar <- par(mfrow=c(1,1),mar=c(2,2,0,0))
plot(path2$s5lon, path2$s5lat, xlim = c(-26, 15), ylim = c(40, 70),asp=1.5,pch=19,col="white")
map("world",add=T,col=8)
for(i in ids) lines(path2b$s5lon[path2b$id==i],path2b$s5lat[path2b$id==i])
points(path2$s5lon, path2$s5lat ,pch=21, bg = cols[path2$doy])
map("world",add=T,col=8)
points(-meta2$nest.lon,meta2$nest.lat, pch = 23, cex = 1.5, bg = "gold") # adding the release location
mapplots::add.pie(x = -24, y = 45, z = rep(1, 12), radius = 4, col = Seasonal_palette(12), init.angle = 270)
par(opar)


#plot each individual
a.id <- unique(path2$animal.id)
for(a in a.id){
  
  png(paste("figures/",a,".png",sep=""),units = "cm",width=15,height=15,res = 500)
  
  ids3 <- unique(twl4.all$id[twl4.all$animal.id==a])
  cols <- Seasonal_palette(366)
  opar <- par(mfrow=c(1,1),mar=c(2,2,0,0))
  plot(path2$s5lon, path2$s5lat, xlim = c(-26, 15), ylim = c(40, 70),asp=1.5,pch=19,col="white")
  map("world",add=T,col=8)
  for(i in ids3) lines(path2$s5lon[path2$id==i],path2$s5lat[path2$id==i])
  for(i in ids3) lines(path2b$s5lon[path2b$id==i],path2b$s5lat[path2b$id==i],lty=2)
  points(path2$s5lon[path2$animal.id==a], path2$s5lat[path2$animal.id==a] ,pch=21+path2$year2-min(twl4.all$year2), bg = cols[path2$doy[path2$animal.id==a]])
  points(-meta2$nest.lon,meta2$nest.lat, pch = 23, cex = 1.5, bg = "gold") # adding the release location
  mapplots::add.pie(x = -21, y = 45, z = rep(1, 12), radius = 4, col = Seasonal_palette(12), init.angle = 270)
  par(opar)
  
  dev.off()
}


cols2 <- colorRampPalette(c(1,2,4,3,7,"white"))(length(levels(path2$animal.id)))
opar <- par(mfrow=c(2,1),mar=c(2,2,0,0))
plot(path2$doy2, path2$s5lat, ylim = c(40, 70),pch=19,col="white")
for(i in ids) lines(path2$doy2[path2$id==i],path2$s5lat[path2$id==i])
points(path2$doy2, path2$s5lat ,pch=21, bg = cols2[path2$animal.id])
abline(v=c(80,264))
plot(path2$doy2, path2$s5lon, ylim = c(-26, 15),pch=19,col="white")
for(i in ids) lines(path2$doy2[path2$id==i],path2$s5lon[path2$id==i])
points(path2$doy2, path2$s5lon ,pch=21, bg = cols2[path2$animal.id])
par(opar)



# calc weekly centroids
path3 <- data.frame(mean.lon = tapply(path2$s5lon,paste(path2$week2,path2$id),mean, na.rm=T),
                    mean.lat = tapply(path2$s5lat,paste(path2$week2,path2$id),mean, na.rm=T),
                    sd.lon   = tapply(path2$s5lon,paste(path2$week2,path2$id),sd, na.rm=T),
                    sd.lat   = tapply(path2$s5lat,paste(path2$week2,path2$id),sd, na.rm=T))
path3$id2 <- rownames(path3)
path3$week2 <- as.numeric(str_split_fixed(path3$id2," ",2)[,1])
path3$week <- path3$week2 + 26
path3$week[path3$week > 53] <- path3$week[path3$week > 53] - 53
path3$id <- str_split_fixed(path3$id2," ",2)[,2]
path3$animal.id <- as.factor(str_split_fixed(path3$id," ",2)[,1])
path3 <- path3[order(path3$week2),]
ids <- unique(path3$id)


cols2 <- colorRampPalette(c(1,2,4,3,7,"white"))(length(levels(path3$animal.id)))
opar <- par(mfrow=c(2,1),mar=c(2,2,0,0))
plot(path3$week2, path3$mean.lat, ylim = c(40, 70),pch=19,col="white")
for(i in ids) lines(path3$week2[path3$id==i],path3$mean.lat[path3$id==i])
points(path3$week2, path3$mean.lat ,pch=21, bg = cols2[path3$animal.id])
plot(path3$week2, path3$mean.lon, ylim = c(-26, 15),pch=19,col="white")
for(i in ids) lines(path3$week2[path3$id==i],path3$mean.lon[path3$id==i])
points(path3$week2, path3$mean.lon ,pch=21, bg = cols2[path3$animal.id])
par(opar)


path3 <- path3[!is.na(path3$mean.lat),]


# plot weekly centroids
cols <- colorRampPalette(c("orange",2,4,3))(50)
cols <- Seasonal_palette(365)
opar <- par(mfrow=c(1,1),mar=c(2,2,0,0))
plot(path3$mean.lon, path3$mean.lat, xlim = c(-26, 15), ylim = c(40, 70),asp=1.5,pch=19,col="white")
map("world",add=T,col=8)
points(path3$mean.lon, path3$mean.lat ,pch=21)
for(i in ids) lines(path3$mean.lon[path3$id==i],path3$mean.lat[path3$id==i])
points(path3$mean.lon, path3$mean.lat ,pch=21, bg = cols[path3$week*7])
map("world",add=T,col=8)
points(-meta2$nest.lon,meta2$nest.lat, pch = 21, cex = 1.5, bg = "gold") # adding the release location
mapplots::add.pie(x = -24, y = 45, z = rep(1, 12), radius = 4, col = Seasonal_palette(12), init.angle = 270)
par(opar)
