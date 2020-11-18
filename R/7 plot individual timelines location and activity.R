library(raster)
source("R/FUNCTIONS.R")

load("output/all twilight transitions of all logger run through twilightEdit and midnight sun cleaned.RData")
load("output/all activity data of all logger.RData")
load("output/meta data of individuals with logger files.RData")
load("output/IH calibration.RData")
load("output/threshold locations all logger.RData")

table(twl4.all$id)
ids <- unique(twl4.all$id)
id <- "339167 17647001 2010"
id <- '339185 17653001 2010'
id <- '339193 18B004001 2011'

for(id in unique(twl4.all$id)){
  
  twl <- twl4.all[twl4.all$id == id,]
  act <- t2.all  [t2.all$id == id,]
  loc <- path2   [path2$id == id,]
  met <- meta2   [meta2$animal.id == twl$animal.id[1] & meta2$logger.id == twl$logger.id[1],][1,]
  loc$nsd <- (pointDistance(cbind(loc$s5lon, loc$s5lat),c(-met$nest.lon, met$nest.lat),lonlat = T)/1000)^2
  act
  
  xlim <- range(loc$doy2)
  m <- data.frame(date=as.Date(0:364,"2010-01-01"))
  m$month <- as.numeric(strftime(m$date,"%m"))
  m$doy <- m$doy2 <- as.numeric(strftime(m$date,"%j"))
  m$doy2 <- m$doy2 - 182
  m$doy2[m$doy2<1] <- m$doy2[m$doy2<1] + 366
  m <- m[order(m$doy2),]
  m$calib <- 0
  m$calib[m$date <= as.Date("2010-05-12") & m$date >= as.Date("2010-04-25")] <- 1
  m2 <- m[!duplicated(m$month),]
  
  png(paste("figures/timeline ",id,".png",sep=""),units = "cm",width=15,height=20,res = 500)
  opar <- par(mfrow=c(4,1),mar=c(0,4,0,0),oma=c(2.5,0.2,1,1))
  plot(loc$doy2,loc$lat,type='o',pch=16,xlim=xlim,ylab="latitude",xaxt="n")
  abline(v=range(m$doy2[m$calib==1]),lty=3,col=3)
  plot(loc$doy2,loc$lon,type='o',pch=16,xlim=xlim,ylab="longitude",xaxt="n")
  abline(v=range(m$doy2[m$calib==1]),lty=3,col=3)
  plot(loc$doy2,loc$nsd,type='o',pch=16,xlim=xlim,ylab="NSD",ylim=c(0,max(loc$nsd,na.rm=T)),xaxt="n")
  abline(v=range(m$doy2[m$calib==1]),lty=3,col=3)
  plot(unique(act$doy2),tapply(act$act, act$doy2,sum)*3/86400,type='o',pch=16,xlim=xlim,ylab="activity",xaxt="n")
  abline(v=range(m$doy2[m$calib==1]),lty=3,col=3)
  axis(1,at=m2$doy2,labels = strftime(m2$date,"%b"))
  par(opar)
  dev.off()
}
