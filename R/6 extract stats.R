library(raster)
library(rptR)
source("R/FUNCTIONS.R")

load("output/all twilight transitions of all logger run through twilightEdit and midnight sun cleaned.RData")
load("output/all activity data of all logger.RData")
load("output/meta data of individuals with logger files.RData")
load("output/IH calibration.RData")
load("output/threshold locations all logger.RData")

table(twl4.all$id)
ids <- unique(twl4.all$id)


for(id in unique(twl4.all$id)){
  
  twl <- twl4.all[twl4.all$id == id,]
  act <- t2.all  [t2.all$id == id,]
  met <- meta2   [meta2$animal.id == twl$animal.id[1] & meta2$logger.id == twl$logger.id[1],][1,]
  loc <- path2   [path2$id == id,]
  loc$lon.diff   <- abs(loc$lon - -met$nest.lon)
  loc$s5lon.diff   <- abs(loc$s5lon - -met$nest.lon)
  loc$date <- as.Date(loc$datetime)
  loc$distance <- pointDistance(cbind(loc$s5lon, loc$s5lat),c(-met$nest.lon, met$nest.lat),lonlat = T)/1000
  loc$nsd <- (pointDistance(cbind(loc$s5lon, loc$s5lat),c(-met$nest.lon, met$nest.lat),lonlat = T)/1000)^2
  loc <- loc[loc$doy2 > 50,]
  loc <- loc[loc$doy2 < 320,]
  
  act2 <- data.frame(sum=tapply(act$act, as.Date(act$datetime), sum)*3/86400)
  act2$date <- as.Date(rownames(act2))
  actf <- act2[act2$sum>0.15,]
  
  # cL     <- changeLight(twl = twl.gl, quantile = 0.96, days = 20)
  # mS <- mergeSites(twl = twl.gl, site = cL$site, 
  #                   distThreshold = 500, 
  #                   degElevation = -calib[1] + 90,         # the HE corrected zero sun elevation angle
  #                   alpha = calib[3:4])                     # mask option
  
  #if(max(loc$s5lon.diff))
  date1.lon <- date2.lon <- date1.act <- date2.act <- NA
  date1.lon = loc$date[loc$s5lon.diff > mean(range(loc$s5lon.diff))][1]
  date2.lon = loc$date[loc$s5lon.diff > mean(range(loc$s5lon.diff))]
  date2.lon <- date2.lon[length(date2.lon)]+1
  
  
  
  
  if(nrow(loc)>60){
    
    # twl.gl <- export2GeoLight(twl) 
    # twl.gl$length  <- as.numeric(difftime(twl.gl$tSecond,twl.gl$tFirst,units="sec")) 
    # twl.gl$migrate <- 0
    # twl.gl$migrate[twl.gl$tSecond <= max(loc$datetime[loc$date==date1.lon]) + (0.8*86400) & 
    #                  twl.gl$tFirst  >= min(loc$datetime[loc$date==date1.lon]) - (0.8*86400)] <- 1
    # t2 <- twl.gl[twl.gl$tSecond <= max(loc$datetime[loc$date==date1.lon]) + (4*86400) & 
    #                twl.gl$tFirst  >= min(loc$datetime[loc$date==date1.lon]) - (4*86400),]
    # t2n <- t2[t2$migrate==0,]
    # act$migrate <- 0
    # act$migrate[act$datetime <= max(t2$tSecond[t2$migrate==1]) & act$datetime >= min(t2$tFirst[t2$migrate==1])] <- 1
    # a2 <- act[act$datetime <= max(t2$tSecond) & act$datetime >= min(t2$tFirst),]
    # 
    # plot(a2$datetime,a2$act,type='l',col=grey(0.8))
    # points(a2$datetime[a2$migrate==1],a2$act[a2$migrate==1],type='l',lwd=2)
    # 
    # fdm <- sum(a2$act[a2$datetime > t2$tFirst[t2$migrate==1 & t2$type==1] & 
    #                     a2$datetime < t2$tSecond[t2$migrate==1 & t2$type==1]])/
    #   (t2$length[t2$migrate==1 & t2$type==1]/3)
    # fnm <- sum(a2$act[a2$datetime > t2$tFirst[t2$migrate==1 & t2$type==2] & 
    #                     a2$datetime < t2$tSecond[t2$migrate==1 & t2$type==2]])/
    #   (t2$length[t2$migrate==1 & t2$type==2]/3)
    # 
    # for(tt in 1:length(t2n$tFirst[t2n$type==1])){
    #   t3 <- t2[t2$type==1,][tt,]
    #   
    #   fdn <- sum(a2$act[a2$datetime > t3$tFirst[t3$type==1] & a2$datetime < t3$tSecond[t3$type==1]])/(t3$length[t3$type==1]/3)
    #   if(tt==1) fdn2 <- fdn else fdn2 <- sum(fdn,fdn2)
    # }
    # fdn2 <- fdn2/tt
    # for(tt in 1:length(t2n$tFirst[t2n$type==2])){
    #   t3 <- t2[t2$type==2,][tt,]
    #   
    #   fnn <- sum(a2$act[a2$datetime > t3$tFirst[t3$type==2] & a2$datetime < t3$tSecond[t3$type==2]])/(t3$length[t3$type==2]/3)
    #   if(tt==1) fnn2 <- fnn else fnn2 <- sum(fnn,fnn2)
    # }
    # fnn2 <- fnn2/tt
    
    
    
    date1.act = as.Date(rownames(actf[1,]))
    date2.act = as.Date(rownames(actf[nrow(actf),]))+1
    
    act3l <- act2[act2$date >= date1.lon & act2$date <= date2.lon,]
    act3a <- act2[act2$date >= date1.act & act2$date <= date2.act,]
    
    da <- data.frame(date1.act, date2.act,
                     date1.lon, date2.lon,
                     activity = actf$sum[1], id, 
                     n.act=nrow(act2),
                     n.lon=nrow(loc),
                     prop.salt.lon = length(act3l$sum[act3l$sum==0])/length(act3l$sum),
                     prop.salt.act = length(act3a$sum[act3a$sum==0])/length(act3a$sum),
                     max.distance = max(loc$distance,na.rm=T)
                     # fracwet.mig.d <- fdm,
                     # fracwet.mig.n <- fnm,
                     # fracwet.nm.d <- fdn2,
                     # fracwet.nm.n <- fnn2
                     # dl.autumn =  (mean(t2$length[t2$type==1 & t2$migrate==1])-mean(t2$length[t2$type==1 & t2$migrate==0]))/mean(t2$length[t2$type==1 & t2$migrate==0]),
                     # nl.autumn =  (mean(t2$length[t2$type==2 & t2$migrate==1])-mean(t2$length[t2$type==2 & t2$migrate==0]))/mean(t2$length[t2$type==2 & t2$migrate==0]),
                     # dl.spring =  (mean(t3$length[t3$type==1 & t3$migrate==1])-mean(t3$length[t3$type==1 & t3$migrate==0]))/mean(t3$length[t3$type==1 & t3$migrate==0]),
                     # nl.spring =  (mean(t3$length[t3$type==2 & t3$migrate==1])-mean(t3$length[t3$type==2 & t3$migrate==0]))/mean(t3$length[t3$type==2 & t3$migrate==0]))
                     # dl.autumn =  (mean(t2$length[t2$type==1 & t2$migrate==1])-mean(t2$length[t2$type==1 & t2$migrate==0])),
                     # nl.autumn =  (mean(t2$length[t2$type==2 & t2$migrate==1])-mean(t2$length[t2$type==2 & t2$migrate==0])),
                     # dl.spring =  (mean(t3$length[t3$type==1 & t3$migrate==1])-mean(t3$length[t3$type==1 & t3$migrate==0])),
                     # nl.spring =  (mean(t3$length[t3$type==2 & t3$migrate==1])-mean(t3$length[t3$type==2 & t3$migrate==0]))
                     )
    
    if(id == unique(twl4.all$id)[1]) da2 <- da else da2 <- rbind(da2,da)
  }
}
da2$animal_id <- str_split_fixed(da2$id," ",2)[,1]
da2$doy1.act <- as.numeric(strftime(da2$date1.act,"%j"))
da2$doy1.lon <- as.numeric(strftime(da2$date1.lon,"%j"))
da2$doy2.act <- as.numeric(strftime(da2$date2.act,"%j"))
da2$doy2.lon <- as.numeric(strftime(da2$date2.lon,"%j"))
da2$doy1.diff <- da2$doy1.lon - da2$doy1.act
da2$doy2.diff <- da2$doy2.act - da2$doy2.lon
da2$date2.lon[da2$doy1.diff < 0] <- da2$doy2.lon[da2$doy1.diff < 0] <- NA
da2$date1.lon[da2$doy1.diff < 0] <- da2$doy1.lon[da2$doy1.diff < 0] <- NA
da2$date2.lon[da2$id=="325233 V185005001 2012"] <- da2$doy2.lon[da2$id=="325233 V185005001 2012"] <- NA
da2$date2.lon[da2$n.lon < 400] <- da2$doy2.lon[da2$n.lon < 400] <- NA
da2$date2.act[da2$n.act < 300] <- da2$doy2.act[da2$n.act < 300] <- NA
da2$y <- 1
da2$winter.act <- as.numeric(da2$date2.act) - as.numeric(da2$date1.act)
da2$winter.lon <- as.numeric(da2$date2.lon) - as.numeric(da2$date1.lon)

plot(da2$doy1.act,da2$doy1.lon)
lines(c(0,1000),c(0,1000))

cor(da2$prop.salt.act[!is.na(da2$prop.salt.lon)], 
    da2$prop.salt.lon[!is.na(da2$prop.salt.lon)])


rpt(max.distance ~ (1 | animal_id), grname = "animal_id", data = da2, 
    datatype = "Gaussian", nboot = 100, npermut = 100)

rpt(prop.salt.lon ~ (1 | animal_id), grname = "animal_id", data = da2, 
    datatype = "Gaussian", nboot = 100, npermut = 100)

rpt(winter.lon ~ (1 | animal_id), grname = "animal_id", data = da2, 
    datatype = "Gaussian", nboot = 100, npermut = 100)

plot(da2$doy1.lon,da2$dl.autumn)
abline(h=0,v=264,lty=3,col=grey(0.5))

plot(da2$doy2.lon,da2$dl.spring)
abline(h=0,v=80,lty=3,col=grey(0.5))

da2$dl.autumn[da2$doy1.lon < 264 & !is.na(da2$doy1.lon)] = (- da2$dl.autumn[da2$doy1.lon < 264 & !is.na(da2$doy1.lon)])
da2$nl.autumn[da2$doy1.lon < 264 & !is.na(da2$doy1.lon)] = (- da2$nl.autumn[da2$doy1.lon < 264 & !is.na(da2$doy1.lon)])

boxplot(da2$dl.autumn)
boxplot(da2$nl.autumn)



png("figures/day and night length difference during migration.png",units = "cm",width=15,height=15,res = 500)
opar <- par(mar=c(4.1,4.1,1,1),mfrow=c(1,1))
plot(da2$dl.autumn,da2$nl.autumn,asp=1,col=2,pch=1, ylim=c(-1.4,1.4),xlim=c(-1.4,1.4),
     ylab="night length difference during migratoin [hours]",xlab="day length difference during migration [hours]")
points(da2$dl.spring,da2$nl.spring,col=4,pch=1)
points(da2$dl.spring[da2$dl.spring>1.2],da2$nl.spring[da2$dl.spring>1.2],col=4,pch=4)
abline(h=0,v=0,lty=3,col=grey(0.5))
lines(c(-1000,1000),c(1000,-1000),lty=2,col=grey(0.5))
legend("topright",legend=c("autumn","spring"),col=c(2,4),pch=1)
da2b <- da2
da2b$dl.autumn <- da2b$dl.spring
da2b$nl.autumn <- da2b$nl.spring
da3 <- rbind(da2,da2b)
l1 <- lm(nl.autumn ~ dl.autumn, data=da3[da3$dl.autumn<1.2,])
text(1.1,-0.7,"R^2 = 0.58")
lines(range(da3$dl.autumn[da3$dl.autumn<1.2],na.rm=T),
      predict(l1,newdata=data.frame(dl.autumn = range(da3$dl.autumn[da3$dl.autumn<1.2],na.rm=T))))
par(opar)
dev.off()


png("figures/proportion of time in fresh water during winter.png",units = "cm",width=10,height=10,res = 500)
opar <- par(mar=c(2,2,2,1))
hist(da2$prop.salt.lon[!duplicated(da2$animal_id)],20,main="proportion in fresh water during winter")
par(opar)
dev.off()



ggplot(da2[!duplicated(da2$animal_id),], aes(x = prop.salt.lon)) +
  geom_dotplot(binaxis = "x", binwidth = 0.01, stackdir = "center")


png("figures/maximum distnace from nesting site.png",units = "cm",width=10,height=10,res = 500)
opar <- par(mar=c(2,2,2,1))
hist(da2$max.distance[!duplicated(da2$animal_id)],main="maximum distnace from nesting site")
par(opar)
dev.off()


png("figures/departure and arrival densities.png",units = "cm",width=15,height=10,res = 500)
dd1 <- density(da2$doy1.lon[!is.na(da2$doy1.lon)])
dd2 <- density(da2$doy2.lon[!is.na(da2$doy2.lon)])
dda1 <- density(da2$doy1.act[!is.na(da2$doy1.act)])
dda2 <- density(da2$doy2.act[!is.na(da2$doy2.act)])
m <- data.frame(dy=0:600,date=as.Date(0:600,"2010-01-01"))
m$year  <- as.numeric(strftime(m$date,"%Y"))
m$month <- as.numeric(strftime(m$date,"%m"))
m2 <- m[!duplicated(paste(m$month,m$year)),]
opar <- par(mfrow=c(1,1),mar=c(2,4,2,0))
plot(dd1$x,dd1$y,type="l",xaxt="n",ylab="density",
     xlim=c(150,550),ylim=c(-0.055,0.055))
abline(h=0)
lines(dd2$x+365,dd2$y)
lines(dda1$x,-dd2$y)
lines(dda2$x+365,-dd2$y)
axis(1,at=m2$dy,labels = strftime(m2$date,"%b"))
axis(3,at=m2$dy,labels = strftime(m2$date,"%b"))
mtext("saltwater switch",1,line=-1.5)
mtext("longitude",3,line=-1.5)
par(opar)
dev.off()


png("figures/winter period length 2.png",units = "cm",width=15,height=10,res = 500)
opar <- par(mfrow=c(1,1),mar=c(4,4,4,1))
boxplot(da2$winter.lon[!duplicated(da2$animal_id)],horizontal = T,xlim=c(0.5,2.5),ylim=c(150,270), xlab="length [days]")
boxplot(da2$winter.act[!duplicated(da2$animal_id)],horizontal = T, at=2,add=T)
axis(2, at=1:2, labels = c("longitude","saltwater switch"))
# axis(3, at=c(0,365*0.25, 365/2,365*0.75,365), labels = c(0,0.25,0.5,0.75,1))
par(opar)
dev.off()



#mode
dd1$x[which.max(dd1$y)]