library(GeoLight)
library(devtools)
devtools::install_github("SLisovski/TwGeos")
devtools::install_github("SLisovski/GeoLocTools")
devtools::install_github("SWotherspoon/SGAT")
library(TwGeos)
library(GeoLocTools)
library(SGAT)
library(readxl)
library(stringr)
library(zoo)
library(geosphere)
library(MASS)
source("R/FUNCTIONS.R")

meta               <- data.frame(read_excel("meta_Podicepsauritus_2009-2015.xls"))
meta$deployed      <- as.POSIXct(strptime(paste(meta$deployed,substr(meta$deploy_time,12,19)), "%Y-%m-%d %H:%M:%S"),tz="UTC")
meta$retrieved     <- as.POSIXct(strptime(paste(meta$recovered,substr(meta$time,12,19)), "%Y-%m-%d %H:%M:%S"),tz="UTC")
meta$calib1.start  <- as.POSIXct(strptime(paste(meta$calib_start1 ,substr(meta$time_start1 ,12,19)), "%Y-%m-%d %H:%M:%S"),tz="UTC")
meta$calib1.end    <- as.POSIXct(strptime(paste(meta$calib_end1  ,substr(meta$time_end1  ,12,19)), "%Y-%m-%d %H:%M:%S"),tz="UTC")
meta$calib3.start  <- as.POSIXct(strptime(paste(meta$calib_start3 ,substr(meta$time_start3  ,12,19)), "%Y-%m-%d %H:%M:%S"),tz="UTC")
meta$calib3.end    <- as.POSIXct(strptime(paste(meta$calib_end3  ,substr(meta$time_end3  ,12,19)), "%Y-%m-%d %H:%M:%S"),tz="UTC")


meta$retrieved[is.na(meta$retrieved)] <- as.POSIXct("2100-01-01")
meta$deployed[is.na(meta$deployed)]   <- as.POSIXct("2000-01-01")



fp <- "raw data"
light.files <- list.files(fp,".lig$$",recursive = T)
act.files   <- list.files(fp,".act$",recursive = T)
for(i in 1:length(light.files)){
  gls <- str_split_fixed(light.files[i],"_",2)[,1]
  meta$light.file[meta$logger.id == gls] <- light.files[i]
}
for(i in 1:length(act.files)){
  gls <- str_split_fixed(act.files[i],"_",2)[,1]
  meta$act.files[meta$logger.id == gls] <- act.files[i]
}
meta2 <- meta[!is.na(meta$light.file),]
meta2 <- meta2[meta2$logger.id != "V185006001",] # remove file without data

#save(meta2,file="output/meta data of individuals with logger files.RData")
load("output/meta data of individuals with logger files.RData")

offset <- 14 # adjusts the y-axis to put night (dark shades) in the middle
threshold = 11


###### calibration at known site - Víkingavatn (N66.10070 W16.84508)
# t1          <- read.table(paste(fp,meta2$light.file[3],sep="/"),skip=1,header = F,sep=",")
# t1$datetime <- as.POSIXct(strptime(t1$V2, "%d/%m/%y %H:%M:%S"),tz="UTC")
# t1$light    <- t1[,4]
# t1          <- t1[t1$datetime > meta2$calib3.start[ii] & t1$datetime < meta2$calib3.end[ii],] 
# t2          <- t1[,c("datetime","light")]
# colnames(t2)<- c("Date","Light")
# 
# lightImage( tagdata = t2, offset = offset, zlim = c(0, max(t2$Light)))
# tsimageDeploymentLines(t2$Date, lon = -meta2$nest.lon[ii], lat = meta2$nest.lat[ii],
#                        offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))
# 
# twl     <- twilightCalc(t1$datetime, t1$light, ask=F, LightThreshold = threshold)
# RT      <- getElevation(twl = twl, known.coord = c(-16.84508, 66.10070))
# RT

# calc all twilight transitions of all logger
for(ii in 1:nrow(meta2)){
 
  cat("\r",ii,"  ")
  t1          <- read.table(paste(fp,meta2$light.file[ii],sep="/"),skip=1,header = F,sep=",")
  # ERROR: last two days data invalid; 1 missing samples
  t1$datetime <- as.POSIXct(strptime(t1$V2, "%d/%m/%y %H:%M:%S"),tz="UTC")
  t1$light    <- t1[,4]
  t1          <- t1[t1$datetime > meta2$deployed[ii] & t1$datetime < meta2$retrieved[ii],] 
  t2          <- t1[,c("datetime","light")]
  colnames(t2)<- c("Date","Light")
  
  
  # lightImage( tagdata = t2, offset = offset, zlim = c(0, max(t2$Light)))
  # tsimageDeploymentLines(t2$Date, lon = -meta2$nest.lon[ii], lat = meta2$nest.lat[ii],
  #                        offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))
  
  twl         <- twilightCalc(t1$datetime, t1$light, ask=F, LightThreshold = threshold)
  twl$Rise    <- T
  twl$Rise[twl$type==2] <- F
  twl$year  <- as.numeric(strftime(twl$tFirst,'%Y'))
  twl$month <- as.numeric(strftime(twl$tFirst,'%m'))
  twl$diff  <- as.numeric(difftime(twl$tSecond, twl$tFirst, units = "hour"))
  twl$logger.id <- meta2$logger.id[ii]
  twl$animal.id <- meta2$animal.id[ii]
  
  if(ii==1) twl.all <- twl else twl.all <- rbind(twl.all,twl)
}
twl.all$year2 <- twl.all$year
twl.all$year2[twl.all$month<7] <- twl.all$year2[twl.all$month<7] - 1
twl.all$id <- paste(twl.all$animal.id,twl.all$logger.id,twl.all$year2)
twl.all <- twl.all[twl.all$id %in% names(table(twl.all$id)[table(twl.all$id)>30]),]
#save(twl.all, file = "output/all twilight transitions of all logger.RData")
load("output/all twilight transitions of all logger.RData")

# load all twilight transitions of all logger, transform to new twGeos format and run through twilightEdit
ids <- unique(twl.all$id)
for(ii in 1:length(ids)){
  
  cat("\r",ii,"  ")
  
  twl <- twl.all[twl.all$id == ids[ii],]
 
  start.cutoff <- twl$tSecond[twl$year==min(twl$year) & twl$diff> 48]
  start.cutoff <- start.cutoff[as.numeric(strftime(start.cutoff,"%m"))<8]
  end.cutoff   <- twl$tSecond[twl$year==max(twl$year) & twl$diff> 48]
  end.cutoff <- end.cutoff[!end.cutoff %in% start.cutoff]
  end.cutoff <- end.cutoff[as.numeric(strftime(end.cutoff,"%m"))>4]
  if(length(start.cutoff)==0) start.cutoff <- twl$tSecond[1]
  if(length(end.cutoff)==0) end.cutoff <- twl$tSecond[nrow(twl)]


  twl.switch <- twl
  twl.switch$tFirst <- twl.switch$tSecond
  twl.switch$Rise[twl.switch$type==2] <- T
  twl.switch$Rise[twl.switch$type==1] <- F

  twl <- rbind(twl,twl.switch)
  twl2 <- twl[twl$tSecond >= max(start.cutoff) & twl$tFirst <= min(end.cutoff),c(1,4)]
  colnames(twl2) <- c("Twilight","Rise")
  twl2 <- twl2[!duplicated(twl2$Twilight),]

  twl3 <- MytwilightEdit(twilights = twl2,
                      offset = offset,
                      window = 10,           # two days before and two days after
                      outlier.mins = 45,    # difference in mins
                      stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                      plot = TRUE)

  twl3$year  <- as.numeric(strftime(twl3$Twilight,'%Y'))
  twl3$month <- as.numeric(strftime(twl3$Twilight,'%m'))
  twl3$year2 <- twl3$year
  twl3$year2[twl3$month<7] <- twl3$year2[twl3$month<7] - 1
  
  twl3$logger.id <- twl$logger.id[1]
  twl3$animal.id <- twl$animal.id[1]
  
  if(ii==1) twl3.all <- twl3 else twl3.all <- rbind(twl3.all,twl3)
}
#save(twl3.all, file = "output/all twilight transitions of all logger run through twilightEdit.RData")
load("output/all twilight transitions of all logger run through twilightEdit.RData")


solarangle = -2.8 #HE ind 3
solarangle = -2.9 #HE ind 16
solarangle = -3.39 # RT ind 3
solarangle = -3.09 # IH ind 16

# calculate time of sunrise and sunset at nesting site throughout the year
twl3.all$nest.rise <- twilight(twl3.all$Twilight, lon = -mean(meta2$nest.lon, na.rm = T), lat =  mean(meta2$nest.lat, na.rm = T), zenith = 90 - solarangle, rise = T)
twl3.all$nest.set  <- twilight(twl3.all$Twilight, lon = -mean(meta2$nest.lon, na.rm = T), lat =  mean(meta2$nest.lat, na.rm = T), zenith = 90 - solarangle, rise = F)
twl3.all$nest.daylength <- as.numeric(difftime(twl3.all$nest.set, twl3.all$nest.rise, units = "hours"))

# remove all twilights removed by twilightEdit
twl4.all <- twl3.all[twl3.all$Deleted==F,]
twl4.all$id <- paste(twl4.all$animal.id,twl4.all$logger.id,twl4.all$year2)
# remove ids with less than 30 twilights
twl4.all <- twl4.all[twl4.all$id %in% names(table(twl4.all$id)[table(twl4.all$id)>30]),]
# remove twilights in May and July with NA or <0 day length at nesting site
twl4.all <- twl4.all[!(twl4.all$month==5 & (is.na(twl4.all$nest.daylength) | twl4.all$nest.daylength<=0)),]
twl4.all <- twl4.all[!(twl4.all$month==7 & (is.na(twl4.all$nest.daylength) | twl4.all$nest.daylength<=0)),]

#save(twl4.all, file = "output/all twilight transitions of all logger run through twilightEdit and midnight sun cleaned.RData")
load("output/all twilight transitions of all logger run through twilightEdit and midnight sun cleaned.RData")

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

