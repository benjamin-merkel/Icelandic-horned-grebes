

load("output/meta data of individuals with logger files.RData")
fp <- "raw data"


# load activity data of all logger
for(ii in 1:nrow(meta2)){
  
  cat("\r",ii,"  ")
  t1          <- read.table(paste(fp,meta2$act.file[ii],sep="/"),skip=1,header = F,sep=",")
  t1$datetime <- as.POSIXct(strptime(t1$V2, "%d/%m/%y %H:%M:%S"),tz="UTC")
  t1$act    <- t1[,4]
  t1          <- t1[t1$datetime > meta2$deployed[ii] & t1$datetime < meta2$retrieved[ii],] 
  t2          <- t1[,c("datetime","act")]
  
  t2$year  <- as.numeric(strftime(t2$datetime,'%Y'))
  t2$month <- as.numeric(strftime(t2$datetime,'%m'))
  t2$day   <- as.numeric(strftime(t2$datetime,'%d'))
  t2$logger.id <- meta2$logger.id[ii]
  t2$animal.id <- meta2$animal.id[ii]
  
  if(ii==1) t2.all <- t2 else t2.all <- rbind(t2.all,t2)
}
t2.all$year2 <- t2.all$year
t2.all$year2[t2.all$month<7] <- t2.all$year2[t2.all$month<7] - 1
t2.all$id <- paste(t2.all$animal.id,t2.all$logger.id,t2.all$year2)
t2.all <- t2.all[t2.all$id %in% names(table(t2.all$id)[table(t2.all$id)>30]),]
t2.all$doy <- as.numeric(strftime(t2.all$datetime, "%j"))
t2.all$doy2 <- t2.all$doy
t2.all$doy2 <- t2.all$doy2 - 182
t2.all$doy2[t2.all$doy2<1] <- t2.all$doy2[t2.all$doy2<1] + 366
save(t2.all, file = "output/all activity data of all logger.RData")
load("output/all activity data of all logger.RData")
