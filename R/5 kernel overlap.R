library(adehabitatHR)
library(raster)
library(ggplot2)
library(RColorBrewer)
library(rgdal)
library(stringr)
source("R/FUNCTIONS.R")

load("output/meta data of individuals with logger files.RData")
load("output/threshold locations all logger.RData")

proj.latlon <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj.aezd   <- CRS("+proj=aeqd  +lat_0=60  +lon_0=-9 +units=km")


land      <- readOGR("D:/map data","ne_10m_land")
grat5     <- readOGR("D:/map data","ne_50m_graticules_10")
ext  <- extent(-30,30,40,70)
ext2 <- extent(-878,1159,-1678,936)
land <- crop(land,ext)
land <- spTransform(land,proj.aezd)
grat5<- spTransform(grat5,proj.aezd)

bath <- raster("D:/map data/ETOPO1_Ice_c_geotiff.tif")
projection(bath) <- proj.latlon
bath <- crop(bath,ext)
bath2 <- projectRaster(bath,res=1,crs = (proj.aezd))
bath3 <- bath2
bath3[bath3 >  0]<-  0
bath3 <- crop(bath3,ext2)

batht <- bath2
batht[batht <  0]<-  NA
batht <- crop(batht,ext2)



loc2              <- path2[!is.na(path2$s5lat),]
coordinates(loc2) <- cbind(loc2$s5lon,loc2$s5lat)
proj4string(loc2) <- proj.latlon
loc2 <- spTransform(loc2,proj.aezd)
loc2 <- loc2[loc2$month %in% c(12,1),]

colony <- meta2[!duplicated(meta2$Colony),]
coordinates(colony) <- cbind(-colony$nest.lon,colony$nest.lat)
proj4string(colony) <- proj.latlon
colony <- spTransform(colony,proj.aezd)


# define kernel grid
test.grid  <- GridTopology(cellcentre.offset=c(-2000,-2000),cellsize=c(5,5),cells.dim=c(200*5,200*5))
test.point <- SpatialPoints(cbind(c(-100),c(-100)))
test.pixel <- SpatialPixels(test.point,proj4string = proj.aezd, round = NULL, grid = test.grid)
test.raster<- raster(extent(test.grid),crs=proj.aezd,resolution=25)


sex <- read.table("sexing.txt", header = T)


## run if you did not already
# ki   <- kernelUD(loc2[,12], grid=test.pixel, h='LSCV')
# for(u in 1:length(unique(loc2$id))){
#   ki1 <- data.frame(h = ki[[u]]@h$h, id = unique(loc2$id)[u])
#   if(u==1) ki2 <- ki1 else ki2 <- rbind(ki2,ki1)
# }
# summary(ki2$h)
# kj   <- kernelUD(loc2[,12], grid=test.pixel, h=17)
# 
# 
# k75  <- getverticeshr(kj,75)
# k75$kernel <- 75
# rownames(k75) <- NULL
# k50  <- getverticeshr(kj,50)
# k50$kernel <- 50
# rownames(k50) <- NULL
# k25  <- getverticeshr(kj,25)
# k25$kernel <- 25
# rownames(k25) <- NULL
# ki <- rbind(k75, k50, k25)
# ki$animal_id <- str_split_fixed(ki$id," ",2)[,1]
# 
# save(ki, file = "output/kernel vertices.RData")
# save(kj, file = "output/kernel UD.RData")
load("output/kernel vertices.RData")
load("output/kernel UD.RData")

ki$sex <- NA
for(i in 1:nrow(ki)) if(length(sex$sex[sex$id == ki$animal_id[i]])) ki$sex[i] <- sex$sex[sex$id == ki$animal_id[i]]

ol <- kerneloverlaphr(kj,method = "BA",percent = 75)
ol[upper.tri(ol,diag = T)] <- NA
rownames(ol) <- str_split_fixed(rownames(ol)," ",2)[,1]
colnames(ol) <- str_split_fixed(colnames(ol)," ",2)[,1]

ol2 <- ol
ol2[!is.na(ol2)] <- 0
for(u in 1:nrow(ol2)) ol2[colnames(ol2) == rownames(ol2)[u],u] <- 1
ol2[upper.tri(ol2,diag = T)] <- NA

ol3 <- data.frame(overlap=c(ol[ol2==1], ol[ol2==0]),
                  identical = c(rep(1,length(ol[ol2==1])), rep(0,length(ol[ol2==0]))))
ol3 <- ol3[!is.na(ol3$overlap),]
ol3$identical <- as.character(ol3$identical)



png("figures/winter overlap 2.png",units = "cm",width=10,height=15,res = 500)
opar <- par(mar=c(4,4,1,1))
boxplot(ol3$overlap~ol3$identical, ylim=c(0,1),xaxt="n",las=1,range=0,
        xlab="",ylab="Winter overlap as Bhattacharyya's affinity")
axis(1,at=1:2,labels = c("between","same"))
axis(1,at=1:2,labels = c("individuals","individual"),line=1,tick = F)
dev.off()







xex <- extent(bath3)[2]-extent(bath3)[1]
yex <- extent(bath3)[4]-extent(bath3)[3]
ratio <- yex/xex

bath.col  <- c(colorRampPalette(c(grey(0.1),rgb(80,90,100,maxColorValue = 255),
                                  rgb(190,240,250,maxColorValue = 255)))(100))

bath.col  <- c(colorRampPalette(c(grey(0.1), grey(0.95)))(100))
ter.col  <- c(colorRampPalette(c('beige', terrain.colors(10)[1]))(100))
ter.col  <- colorRampPalette(c('beige', "khaki4"))(100)
#ter.col <- terrain.colors(100,rev = T)
bird.col <- rainbow(19)#hcl.colors(19, palette = "viridis")
f.col <- colorRampPalette(brewer.pal(9, "OrRd")[3:9])(13)
m.col <- brewer.pal(8, "Blues")[3:8]


png("figures/winter 75 kernels 4.png",res=700,units="cm",width=20,height=20*1.27)
opar=par(mar=c(rep(1,4)))
image(bath3,asp=1,col=bath.col,ann=F,axes=F)
#plot(bath3,col=bath.col,maxpixels=1e10,interpolate=T,legend=F,add=T)
contour(-bath3,levels=c(500,200),add=T, col=grey(0.6),lty=1,lwd=0.5,labels=c("500 m","200 m"))
# plot(land,col=ter.col[1],border=grey(0.55),lwd=1,add=T)
# plot(batht,col=ter.col,add=T,legend=F)
# plot(land,col="transparent",border=grey(0.55),lwd=1,add=T)
plot(land,col=grey(1),border=grey(0.35),lwd=0.8,add=T)
plot(grat5,col=grey(0.8),lwd=1.3,add=T)
ids <- unique(ki$animal_id[ki$sex=="m"])
ids <- ids[!is.na(ids)]
for(id in 1:length(ids)){
  ki2 <- ki[ki$animal_id == sort(ids, decreasing = F)[id] & ki$kernel == 75,]
  for(i2 in 1:nrow(ki2)) plot(ki2[i2,],lwd=2,lty=i2,add=T,border=m.col[id],col=adjustcolor(m.col[id],0.3))
}
ids <- unique(ki$animal_id[ki$sex=="f"])
ids <- ids[!is.na(ids)]
for(id in 1:length(ids)){
  ki2 <- ki[ki$animal_id == sort(ids, decreasing = T)[id] & ki$kernel == 75,]
  for(i2 in 1:nrow(ki2)) plot(ki2[i2,],lwd=2,lty=i2,add=T,border=f.col[id],col=adjustcolor(f.col[id],0.3))
}
plot(colony, pch = 23, cex = 1.5, bg = "gold",add=T) # adding the release location

text(coordinates(colony)[,1],coordinates(colony)[,2],c("Víkingavatn", "Ástjorn"), 
     adj=c(0,1), pos=c(4,3), offset=1,cex=1.5)

axis(2,at = c(-1020, 100), tick = F, line = -1, cex.axis=0.7,
     labels = c(expression(paste("50",degree,"N")), 
                expression(paste("60",degree,"N"))))
axis(3,at = c(-470, -50, 380, 800), tick = F, line = -1, cex.axis=0.7,
     labels = c(expression(paste("20",degree,"W")), 
                expression(paste("10",degree,"W")), 
                expression(paste("0",degree,"")), 
                expression(paste("10",degree,"E"))))

box()
par(opar)
dev.off()

