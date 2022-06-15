#Rather than plotting pies directly on the map, this version plots points on the
#map and pies separately, which will be connected to the map in GIMP.
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(mapplots) #for pie charts
library(scales)   #for transparency
library(raster)

#set working directory here

#Load shape file for states (Mostly from Kat's work)
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
states=readShapePoly("Mexstates/mexstates.shp",proj4string=crswgs84,verbose=TRUE)

hybs = as.matrix(read.table("ZeaPopInfo.txt", sep="\t"))
Hybs = data.frame(hybs)

plot.new()
#plot big
##plot background
pdf("HybridsBig.pdf")
map("worldHires","Mexico", col="black", fill=FALSE, lwd=0.7)
plot(states, lwd=0.7, add=TRUE)
#plot points
points(Hybs$V6, Hybs$V5, pch=16, cex=1.5)
dev.off()


#plot small
pdf("HybridsSmall.pdf")
map("worldHires","Mexico", col="black", fill=FALSE, lwd=0.7, xlim=c(-103,-82), ylim=c(15,22))
plot(states, lwd=0.6, add=TRUE)
#plot points
points(Hybs$V6, Hybs$V5, pch=16, cex=0.5)
#add legend
legend(list(x=-96,y=21.2), c("Parv attr","Mex attr","Maize attr"), col=c("blue2","firebrick4","gray38"), pch=c(19,19,19), lwd=2, lty=0, y.intersp=0.7)
dev.off()


#plot all pies
plot.new()
for (cnt in seq(1,nrow(Hybs))) {
  slices = c(as.numeric(as.character(strsplit(Hybs[,10],"_")[[cnt]][1])),as.numeric(as.character(strsplit(Hybs[,10],"_")[[cnt]][2])),as.numeric(as.character(strsplit(Hybs[,10],"_")[[cnt]][3])))
  pdfName = sprintf("pop%s_pie.pdf",trimws(Hybs$V1[cnt]))
  pdf(pdfName)
  pie(slices, labels=NA, main=NA, col=c("blue2","firebrick4","gray38"))
  title(trimws(Hybs$V1[cnt]), line=-27, cex.main=5.5)
  dev.off()
}

