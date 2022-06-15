library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales)   #for transparency
library(raster)

#Load shape file for states 
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
states=readShapePoly("mexstates.shp",proj4string=crswgs84,verbose=TRUE)

# Load tiles of elevation graphics
block1 <- getData('SRTM', lon=-105, lat=25)
block2 <- getData('SRTM', lon=-100, lat=20)
block3 <- getData('SRTM', lon=-95, lat=20)
block4 <- getData('SRTM', lon=-90, lat=20)
block5 <- getData('SRTM', lon=-110, lat=25)
block6 <- getData('SRTM', lon=-100, lat=25)
block7 <- getData('SRTM', lon=-95, lat=25)
block8 <- getData('SRTM', lon=-90, lat=25)
block9 <- getData('SRTM', lon=-120, lat=30)
block10 <- getData('SRTM', lon=-115, lat=30)
block11 <- getData('SRTM', lon=-105, lat=20)
block12 <- getData('SRTM', lon=-110, lat=30)
block13 <- getData('SRTM', lon=-100, lat=30)
block14 <- getData('SRTM', lon=-115, lat=35)
block15 <- getData('SRTM', lon=-110, lat=35)
block16 <- getData('SRTM', lon=-115, lat=25)
block17 <- getData('SRTM', lon=-105, lat=30)
block18 <- getData('SRTM', lon=-120, lat=35)
block19 <- getData('SRTM', lon=-105, lat=35)

#Mash tiles into one graphic
altitudes <- mosaic(block1, block2, block3, block4, block5, block6, block7, block8, block9, block10, block11, block12, block13, block14, block15, block16, block17, block18, block19, fun=mean)

#create my own color palette which is greyscale
greyscale = c("gray100", "gray99", "gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray87", "gray86", "gray85", "gray84", "gray83", "gray82", "gray81", "gray80", "gray79", "gray78", "gray77", "gray76", "gray75", "gray74", "gray73", "gray72", "gray71", "gray70", "gray69", "gray68", "gray67", "gray66", "gray65", "gray64", "gray63", "gray62", "gray61", "gray60", "gray59", "gray58", "gray57", "gray56", "gray55", "gray54", "gray53", "gray52", "gray51", "gray50", "gray49", "gray48", "gray47", "gray46", "gray45", "gray44", "gray43", "gray42", "gray41", "gray40", "gray39", "gray38", "gray37", "gray36", "gray35", "gray34", "gray33", "gray32", "gray31", "gray30", "gray29", "gray28", "gray27", "gray26", "gray25", "gray24", "gray23", "gray22", "gray21", "gray20", "gray19", "gray18", "gray17", "gray16", "gray15", "gray14", "gray13", "gray12", "gray11", "gray10", "gray9", "gray8", "gray7", "gray6", "gray5", "gray4", "gray3", "gray2", "gray1", "gray0")

data = read.table("55kTeoSNPsNew_v4.txt", header=FALSE, sep='\t')
data = data.frame(data) #convert data to a dataframe otherwise points gets mad

#set shapes and colors
cnt = 1
shapes = c()
colors = c()
for (rowNum in 1:nrow(data)){
  tax = data[rowNum,3]
  hyb = data[rowNum,10]
  if (tax == "Zea_mays_parviglumis"){
    shapes = c(shapes,15)
    if (hyb == "parvHC"){
      colors = c(colors, "blue2")
    }		
    else if (hyb == "hybrid"){
      colors = c(colors, "purple1")
    }
    else{
      colors = c(colors, "cyan") 
    }
  }
  
  else if (tax == "Zea_mays_mexicana"){
    shapes = c(shapes,17)
    if (hyb == "mexHC"){
      colors = c(colors, "firebrick4") 
    }		
    else if (hyb == "hybrid"){
      colors = c(colors, "purple1") 
    }
    else{ 
      colors = c(colors, "orangered") 
    }
  }
  
  else{
    print("ERROR1 !")	
  }
  
  cnt = cnt + 1
}

#Make plots
plot.new()
pdf("55kOnMexBig.pdf")
#plot Mexico
map("worldHires","Mexico", col="black", fill=FALSE, lwd=0.7, xlim=c(-117,-83))
#plot altitude
plot(altitudes, add=TRUE, col=greyscale)
#plot states
plot(states, lwd=0.5, add=TRUE)
#plot samples
points(as.numeric(as.character(data$V5)), as.numeric(as.character(data$V4)), pch=shapes, col=colors, cex=.7) 
dev.off()

pdf("55kOnMexSmall.pdf")
#plot Mexico
map("worldHires","Mexico", col="black", fill=FALSE, lwd=0.7, xlim=c(-107,-96), ylim=c(15,21))
#plot altitude
plot(altitudes, add=TRUE, col=greyscale)
#plot states
plot(states, lwd=0.5, add=TRUE)
#plot samples
points(as.numeric(as.character(data$V5)), as.numeric(as.character(data$V4)), pch=shapes, col=colors, cex=.75) 
#add legend
legend(list(x=-104,y=17), c("parvHC","parvAmb", "parvHyb","mexHC","mexAmb","mexHyb"), col=c("blue2","cyan","purple1","firebrick4","orangered","purple1"), pch=c(15,15,15,17,17,17), lwd=2, lty=0, y.intersp=1.3)
dev.off()