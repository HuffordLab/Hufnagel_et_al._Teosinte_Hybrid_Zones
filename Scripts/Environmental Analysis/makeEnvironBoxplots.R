#used https://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf 
#for guidance. This script makes boxplots for altitude, precipitation, 
#and temperature for parv, mex, and hybrid groups
#Created by David E. Hufnagel
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales)   #for transparency
library(dismo)

#import coord data
#set working directory here
data = read.table("ZeaAllInfo.pmz", header=FALSE, sep='\t')
dataHyb = data[,c(6,8,7,11)]  #pop, lon, lat, alt
dataPM = data[,c(2,5,8,7,11)] #tax, hyb, lon, lat, alt
dataP = subset(dataPM, V2=="Zea_mays_parviglumis" & V5!="hybrid")
dataM = subset(dataPM, V2=="Zea_mays_mexicana" & V5!="hybrid")
dataCP = subset(dataHyb, V6=="Central_Plateau")
dataSG = subset(dataHyb, V6=="South_Guerrero")
dataCB = subset(dataHyb, V6=="Central_Balsas")

#Seperate out groups to compare across environmental variables
dataP2 = dataP[,c(3,4)]
dataM2 = dataM[,c(3,4)]
dataCP2 = dataCP[,c(2,3)]
dataSG2 = dataSG[,c(2,3)]
dataCB2 = dataCB[,c(2,3)]

#import altitudes
altP = dataP[,5]
altM = dataM[,5]
altCP = dataCP[,4]
altCB = dataCB[,4]
altSG = dataSG[,4]


#extract environmental data
rastTemp = raster("bio_22/bio1_22.bil")
rastPrecip = raster("bio_22/bio12_22.bil")

tempP = extract(rastTemp, dataP2)
tempM = extract(rastTemp, dataM2)
tempCP = extract(rastTemp, dataCP2)
tempSG = extract(rastTemp, dataSG2)
tempCB = extract(rastTemp, dataCB2)

precipP = extract(rastPrecip, dataP2)
precipM = extract(rastPrecip, dataM2)
precipCP = extract(rastPrecip, dataCP2)
precipSG = extract(rastPrecip, dataSG2)
precipCB = extract(rastPrecip, dataCB2)

#Set titles 
titles = c(rep("Parviglumis",length(tempP)), rep("Mexicana",length(tempM)), rep("Hyb CPG",length(tempCP)), rep("Hyb CBG",length(tempCB)), rep("Hyb SGG",length(tempSG)))


#make box plots for altitude
alts = as.integer(c(altP, altM, altCP, altCB, altSG))
dev.new()
pdf("altitudePlotColored.pdf")
toPlot = data.frame(titles, alts)
titles_ordered = factor(toPlot$titles, c("Parviglumis", "Mexicana", "Hyb CPG", "Hyb CBG", "Hyb SGG"))
par(mar=c(2,3.8,2.2,0.7), mgp=c(2.3,0.7,0))
boxplot(toPlot$alts ~ titles_ordered, main="Altitude", ylab="Altitude (meters above sea level)", col=c("blue2", "firebrick4", "yellow1", "orange4","wheat3"), xlab="", cex.lab=1.3, cex.axis=1.2, cex.main=1.6)
dev.off()

#make box plots for annual precipitation
dev.new()
pdf("precipPlotColored.pdf")
precips=as.integer(c(precipP, precipM, precipCP, precipCB, precipSG))
toPlot = data.frame(titles, precips)
par(mar=c(2,3.8,2.2,0.7), mgp=c(2.3,0.7,0))
boxplot(toPlot$precips ~ titles_ordered, main="Annual Precipitation", ylab="Annual Precipitation (mm)", col=c("blue2", "firebrick4", "yellow1", "orange4","wheat3"), xlab="", cex.lab=1.3, cex.axis=1.2, cex.main=1.6)
dev.off()

#make box plots for annual mean temperature
dev.new()
pdf("tempPlotColored.pdf")
temps=as.integer(c(tempP, tempM, tempCP, tempCB, tempSG))
toPlot = data.frame(titles, temps)
par(mar=c(2,3.8,2.2,0.7), mgp=c(2.3,0.7,0))
boxplot(toPlot$temps ~ titles_ordered, main="Annual Mean Temperature", ylab="Annual Mean Temperature (deg Celcius * 10)", col=c("blue2", "firebrick4", "yellow1", "orange4","wheat3"), xlab="", cex.lab=1.3, cex.axis=1.2, cex.main=1.6)
dev.off()

