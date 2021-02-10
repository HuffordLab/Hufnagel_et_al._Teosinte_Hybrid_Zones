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
dataHyb = data[,c(6,8,7)]  #pop, lon, lat
dataPM = data[,c(2,5,8,7)] #tax, hyb, lon, lat
dataP = subset(dataPM, V2=="Zea_mays_parviglumis" & V5!="hybrid")
dataM = subset(dataPM, V2=="Zea_mays_mexicana" & V5!="hybrid")
dataCP = subset(dataHyb, V6=="Central_Plateau")
dataSG = subset(dataHyb, V6=="South_Guerrero")
dataCB = subset(dataHyb, V6=="Central_Balsas")

#Seperate out groups to compare across environmental variables
dataP2 = dataP[,c(-1,-2)]
dataM2 = dataM[,c(-1,-2)]
dataCP2 = dataCP[,-1]
dataSG2 = dataSG[,-1]
dataCB2 = dataCB[,-1]

#import altitudes
dataPMAlt = data[,c(2,5,11)] #tax, hyb, alt
dataHybAlt = data[,c(6,11)]  #pop, alt
altP = subset(dataPMAlt, V2=="Zea_mays_parviglumis" & V5!="hybrid")
altM = subset(dataPMAlt, V2=="Zea_mays_mexicana" & V5!="hybrid")
altCP = subset(dataHybAlt, V6=="Central_Plateau")
altSG = subset(dataHybAlt, V6=="South_Guerrero")
altCB = subset(dataHybAlt, V6=="Central_Balsas")

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

#make box plots for altitude
titlesP = rep("Parviglumis",nrow(altP))
titlesM = rep("Mexicana",nrow(altM))
titlesCP = rep("Hyb CP",nrow(altCP))
titlesSG = rep("Hyb SG",nrow(altSG))
titlesCB = rep("Hyb CB",nrow(altCB))

altPboxData = matrix(c(titlesP,altP[,3]), ncol=2)
altMboxData = matrix(c(titlesM,altM[,3]), ncol=2)
altCPboxData = matrix(c(titlesCP,altCP[,2]), ncol=2)
altSGboxData = matrix(c(titlesSG,altSG[,2]), ncol=2)
altCBboxData = matrix(c(titlesCB,altCB[,2]), ncol=2)

altAllData = rbind(altPboxData, altMboxData, altCPboxData, altSGboxData, altCBboxData)
plot.new()
pdf("altitudePlotColored.pdf")
boxplot(as.numeric(altAllData[,2]) ~ altAllData[,1], main="Altitudes", ylab="Altitude (meters above sea level)", col=c("orange4","yellow1","wheat3","firebrick4","blue2"), xlab="")
dev.off()

#make box plots for annual precipitation
prePboxData = matrix(c(titlesP,precipP), ncol=2)
preMboxData = matrix(c(titlesM,precipM), ncol=2)
preCPboxData = matrix(c(titlesCP,precipCP), ncol=2)
preSGboxData = matrix(c(titlesSG,precipSG), ncol=2)
preCBboxData = matrix(c(titlesCB,precipCB), ncol=2)
preAllData = rbind(prePboxData, preMboxData, preCPboxData, preSGboxData, preCBboxData)
pdf("precipPlotColored.pdf")
boxplot(as.numeric(preAllData[,2]) ~ preAllData[,1], main="Annual Precipitation", ylab="Annual Precipitation (mm)", col=c("orange4","yellow1","wheat3","firebrick4","blue2"), xlab="")
dev.off()

#make box plots for annual mean temperature
tempPboxData = matrix(c(titlesP,tempP), ncol=2)
tempMboxData = matrix(c(titlesM,tempM), ncol=2)
tempCPboxData = matrix(c(titlesCP,tempCP), ncol=2)
tempSGboxData = matrix(c(titlesSG,tempSG), ncol=2)
tempCBboxData = matrix(c(titlesCB,tempCB), ncol=2)
tempAllData = rbind(tempPboxData, tempMboxData, tempCPboxData, tempSGboxData, tempCBboxData)
pdf("tempPlotColored.pdf")
boxplot(as.numeric(tempAllData[,2]) ~ tempAllData[,1], main="Annual Mean Temperature", ylab="Annual Mean Temperature (deg Celcius * 10)", col=c("orange4","yellow1","wheat3","firebrick4","blue2"), xlab="")
dev.off()

