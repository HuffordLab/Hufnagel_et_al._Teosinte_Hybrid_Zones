#used https://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf 
#for guidance. This script makes boxplots for precipitation and 
#temperature for parv, mex, and hybrid populations in the GAH project
#Created by David E. Hufnagel on Jun 20, 2019
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales)   #for transparency
library(raster)

#import coord data
data = read.table("55kTeoSNPsNew_v4.txt", header=FALSE, sep='\t')
dataHyb = data[,c(2,5,4,6)]  #pop, lon, lat, alt
dataPM = data[,c(3,10,5,4,6)] #tax, hyb, lon, lat, alt

#Seperate out groups to compare across environmental variables
dataP = subset(dataPM, V3=="Zea_mays_parviglumis" & V10!="hybrid")
dataM = subset(dataPM, V3=="Zea_mays_mexicana" & V10!="hybrid")
dataTelo = subset(dataHyb, V2=="Teloloapan")
dataAlcho = subset(dataHyb, V2=="Alcholoa")
dataAhua = subset(dataHyb, V2=="Ahuacatitlan")
dataPasoM = subset(dataHyb, V2=="Paso_Morelos")
dataChilp = subset(dataHyb, V2=="Chilpancingo")
dataMoch = subset(dataHyb, V2=="Mochitlan")

dataP2 = dataP[,c(3,4)]
dataM2 = dataM[,c(3,4)]
dataTelo2 = dataTelo[,c(2,3)]
dataAlcho2 = dataAlcho[,c(2,3)]
dataAhua2 = dataAhua[,c(2,3)]
dataPasoM2 = dataPasoM[,c(2,3)]
dataChilp2 = dataChilp[,c(2,3)]
dataMoch2 = dataMoch[,c(2,3)]

#extract environmental data
rastPrecip = raster("bio_22/bio12_22.bil")
rastTemp = raster("bio_22/bio1_22.bil")

altP = dataP[,5]
altM = dataM[,5]
altTelo = dataTelo[,4]
altAlcho = dataAlcho[,4]
altAhua = dataAhua[,4]
altPasoM = dataPasoM[,4]
altChilp = dataChilp[,4]
altMoch = dataMoch[,4]

precipP = extract(rastPrecip, dataP2)
precipM = extract(rastPrecip, dataM2)
precipTelo = extract(rastPrecip, dataTelo2)
precipAlcho = extract(rastPrecip, dataAlcho2)
precipAhua = extract(rastPrecip, dataAhua2)
precipPasoM = extract(rastPrecip, dataPasoM2)
precipChilp = extract(rastPrecip, dataChilp2)
precipMoch = extract(rastPrecip, dataMoch2)

tempP = extract(rastTemp, dataP2)
tempM = extract(rastTemp, dataM2)
tempTelo = extract(rastTemp, dataTelo2)
tempAlcho = extract(rastTemp, dataAlcho2)
tempAhua = extract(rastTemp, dataAhua2)
tempPasoM = extract(rastTemp, dataPasoM2)
tempChilp = extract(rastTemp, dataChilp2)
tempMoch = extract(rastTemp, dataMoch2)

#Set titles 
titles = c(rep("Parviglumis",length(tempP)), rep("Mexicana",length(tempM)), rep("CBG_Teloloapan",length(tempTelo)), rep("CBG_Alcholoa",length(tempAlcho)), rep("CBG_Ahuacatitlan",length(tempAhua)), rep("EBG_Paso_Morelos",length(tempPasoM)), rep("SGG_Chilpancingo",length(tempChilp)), rep("SGG_Mochitlan",length(tempMoch)))

#Plot for altitude
dev.new()
pdf("altBoxplotGAH.pdf")
alts = as.integer(c(altP, altM, altTelo, altAlcho, altAhua, altPasoM, altChilp, altMoch))
toPlot = data.frame(titles, alts)
titles_ordered = factor(toPlot$titles, c("Parviglumis", "Mexicana", "CBG_Teloloapan", "CBG_Alcholoa", "CBG_Ahuacatitlan", "EBG_Paso_Morelos", "SGG_Chilpancingo", "SGG_Mochitlan"))
par(mar=c(9,4.1,2.5,1), mgp=c(2.8,0.7,0))
boxplot(toPlot$alts ~ titles_ordered, las=2, col=c("blue2", "firebrick4", "darkorchid1", "darkorchid1", "darkorchid1", "darkorange2", "orchid1", "orchid1"), main="Altitude", ylab="Altitude (meters above sea level)", xlab="", cex.lab=1.5, cex.axis=1.0, cex.main=1.6)
dev.off()


#Plot for precipitation
dev.new()
pdf("precipBoxplotGAH.pdf")
precips=as.integer(c(precipP, precipM, precipTelo, precipAlcho, precipAhua, precipPasoM, precipChilp, precipMoch))
toPlot = data.frame(titles, precips)
titles_ordered = factor(toPlot$titles, c("Parviglumis", "Mexicana", "CBG_Teloloapan", "CBG_Alcholoa", "CBG_Ahuacatitlan", "EBG_Paso_Morelos", "SGG_Chilpancingo", "SGG_Mochitlan"))
par(mar=c(9,4.1,2.5,1), mgp=c(2.8,0.7,0))
boxplot(toPlot$precips ~ titles_ordered, las=2, col=c("blue2", "firebrick4", "darkorchid1", "darkorchid1", "darkorchid1", "darkorange2", "orchid1", "orchid1"), main="Annual Precipitation", ylab="Annual Precipitation (mm)", xlab="", cex.lab=1.5, cex.axis=1.0, cex.main=1.6)
dev.off()


#Plot for temperature
dev.new()
pdf("tempBoxplotGAH.pdf")
temps=as.integer(c(tempP, tempM, tempTelo, tempAlcho, tempAhua, tempPasoM, tempChilp, tempMoch))
toPlot = data.frame(titles, temps)
par(mar=c(9,4.1,2.5,1), mgp=c(2.6,0.7,0))
boxplot(toPlot$temps ~ titles_ordered, las=2, col=c("blue2", "firebrick4", "darkorchid1", "darkorchid1", "darkorchid1", "darkorange2", "orchid1", "orchid1"), main="Annual Mean Temperature", ylab="Annual Mean Temperature (deg Celcius * 10)", xlab="", cex.lab=1.5, cex.axis=1.0, cex.main=1.6)
dev.off()

