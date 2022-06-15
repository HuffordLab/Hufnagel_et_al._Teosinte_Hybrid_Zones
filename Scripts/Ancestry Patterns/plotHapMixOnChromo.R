#This script is designed to plot HAPMIX ancestry tracts using prepared 
#  input files.  
#Created by David E. Hufnagel on Jun 3, 2019

#import data files
#set directory here
cbTeloData = read.table("cbTeloHapmixProbs.txt")
cbR96Data = read.table("cbR96HapmixProbs.txt")
cbAlchoData = read.table("cbAlchoHapmixProbs.txt")
ebPasoData = read.table("ebPasoHapmixProbs.txt")
sgMochiData = read.table("sgMochiHapmixProbs.txt")
sgChilpData = read.table("sgChilpHapmixProbs.txt")


#Get centromere coords
#set working directory here
centCoordData = read.table("chosenCentCoords2.txt", skip=2)
mod = 30
c1Start = centCoordData[1, 2] - mod
c1Stop = centCoordData[1, 3] + mod
c2Start = centCoordData[2, 2] - mod
c2Stop = centCoordData[2, 3] + mod
c3Start = centCoordData[3, 2] - mod
c3Stop = centCoordData[3, 3] + mod
c4Start = centCoordData[4, 2] - mod
c4Stop = centCoordData[4, 3] + mod
c5Start = centCoordData[5, 2] - mod
c5Stop = centCoordData[5, 3] + mod
c6Start = centCoordData[6, 2] - mod
c6Stop = centCoordData[6, 3] + mod
c7Start = centCoordData[7, 2] - mod
c7Stop = centCoordData[7, 3] + mod
c8Start = centCoordData[8, 2] - mod
c8Stop = centCoordData[8, 3] + mod
c9Start = centCoordData[9, 2] - mod
c9Stop = centCoordData[9, 3] + mod
c10Start = centCoordData[10, 2] - mod
c10Stop = centCoordData[10, 3] + mod

#subset Data by chromosome
cbTelo1 = subset(cbTeloData, cbTeloData$V1 == 1)
cbTelo2 = subset(cbTeloData, cbTeloData$V1 == 2)
cbTelo3 = subset(cbTeloData, cbTeloData$V1 == 3)
cbTelo4 = subset(cbTeloData, cbTeloData$V1 == 4)
cbTelo5 = subset(cbTeloData, cbTeloData$V1 == 5)
cbTelo6 = subset(cbTeloData, cbTeloData$V1 == 6)
cbTelo7 = subset(cbTeloData, cbTeloData$V1 == 7)
cbTelo8 = subset(cbTeloData, cbTeloData$V1 == 8)
cbTelo9 = subset(cbTeloData, cbTeloData$V1 == 9)
cbTelo10 = subset(cbTeloData, cbTeloData$V1 == 10)

cbR961 = subset(cbR96Data, cbR96Data$V1 == 1)
cbR962 = subset(cbR96Data, cbR96Data$V1 == 2)
cbR963 = subset(cbR96Data, cbR96Data$V1 == 3)
cbR964 = subset(cbR96Data, cbR96Data$V1 == 4)
cbR965 = subset(cbR96Data, cbR96Data$V1 == 5)
cbR966 = subset(cbR96Data, cbR96Data$V1 == 6)
cbR967 = subset(cbR96Data, cbR96Data$V1 == 7)
cbR968 = subset(cbR96Data, cbR96Data$V1 == 8)
cbR969 = subset(cbR96Data, cbR96Data$V1 == 9)
cbR9610 = subset(cbR96Data, cbR96Data$V1 == 10)

cbAlcho1 = subset(cbAlchoData, cbAlchoData$V1 == 1)
cbAlcho2 = subset(cbAlchoData, cbAlchoData$V1 == 2)
cbAlcho3 = subset(cbAlchoData, cbAlchoData$V1 == 3)
cbAlcho4 = subset(cbAlchoData, cbAlchoData$V1 == 4)
cbAlcho5 = subset(cbAlchoData, cbAlchoData$V1 == 5)
cbAlcho6 = subset(cbAlchoData, cbAlchoData$V1 == 6)
cbAlcho7 = subset(cbAlchoData, cbAlchoData$V1 == 7)
cbAlcho8 = subset(cbAlchoData, cbAlchoData$V1 == 8)
cbAlcho9 = subset(cbAlchoData, cbAlchoData$V1 == 9)
cbAlcho10 = subset(cbAlchoData, cbAlchoData$V1 == 10)

ebPaso1 = subset(ebPasoData, ebPasoData$V1 == 1)
ebPaso2 = subset(ebPasoData, ebPasoData$V1 == 2)
ebPaso3 = subset(ebPasoData, ebPasoData$V1 == 3)
ebPaso4 = subset(ebPasoData, ebPasoData$V1 == 4)
ebPaso5 = subset(ebPasoData, ebPasoData$V1 == 5)
ebPaso6 = subset(ebPasoData, ebPasoData$V1 == 6)
ebPaso7 = subset(ebPasoData, ebPasoData$V1 == 7)
ebPaso8 = subset(ebPasoData, ebPasoData$V1 == 8)
ebPaso9 = subset(ebPasoData, ebPasoData$V1 == 9)
ebPaso10 = subset(ebPasoData, ebPasoData$V1 == 10)

sgMochi1 = subset(sgMochiData, sgMochiData$V1 == 1)
sgMochi2 = subset(sgMochiData, sgMochiData$V1 == 2)
sgMochi3 = subset(sgMochiData, sgMochiData$V1 == 3)
sgMochi4 = subset(sgMochiData, sgMochiData$V1 == 4)
sgMochi5 = subset(sgMochiData, sgMochiData$V1 == 5)
sgMochi6 = subset(sgMochiData, sgMochiData$V1 == 6)
sgMochi7 = subset(sgMochiData, sgMochiData$V1 == 7)
sgMochi8 = subset(sgMochiData, sgMochiData$V1 == 8)
sgMochi9 = subset(sgMochiData, sgMochiData$V1 == 9)
sgMochi10 = subset(sgMochiData, sgMochiData$V1 == 10)

sgChilp1 = subset(sgChilpData, sgChilpData$V1 == 1)
sgChilp2 = subset(sgChilpData, sgChilpData$V1 == 2)
sgChilp3 = subset(sgChilpData, sgChilpData$V1 == 3)
sgChilp4 = subset(sgChilpData, sgChilpData$V1 == 4)
sgChilp5 = subset(sgChilpData, sgChilpData$V1 == 5)
sgChilp6 = subset(sgChilpData, sgChilpData$V1 == 6)
sgChilp7 = subset(sgChilpData, sgChilpData$V1 == 7)
sgChilp8 = subset(sgChilpData, sgChilpData$V1 == 8)
sgChilp9 = subset(sgChilpData, sgChilpData$V1 == 9)
sgChilp10 = subset(sgChilpData, sgChilpData$V1 == 10)


#Divide coordinates by 1,000,000 to present them in terms of Mb
cbTelo1Mb = cbTelo1$V2/1000000;cbTelo2Mb = cbTelo2$V2/1000000
cbTelo3Mb = cbTelo3$V2/1000000;cbTelo4Mb = cbTelo4$V2/1000000
cbTelo5Mb = cbTelo5$V2/1000000;cbTelo6Mb = cbTelo6$V2/1000000
cbTelo7Mb = cbTelo7$V2/1000000;cbTelo8Mb = cbTelo8$V2/1000000
cbTelo9Mb = cbTelo9$V2/1000000;cbTelo10Mb = cbTelo10$V2/1000000

cbR961Mb = cbR961$V2/1000000;cbR962Mb = cbR962$V2/1000000
cbR963Mb = cbR963$V2/1000000;cbR964Mb = cbR964$V2/1000000
cbR965Mb = cbR965$V2/1000000;cbR966Mb = cbR966$V2/1000000
cbR967Mb = cbR967$V2/1000000;cbR968Mb = cbR968$V2/1000000
cbR969Mb = cbR969$V2/1000000;cbR9610Mb = cbR9610$V2/1000000

cbAlcho1Mb = cbAlcho1$V2/1000000;cbAlcho2Mb = cbAlcho2$V2/1000000
cbAlcho3Mb = cbAlcho3$V2/1000000;cbAlcho4Mb = cbAlcho4$V2/1000000
cbAlcho5Mb = cbAlcho5$V2/1000000;cbAlcho6Mb = cbAlcho6$V2/1000000
cbAlcho7Mb = cbAlcho7$V2/1000000;cbAlcho8Mb = cbAlcho8$V2/1000000
cbAlcho9Mb = cbAlcho9$V2/1000000;cbAlcho10Mb = cbAlcho10$V2/1000000

ebPaso1Mb = ebPaso1$V2/1000000;ebPaso2Mb = ebPaso2$V2/1000000
ebPaso3Mb = ebPaso3$V2/1000000;ebPaso4Mb = ebPaso4$V2/1000000
ebPaso5Mb = ebPaso5$V2/1000000;ebPaso6Mb = ebPaso6$V2/1000000
ebPaso7Mb = ebPaso7$V2/1000000;ebPaso8Mb = ebPaso8$V2/1000000
ebPaso9Mb = ebPaso9$V2/1000000;ebPaso10Mb = ebPaso10$V2/1000000

sgMochi1Mb = sgMochi1$V2/1000000;sgMochi2Mb = sgMochi2$V2/1000000
sgMochi3Mb = sgMochi3$V2/1000000;sgMochi4Mb = sgMochi4$V2/1000000
sgMochi5Mb = sgMochi5$V2/1000000;sgMochi6Mb = sgMochi6$V2/1000000
sgMochi7Mb = sgMochi7$V2/1000000;sgMochi8Mb = sgMochi8$V2/1000000
sgMochi9Mb = sgMochi9$V2/1000000;sgMochi10Mb = sgMochi10$V2/1000000

sgChilp1Mb = sgChilp1$V2/1000000;sgChilp2Mb = sgChilp2$V2/1000000
sgChilp3Mb = sgChilp3$V2/1000000;sgChilp4Mb = sgChilp4$V2/1000000
sgChilp5Mb = sgChilp5$V2/1000000;sgChilp6Mb = sgChilp6$V2/1000000
sgChilp7Mb = sgChilp7$V2/1000000;sgChilp8Mb = sgChilp8$V2/1000000
sgChilp9Mb = sgChilp9$V2/1000000;sgChilp10Mb = sgChilp10$V2/1000000


#Plot data 
options(scipen=999) #disables scientific notation

#plot data for chromosome 1
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch1.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo1$V3~cbTelo1Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 1", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo1Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c1Start, col="red") #centromere range start
abline(v=c1Stop, col="red") #centromere range stop

plot(cbR961$V3~cbR961Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR961Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c1Start, col="red") #centromere range start
abline(v=c1Stop, col="red") #centromere range stop

plot(cbAlcho1$V3~cbAlcho1Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho1Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c1Start, col="red") #centromere range start
abline(v=c1Stop, col="red") #centromere range stop

plot(ebPaso1$V3~ebPaso1Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso1Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c1Start, col="red") #centromere range start
abline(v=c1Stop, col="red") #centromere range stop

plot(sgMochi1$V3~sgMochi1Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi1Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c1Start, col="red") #centromere range start
abline(v=c1Stop, col="red") #centromere range stop

plot(sgChilp1$V3~sgChilp1Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp1Mb)*1.02, digits = -0.5), 10), padj=-0.8)
#end with an empty plot for margins
abline(1, 0, col="blue")
abline(v=c1Start, col="red") #centromere range start
abline(v=c1Stop, col="red") #centromere range stop
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 2
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch2.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo2$V3~cbTelo2Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 2", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo2Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c2Start, col="red") #centromere range start
abline(v=c2Stop, col="red") #centromere range stop

plot(cbR962$V3~cbR962Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR962Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c2Start, col="red") #centromere range start
abline(v=c2Stop, col="red") #centromere range stop

plot(cbAlcho2$V3~cbAlcho2Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho2Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c2Start, col="red") #centromere range start
abline(v=c2Stop, col="red") #centromere range stop

plot(ebPaso2$V3~ebPaso2Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso2Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c2Start, col="red") #centromere range start
abline(v=c2Stop, col="red") #centromere range stop

plot(sgMochi2$V3~sgMochi2Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi2Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c2Start, col="red") #centromere range start
abline(v=c2Stop, col="red") #centromere range stop

plot(sgChilp2$V3~sgChilp2Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp2Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c2Start, col="red") #centromere range start
abline(v=c2Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 3
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch3.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo3$V3~cbTelo3Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 3", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo3Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c3Start, col="red") #centromere range start
abline(v=c3Stop, col="red") #centromere range stop

plot(cbR963$V3~cbR963Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR963Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c3Start, col="red") #centromere range start
abline(v=c3Stop, col="red") #centromere range stop

plot(cbAlcho3$V3~cbAlcho3Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho3Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c3Start, col="red") #centromere range start
abline(v=c3Stop, col="red") #centromere range stop

plot(ebPaso3$V3~ebPaso3Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso3Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c3Start, col="red") #centromere range start
abline(v=c3Stop, col="red") #centromere range stop

plot(sgMochi3$V3~sgMochi3Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi3Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c3Start, col="red") #centromere range start
abline(v=c3Stop, col="red") #centromere range stop

plot(sgChilp3$V3~sgChilp3Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp3Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c3Start, col="red") #centromere range start
abline(v=c3Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 4
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch4.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo4$V3~cbTelo4Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 4", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo4Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c4Start, col="red") #centromere range start
abline(v=c4Stop, col="red") #centromere range stop

plot(cbR964$V3~cbR964Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR964Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c4Start, col="red") #centromere range start
abline(v=c4Stop, col="red") #centromere range stop

plot(cbAlcho4$V3~cbAlcho4Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho4Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c4Start, col="red") #centromere range start
abline(v=c4Stop, col="red") #centromere range stop

plot(ebPaso4$V3~ebPaso4Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso4Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c4Start, col="red") #centromere range start
abline(v=c4Stop, col="red") #centromere range stop

plot(sgMochi4$V3~sgMochi4Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi4Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c4Start, col="red") #centromere range start
abline(v=c4Stop, col="red") #centromere range stop

plot(sgChilp4$V3~sgChilp4Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp4Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c4Start, col="red") #centromere range start
abline(v=c4Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 5
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch5.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo5$V3~cbTelo5Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 5", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo5Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c5Start, col="red") #centromere range start
abline(v=c5Stop, col="red") #centromere range stop

plot(cbR965$V3~cbR965Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR965Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c5Start, col="red") #centromere range start
abline(v=c5Stop, col="red") #centromere range stop

plot(cbAlcho5$V3~cbAlcho5Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho5Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c5Start, col="red") #centromere range start
abline(v=c5Stop, col="red") #centromere range stop

plot(ebPaso5$V3~ebPaso5Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso5Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c5Start, col="red") #centromere range start
abline(v=c5Stop, col="red") #centromere range stop

plot(sgMochi5$V3~sgMochi5Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi5Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c5Start, col="red") #centromere range start
abline(v=c5Stop, col="red") #centromere range stop

plot(sgChilp5$V3~sgChilp5Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp5Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c5Start, col="red") #centromere range start
abline(v=c5Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 6
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch6.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo6$V3~cbTelo6Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 6", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo6Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c6Start, col="red") #centromere range start
abline(v=c6Stop, col="red") #centromere range stop

plot(cbR966$V3~cbR966Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR966Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c6Start, col="red") #centromere range start
abline(v=c6Stop, col="red") #centromere range stop

plot(cbAlcho6$V3~cbAlcho6Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho6Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c6Start, col="red") #centromere range start
abline(v=c6Stop, col="red") #centromere range stop

plot(ebPaso6$V3~ebPaso6Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso6Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c6Start, col="red") #centromere range start
abline(v=c6Stop, col="red") #centromere range stop

plot(sgMochi6$V3~sgMochi6Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi6Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c6Start, col="red") #centromere range start
abline(v=c6Stop, col="red") #centromere range stop

plot(sgChilp6$V3~sgChilp6Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp6Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c6Start, col="red") #centromere range start
abline(v=c6Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 7
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch7.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo7$V3~cbTelo7Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 7", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo7Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c7Start, col="red") #centromere range start
abline(v=c7Stop, col="red") #centromere range stop

plot(cbR967$V3~cbR967Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR967Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c7Start, col="red") #centromere range start
abline(v=c7Stop, col="red") #centromere range stop

plot(cbAlcho7$V3~cbAlcho7Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho7Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c7Start, col="red") #centromere range start
abline(v=c7Stop, col="red") #centromere range stop

plot(ebPaso7$V3~ebPaso7Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso7Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c7Start, col="red") #centromere range start
abline(v=c7Stop, col="red") #centromere range stop

plot(sgMochi7$V3~sgMochi7Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi7Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c7Start, col="red") #centromere range start
abline(v=c7Stop, col="red") #centromere range stop

plot(sgChilp7$V3~sgChilp7Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp7Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c7Start, col="red") #centromere range start
abline(v=c7Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 8
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch8.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo8$V3~cbTelo8Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 8", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo8Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c8Start, col="red") #centromere range start
abline(v=c8Stop, col="red") #centromere range stop

plot(cbR968$V3~cbR968Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR968Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c8Start, col="red") #centromere range start
abline(v=c8Stop, col="red") #centromere range stop

plot(cbAlcho8$V3~cbAlcho8Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho8Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c8Start, col="red") #centromere range start
abline(v=c8Stop, col="red") #centromere range stop

plot(ebPaso8$V3~ebPaso8Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso8Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c8Start, col="red") #centromere range start
abline(v=c8Stop, col="red") #centromere range stop

plot(sgMochi8$V3~sgMochi8Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi8Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c8Start, col="red") #centromere range start
abline(v=c8Stop, col="red") #centromere range stop

plot(sgChilp8$V3~sgChilp8Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp8Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c8Start, col="red") #centromere range start
abline(v=c8Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 9
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch9.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo9$V3~cbTelo9Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 9", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo9Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c9Start, col="red") #centromere range start
abline(v=c9Stop, col="red") #centromere range stop

plot(cbR969$V3~cbR969Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR969Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c9Start, col="red") #centromere range start
abline(v=c9Stop, col="red") #centromere range stop

plot(cbAlcho9$V3~cbAlcho9Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho9Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c9Start, col="red") #centromere range start
abline(v=c9Stop, col="red") #centromere range stop

plot(ebPaso9$V3~ebPaso9Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso9Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c9Start, col="red") #centromere range start
abline(v=c9Stop, col="red") #centromere range stop

plot(sgMochi9$V3~sgMochi9Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi9Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c9Start, col="red") #centromere range start
abline(v=c9Stop, col="red") #centromere range stop

plot(sgChilp9$V3~sgChilp9Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp9Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c9Start, col="red") #centromere range start
abline(v=c9Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()


#plot data for chromosome 10
dev.new(width=13, height=8)
pdf("HAPMIXplots_Ch10.pdf")
par(mfrow=c(8,1), mar=c(.8,3.5,1,1), mgp=c(2,1,0))

#start with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(cbTelo10$V3~cbTelo10Mb, ylim=c(0,2), type="l", ylab="CBG_Telo", xlab="", xaxt="n")
title("Ancestry Across Chromosome 10", cex.main=1.8, line=-6.9, outer=TRUE)
axis(side=1, at=seq(0, round(max(cbTelo10Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c10Start, col="red") #centromere range start
abline(v=c10Stop, col="red") #centromere range stop

plot(cbR9610$V3~cbR9610Mb, ylab="CBG_Ahua", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbR9610Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c10Start, col="red") #centromere range start
abline(v=c10Stop, col="red") #centromere range stop

plot(cbAlcho10$V3~cbAlcho10Mb, ylab="CBG_Alcho", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(cbAlcho10Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c10Start, col="red") #centromere range start
abline(v=c10Stop, col="red") #centromere range stop

plot(ebPaso10$V3~ebPaso10Mb, ylab="EBG_PasoM", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(ebPaso10Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c10Start, col="red") #centromere range start
abline(v=c10Stop, col="red") #centromere range stop

plot(sgMochi10$V3~sgMochi10Mb, ylab="SGG_Mochi", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgMochi10Mb)*1.02, digits = -1), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c10Start, col="red") #centromere range start
abline(v=c10Stop, col="red") #centromere range stop

plot(sgChilp10$V3~sgChilp10Mb, ylab="SGG_Chilp", ylim=c(0,2), type="l", xlab="", xaxt="n")
axis(side=1, at=seq(0, round(max(sgChilp10Mb)*1.02, digits = -0.5), 10), padj=-0.8)
abline(1, 0, col="blue")
abline(v=c10Start, col="red") #centromere range start
abline(v=c10Stop, col="red") #centromere range stop
#end with an empty plot for margins
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Physical SNP Coordinates (Mb)", cex.lab=1.4, line=-4.9)
dev.off()

