library(Hmisc) #for adding minor tick marks
#set working directory here

#import data
data = read.table("ZeaAllInfo.pmz", sep="\t") 

#determine altitudes
alts = data$V11
parvAlts = sort(subset(alts, data$V2=="Zea_mays_parviglumis"), na.last=T)
mexAlts = sort(subset(alts, data$V2=="Zea_mays_mexicana"), na.last=T)

maizeData = subset(data, data$V2=="Zea_mays_mays")
mexMaizeData = subset(maizeData, maizeData$V9=="Mexico")
maizeAlts = sort(mexMaizeData$V11, na.last=T)

#plot altitudes
pdf("parvAltPlotOriginal.pdf")
plot(parvAlts, type="l", ylab="altitude (m)", ylim=c(0,2600))
minor.tick(ny=5)
dev.off()

pdf("mexAltPlotOriginal.pdf")
plot(mexAlts, type="l", ylab="altitude (m)", ylim=c(0,2600))
minor.tick(ny=5)
dev.off()

pdf("maizeAltPlotOriginal.pdf")
plot(maizeAlts, type="l", ylab="altitude (m)", ylim=c(0,2600))
minor.tick(ny=5)
dev.off()


