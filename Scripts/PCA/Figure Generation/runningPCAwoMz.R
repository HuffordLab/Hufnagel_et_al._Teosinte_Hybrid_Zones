#Parts were taken from Ranal.R written by Joost Van Heervaarden.
#Meant for making PCA plots without maize
#Created by David E. Hufnagel

#set working directory here
data = read.table("ZeaAllInfo.pmz.PCAmatrix.noMz", row.names=1, skip=1) #SNP file used to subset
markers = strsplit(strsplit(readLines("ZeaAllInfo.pmz", n=2)[2][1], "data):")[[1]][2], ",")[[1]]
info = read.table("ZeaAllInfo.pmz", sep="\t", row.names=1)              #big info file

#defining a function for all SNPs of a marker.  Homozygous high freq allele returns 2, heterozygous returns 1, homozygous alternate allele returns 0, any missing data returns NA.
fcdat<-function(x){
	x<-as.character(x)                 
	vec<-unlist(strsplit(x,split="/")) 
	avec<-vec[seq(1,length(vec),2)]    
	bvec<-vec[seq(2,length(vec),2)]    
	nucfreq<-sort(table(c(avec,bvec)),decreasing=TRUE) 
	nucfreq<-nucfreq[names(nucfreq)!="?"] 
	base1=names(nucfreq)[1]
	base2=names(nucfreq)[2]
	avec[avec=="?"]=NA 
	bvec[bvec=="?"]=NA
	avec[avec==base1&is.na(avec)==FALSE]<-"S" 
	avec[avec==base2&is.na(avec)==FALSE]<-"N" 
	bvec[bvec==base1&is.na(bvec)==FALSE]<-"S"
	bvec[bvec==base2&is.na(bvec)==FALSE]<-"N"
	avec[avec=="S"]<-1 
	avec[avec=="N"]<-0
	bvec[bvec=="S"]<-1
	bvec[bvec=="N"]<-0
	vec<-as.numeric(avec)+as.numeric(bvec)
	return(vec)
}

SubSampVars = function(inds, allVar){
	newVar = matrix(NA, nrow=length(inds), ncol=ncol(allVar))
	cnt = 1
	for(ind in inds){
		newVar[cnt,] = c(allVar[ind,])
		cnt = cnt + 1
	}
	return(newVar)
}

#apply the function and set the rownames back to the individual names
data2 = apply(data,2,fcdat)
rownames(data2) = rownames(data)

#make vectors of column means and variance factors, sqrt(p*(1-p))
mv = apply(data2,2,mean,na.rm=TRUE)
p = mv/2
vv = sqrt(p*(1-p))
#take out monomorphic loci
data2 = data2[,vv>0]
markers2 = markers[vv>0]

#Set colors to define populations by making a matrix of parv and mex names as rownames and colors in the only column.
spec = info$V2
hyb = info$V5
pop = info$V6
lat = info$V7
lon = info$V8
cntry = info$V9
names = rownames(info)
cols = c() #colors
shapes = c()
for(i in seq(1:nrow(info))){
	cc = "black" #color code
	sh = 11 #shape code
	if(spec[i] == "Zea_mays_parviglumis" && !is.na(spec[i])){ #Dont accept missing values as a match
		sh = 0
		if(hyb[i] == "parvHC"){
			cc = "blue2"
		}
		else if(hyb[i] == "other"){
			cc = "cyan"
		}
	}
	else if(spec[i] == "Zea_mays_mexicana" && !is.na(spec[i])){
		sh = 2
		if(hyb[i] == "mexHC"){
			cc = "firebrick4"
		}
		else if(hyb[i] == "other"){
			cc = "firebrick1"
		}
	}
	
	#If hybrid reset to a hybrid group value
	if((hyb[i] == "hybrid" && !is.na(hyb[i]))){ 
		sh = 8
		if(pop[i] == "Central_Plateau" && !is.na(pop[i])){
			cc = "yellow1"#cols = c(cols, "red") #CP hybrids
		}
		else if(pop[i] == "Central_Balsas" && !is.na(pop[i])){
			cc = "orange4"#cols = c(cols, "red4") #NB hybrids
		}
		else if(pop[i] == "South_Guerrero" && !is.na(pop[i])){
			cc = "wheat3"#cols = c(cols, "orange2") #SB hybrids
		}
		else if(lat[i] > 18.2 && lat[i] < 18.3 && lon[i] > -99.2 && lon[i] < -99.1 && !is.na(lat[i]) && !is.na(lon[i])){
			cc = "seagreen"#cols = c(cols, "violetred1") #Ahuacatitlan
		}
		else{ #non-grouped hybrids
			cc = "yellowgreen"
		}
	}
	else if((hyb[i] == "other" && !is.na(hyb[i]))){
		
	}
	
	#put cc and sh in list
	if(cc != "black"){
		cols = c(cols,cc)
		names(cols)[length(cols)] = names[i]
		
		shapes = c(shapes,as.numeric(sh))
		names(shapes)[length(shapes)] = names[i]
	}
}

#adding the colors and shapes to the data so they are all in 1 matrix
colM = matrix(cols, length(cols), 1)
shapesM = matrix(shapes, length(shapes), 1)
rownames(colM) = names(cols)
rownames(shapesM) = names(shapes)
data3 = matrix(NA, nrow=nrow(colM), ncol=ncol(data2)+2) #create an empty list of the proper size
cnt = 1
for(i in seq(1:nrow(data2))){
	dName = rownames(data2)[i]
	if(dName %in% rownames(colM)){
		data3[cnt,] = c(data2[i,], shapesM[dName,], colM[dName,])
		cnt = cnt + 1
	}
}
rownames(data3) = names(cols)

numDataCols = ncol(data3)-2
data3num = apply(data3[,1:numDataCols], 2, as.numeric, na.rm=TRUE)
mv = apply(data3num,2,mean,na.rm=TRUE) #ignore last column
p = mv/2
vv = sqrt(p*(1-p))

#Handle missing data
StNor = t((t(data3num)-mv)/vv) #Converting each marker's data to a standard normal distribution
StNor[is.na(StNor)]<-0 #Set missing data to the mean of the standard normal (0)

#Run PCA, keep any values above 1% total variance explained
newpc = prcomp(StNor)

################
# Plot results (% variance explaned was determined by examining newpc and hardcoded in after)
pdf("PC1vPC2woMz.pdf")
plot(newpc$x[,2]~newpc$x[,1],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]), ylab="PC 2   (2.10% var. explained)",xlab="PC 1   (6.48% var. explained)")
dev.off()

pdf("PC3vPC4woMz.pdf")
plot(newpc$x[,4]~newpc$x[,3],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]) ,ylab="PC 4   (1.08% var. explained)",xlab="PC 3   (1.28% var. explained)")
dev.off()

pdf("PC5vPC6woMz.pdf")
plot(newpc$x[,6]~newpc$x[,5],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]) ,ylab="PC 6   (0.91% var. explained)",xlab="PC 5   (0.98% var. explained)")
dev.off()

pdf("PC7vPC8woMz.pdf")
plot(newpc$x[,8]~newpc$x[,7],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]) ,ylab="PC 8   (0.76% var. explained)",xlab="PC 7   (0.86% var. explained)")
dev.off()

pdf("PC9vPC10woMz.pdf")
plot(newpc$x[,10]~newpc$x[,9],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]) ,ylab="PC 10   (0.69% var. explained)",xlab="PC 9   (0.69% var. explained)")
dev.off()

#make biplot related tables
library(factoextra)
vars = get_pca_var(newpc)

vecLenTable = matrix(data=NA, nrow=nrow(vars$coord), ncol=5)
vecLenTableSimple = matrix(data=NA, nrow=nrow(vars$coord), ncol=4)
for(rowNum in 1:nrow(vars$coord)){
	xCoord = vars$coord[rowNum,1]
	yCoord = vars$coord[rowNum,2]
	vecLen = sqrt(xCoord^2+yCoord^2)

	vecLenTable[rowNum,1] = vecLen
	vecLenTable[rowNum,2] = markers[rowNum]
	vecLenTable[rowNum,3] = rowNum
	vecLenTable[rowNum,4] = xCoord
	vecLenTable[rowNum,5] = yCoord
	
	vecLenTableSimple[rowNum,1] = vecLen
	vecLenTableSimple[rowNum,2] = rowNum
	vecLenTableSimple[rowNum,3] = xCoord
	vecLenTableSimple[rowNum,4] = yCoord
}

vecLenTable2 = vecLenTable[order(vecLenTable[,1]),] #sorted by vector length
vecLenTableSimple2 = vecLenTableSimple[order(vecLenTableSimple[,1]),] #sorted by vector length

#make subsets of top variables
n = nrow(vecLenTableSimple2)
varSamp100inds = vecLenTableSimple2[(n-99):n,][,2]
varSamp50inds = vecLenTableSimple2[(n-49):n,][,2]
varSamp25inds = vecLenTableSimple2[(n-24):n,][,2]
varSamp10inds = vecLenTableSimple2[(n-9):n,][,2]

subVar100 = SubSampVars(varSamp100inds, vars$coord)
subVar50 = SubSampVars(varSamp50inds, vars$coord)
subVar25 = SubSampVars(varSamp25inds, vars$coord)
subVar10 = SubSampVars(varSamp10inds, vars$coord)

#use vecLenTable2 to export the full sized table
write.table(vecLenTable2, file="ZeaAllInfo.pmz.noMz.pcaMarkerVecLens", sep="\t", row.names=FALSE, col.names=FALSE)

#Make a biplot with all variables
pdf("PC1vPC2woMz_biplotAll.pdf")
biplot(x=newpc$x, y=vars$coord, expand=3, xlim=c(-36,50), ylim=c(-33,23))
dev.off()

#Make a biplot with top 100 variables
pdf("PC1vPC2woMz_biplot100.pdf")
biplot(x=newpc$x, y=subVar100, expand=3, xlim=c(-36,50), ylim=c(-33,23), ylabs=varSamp100inds)
dev.off()

#Make a biplot with top 50 variables
pdf("PC1vPC2woMz_biplot50.pdf")
biplot(x=newpc$x, y=subVar50, expand=3, xlim=c(-36,50), ylim=c(-33,23), ylabs=varSamp50inds)
dev.off()

#Make a biplot with top 25 variables
pdf("PC1vPC2woMz_biplot25.pdf")
biplot(x=newpc$x, y=subVar25, expand=3, xlim=c(-36,50), ylim=c(-33,23), ylabs=varSamp25inds)
dev.off()

#Make a biplot with top 10 variables
pdf("PC1vPC2woMz_biplot10.pdf")
biplot(x=newpc$x, y=subVar10, expand=3, xlim=c(-36,50), ylim=c(-33,23), ylabs=varSamp10inds)
dev.off()


##Determine which top 100 markers are contributing primarily to the differentiation of CB and thase that are contributing primarily to the differentiation of parv, mex, and hybrids all with respect to PCs 1 and 2.
#How to: Use coordinates for top 100 markers to calculate the  angle of the vector against th horizontal line.  Use this info to make a call about primary contribution and output info in a new table along with all the info in the first table
library(NISTunits)
top100info = matrix(data=NA, nrow=100, ncol=7)
for(i in 1:nrow(subVar100)){
	xCoord = subVar100[i,1]
	yCoord = subVar100[i,2]
	vecLen = sqrt(xCoord^2+yCoord^2)
	num = varSamp100inds[i]
	marker = markers[num]
	angle = NISTradianTOdeg(atan(yCoord/xCoord))
	if(angle < 44.0){
		primCont = "parv_mex_hyb"
	}
	else if(angle > 46.0){
		primCont = "CB_others"
	}
	else{
		primCont = "Ambiguous"
	}
	
	top100info[i,1] = vecLen
	top100info[i,2] = marker
	top100info[i,3] = num
	top100info[i,4] = xCoord
	top100info[i,5] = yCoord
	top100info[i,6] = angle
	top100info[i,7] = primCont
}
colnames(top100info) = c("vectorMag","markerName","markerNum","xCoord","yCoord","angle","primaryCont")

#use top100info to export the full sized table
write.table(top100info, file="ZeaAllInfo.pmz.noMz.PC1v2_top100.pcaMarkerVecLens", sep="\t", row.names=FALSE, col.names=TRUE)



