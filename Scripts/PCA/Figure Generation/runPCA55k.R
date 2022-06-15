#Created by David E. Hufnagel on Feb 5, 2019
#Designed to make PCA plots using the input file, "55kTeoSNPsNew_v3.txt"
#Parts were taken from Ranal.R written by Joost Van Heervaarden.
#Derived directly from 19_1_5_runPCA.R
#Updated on July 10, 2019 to add legends

#set working directory here
data = read.table("55kTeoSNPsNew_v4.txt.PCAmatrix.noMz", sep="\t", row.names=1, skip=1) #SNP file used to subset
info = read.table("55kTeoSNPsNew_v4.txt", sep="\t", row.names=1, skip=1)

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

#apply the function and set the rownames back to the individual names
data2 = apply(data,2,fcdat)
rownames(data2) = rownames(data)

#make vectors of column means and variance factors, sqrt(p*(1-p))
mv = apply(data2,2,mean,na.rm=TRUE)
p = mv/2
vv = sqrt(p*(1-p))
#take out monomorphic loci
data2 = data2[,vv>0]

#Set colors and shapes
colors = c()
shapes = c()
for (rowNum in 1:nrow(info)){
	tax = info[rowNum,3]
	hybStat = info[rowNum,9]
	hybGroup = info[rowNum,10]

	if (hybStat == "hybrid"){
		if (hybGroup == "Central_Balsas"){
			colors = c(colors, "darkorchid1")
		}
		else if(hybGroup == "East_Balsas"){
			colors = c(colors, "darkorange2")
		}
		else if(hybGroup == "South_Guerrero"){
			colors = c(colors, "orchid1")
		}
		else{
			colors = c(colors, "darkorchid4")
		}
		shapes = c(shapes, 8)
	}
	
	else if (tax == "Zea_mays_parviglumis"){
		if (hybStat == "parvHC"){
			colors = c(colors, "dodgerblue4")
			shapes = c(shapes, 1)
		}
		else if (hybStat == "parvAmb"){
			colors = c(colors, "dodgerblue2")
			shapes = c(shapes, 1)	
		}
		else{
			print("ERROR 1!")
			break
		}
	}
	
	else if (tax == "Zea_mays_mexicana"){
		if (hybStat == "mexHC"){
			colors = c(colors, "firebrick4")
			shapes = c(shapes, 2)			
		}
		else if (hybStat == "mexAmb"){
			colors = c(colors, "firebrick2")
			shapes = c(shapes, 2)	
		}
		else{
			print("ERROR 2!")
			break
		}
	}	
	else{
		print("ERROR 3!")
		break
	}
}
names(colors) = rownames(info)
names(shapes) = rownames(info)

#adding the colors and shapes to the data so they are all in 1 matrix
colM = matrix(colors, length(colors), 1)
shapesM = matrix(shapes, length(shapes), 1)
rownames(colM) = names(colors)
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
rownames(data3) = names(colors)

data3num = apply(data3, 2, as.numeric, na.rm=TRUE)[,1:ncol(data3)-1]
mv = apply(data3num,2,mean,na.rm=TRUE) #ignore last column
p = mv/2
vv = sqrt(p*(1-p))

#Handle missing data
StNor = t((t(data3num)-mv)/vv) #Converting each marker's data to a standard normal distribution
StNor[is.na(StNor)]<-0 #Set missing data to the mean of the standard normal (0)

#Run PCA, keep any values above 1% total variance explained
newpc = prcomp(StNor)

################
#summary(newpc)
#Importance of components:
#                            PC1      PC2     PC3      PC4      PC5      PC6
#Standard deviation     82.53517 49.87442 46.3761 41.86053 41.36940 37.46735
#Proportion of Variance  0.07982  0.02915  0.0252  0.02053  0.02005  0.01645
#Cumulative Proportion   0.07982  0.10896  0.1342  0.15470  0.17475  0.19120
#                            PC7      PC8      PC9     PC10     PC11
#Standard deviation     35.15125 32.21494 31.25681 30.14623 28.54863
#Proportion of Variance  0.01448  0.01216  0.01145  0.01065  0.00955
#Cumulative Proportion   0.20567  0.21783  0.22928  0.23993  0.24948

################
# Plot results
plot.new()
pdf("PC1vPC2_55k2.pdf")
par(mgp=c(2.2,1,0))
plot(newpc$x[,2]~newpc$x[,1],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]), ylab="PC 2 (2.9% variance explained)",xlab="PC 1 (8.0% variance explained)", cex.lab=1.5)
legend(list(x=-1, y=170), c("parvHC","parvAmb", "CB Hybrids","EB Hybrids", "SG Hybrids", "mexHC","mexAmb"), col=c("dodgerblue4","dodgerblue2","darkorchid1","darkorange2","orchid1","firebrick4", "firebrick2"), pch=c(15,15,8,8,8,17,17), lwd=2, lty=0, y.intersp=0.9, cex=1.5)
dev.off()

pdf("PC3vPC4_55k2.pdf")
plot(newpc$x[,4]~newpc$x[,3],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]), ylab="PC 4 (2.1% variance explained)",xlab="PC 3 (2.5% variance explained)", cex.lab=1.5)
legend(list(x=98, y=191), c("parvHC","parvAmb", "CB Hybrids","EB Hybrids", "SG Hybrids", "mexHC","mexAmb"), col=c("dodgerblue4","dodgerblue2","darkorchid1","darkorange2","orchid1","firebrick4", "firebrick2"), pch=c(15,15,8,8,8,17,17), lwd=2, lty=0, y.intersp=0.9, cex=1.5)
dev.off()

pdf("PC5vPC6_55k2.pdf")
plot(newpc$x[,6]~newpc$x[,5],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]), ylab="PC 6 (1.6% variance explained)",xlab="PC 5 (2.0% variance explained)", cex.lab=1.5)
legend(list(x=35, y=174), c("parvHC","parvAmb", "CB Hybrids","EB Hybrids", "SG Hybrids", "mexHC","mexAmb"), col=c("dodgerblue4","dodgerblue2","darkorchid1","darkorange2","orchid1","firebrick4", "firebrick2"), pch=c(15,15,8,8,8,17,17), lwd=2, lty=0, y.intersp=0.9, cex=1.5)
dev.off()

pdf("PC7vPC8_55k2.pdf")
plot(newpc$x[,8]~newpc$x[,7],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]), ylab="PC 8 (1.2% variance explained)",xlab="PC 7 (1.4% variance explained)", cex.lab=1.5)
legend(list(x=150, y=-43), c("parvHC","parvAmb", "CB Hybrids","EB Hybrids", "SG Hybrids", "mexHC","mexAmb"), col=c("dodgerblue4","dodgerblue2","darkorchid1","darkorange2","orchid1","firebrick4", "firebrick2"), pch=c(15,15,8,8,8,17,17), lwd=2, lty=0, y.intersp=0.9, cex=1.5)
dev.off()

pdf("PC9vPC10_55k2.pdf")
plot(newpc$x[,10]~newpc$x[,9],col=data3[,ncol(data3)], pch=as.numeric(data3[,ncol(data3)-1]), ylab="PC 10 (1.1% variance explained)",xlab="PC 9 (1.1% variance explained)", cex.lab=1.5)
legend(list(x=107, y=-64), c("parvHC","parvAmb", "CB Hybrids","EB Hybrids", "SG Hybrids", "mexHC","mexAmb"), col=c("dodgerblue4","dodgerblue2","darkorchid1","darkorange2","orchid1","firebrick4", "firebrick2"), pch=c(15,15,8,8,8,17,17), lwd=2, lty=0, y.intersp=0.9, cex=1.5)
dev.off()





