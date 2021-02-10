#This script is designed to generate linear regression plots of 
#altitude verses  inv4m inverison SNP state. 1 plot per SNP type.
#Created by David E. Hufnagel

#set working directory here
data = read.table("inv4mSNPbyAcc.txt", sep="\t", header=T)

alt = data$altitude
c = data$PZD00030.1_cFreq
t = data$PZD00030.1_tFreq
n = data$PZD00030.1_nFreq
color = as.character(data$hybStatus)

#This function creates a y=mx+b style formula from an lm summary object. 
#my version
GetForm = function (model) {
	sprintf("y = %.4f * %s + %.4f", coefficients(model)[-1], names(coefficients(model)[-1]), coefficients(model)[1]) 
}

#plot C
pdf("altVcSNP.pdf")
plot(alt,c, col=color, lwd=1.5, xlab="altitude (m)", ylab="frequency of C allele for PZD00030.1")
fit = lm(c~alt);abline(fit, col="black", lwd=0.6) #do a linear regression
form = GetForm(fit);text(0.8,2600, labels=form) #fit line formula
dev.off()

#plot T
pdf("altVtSNP.pdf")
plot(alt,t, col=color, lwd=1.5, xlab="altitude (m)", ylab="frequency of T allele for PZD00030.1")
fit = lm(t~alt);abline(fit, col="black", lwd=0.6) #do a linear regression
form = GetForm(fit);text(0.8,420, labels=form) #fit line formula
dev.off()

#plot N
pdf("altVnSNP.pdf")
plot(alt,n, col=color, lwd=1.5, xlab="altitude (m)", ylab="frequency of missing data 'allele' for PZD00030.1")
fit = lm(n~alt);abline(fit, col="black", lwd=0.6) #do a linear regression
form = GetForm(fit);text(0.8,2600, labels=form) #fit line formula
dev.off()


