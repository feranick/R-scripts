##########################################################
#
# Cluster analysis of Raman spectral maps
#
# Version 3-20160909b
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

##########################################################
# input file is direclty from the data sheet
##########################################################
inputFile<- "Cluster_matrix.txt"

############################
# Script parameters
############################
maxClust = 6

labSpec=T	# Set to true if maps acquired with LabSpec
			# if F, also set csvAsIn = T
			# if T, also set csvAsIn = F
			
cvsAsIn=F	# Set to false is normal txt input 
csvAsOut=F  # Set to false for normal txt output

skimData=F

dimPlot=8
normcoord=F
limClust=T
plotClust=T

dataSetForThreshold = 1   # dataset to be used for extraction of data from phases (1: HC)

############################
# Load Libraries 
############################
library(mclust);library(ellipse);library(Hmisc);library(pixmap)
library(matlab); library(akima)
palette=(c("black","red","blue","magenta","green", "yellow","brown"))

############################
# File load and handling
############################


if(cvsAsIn==F) {
	rootName=gsub(".txt","",inputFile)} else{
		rootName=gsub(".csv","",inputFile)}

if(csvAsOut==T){
    outputFile=paste(rootName,"-clan.csv",sep="")} else {
        
        outputFile=paste(rootName,"-clan.txt",sep="")}

# Get and Set current working directory
(WD <- getwd())
if (!is.null(WD)) setwd(WD)
print(WD)

# Read Matrix From File
if(cvsAsIn==F) {
	m=read.table(inputFile, header = FALSE, fill = TRUE)
	} else {
		m=read.csv(inputFile, header = FALSE)}

Par = matrix(NA, ncol(m)-1, 1)
for(i in 1:ncol(m)-1)
    	{if(cvsAsIn==F) 
    		{Par[i]=as.character(m[1,i])}    			
    			else 
    			{Par[i]=as.character(m[1,i+1])}
	}	

numPar = length(Par)-2
print("Number of data sets (numPar): ")
print(numPar)

if(cvsAsIn==F) {
	m=read.table(inputFile, header = FALSE, fill = TRUE)
	y <- matrix(scan(inputFile, n = (nrow(m))*(ncol(m)), what = double(), skip = 1), nrow(m)-1, ncol(m), byrow = TRUE)} else {
	m=read.csv(inputFile, header = FALSE, skip = 1)
	#y<-asnumeric(as.character(m))}
	y<-m}

############################
# Remove index column	
############################
y<-y[,-1]

	
############################
# Removing spurious data
############################
skimFunction <- function(i1,i2){
	n<-1
	while(n==1)
		{	dev.new(width=dimPlot, height=dimPlot)
		cat("Skimming data from: x =",Par[i1],"; y =",Par[i2],"\n")	
		plot(y[,i1],y[,i2], xlab=Par[i1], ylab=Par[i2], main="Press the mouse right button to refresh after a selection,\n or to continue if no point is selected")
		ind<-identify(y[,i1],y[,i2],tolerance = 0.1,labels="M", plot = TRUE)
				
		if(!length(ind)) break
		y<-y[-ind,]
		dev.off()
		}
	plot(y[,i1],y[,i2], xlab=Par[i1], ylab=Par[i2], main="This data will be analyzed")
	return(y)
}

if(skimData==T){
	y<-skimFunction(1,2)
}

	#aspratio <- length(unique(Ym))/length(unique(Xm))
	aspratio <- 1	
	
############################
# Setup coordinates	
############################
if(labSpec == T) {
	Xm <- y[,numPar+2]
	Ym <- -y[,numPar+1]
	} else {
		Xm <- y[,numPar+1]
		Ym <- y[,numPar+2]}
	
############################
# Plotting plain Raman data
############################
	
dataFile<-paste(rootName,"-data.pdf",sep="")
pdf(file=dataFile, width=dimPlot, height=dimPlot, onefile=T)

for (i in 2:numPar){
	plot(y[,1],y[,i], xlab=Par[1], ylab=Par[i])}
dev.off()

############################
# Plotting plain Raman maps
############################

if(normcoord==T){
	X<-(Xm-min(Xm))
	Y<-(Ym-min(Ym))
}else{
		X=Xm;
		Y=Ym;}	

dataFile<-paste(rootName,"-data-maps.pdf",sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

for (i in 1:numPar){
	image(interp(X,Y,y[,i],xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=Par[i])

	#legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase,pch = 1:numPhase%%10+15, cex = 1.5)
	#grid(gridx,gridy, lwd = 1,equilogs =FALSE)
	}

dev.off()


############################
# Clustering
############################
dataset <-matrix(NA,length(y[,1]),0)
for (i in 1:numPar) {
	dataset<-cbind(dataset,y[,i])
}
elements<-dataset


if(limClust==F){
	print("Cluster analysis in progress: using unlimited number of clusters...")	
	dataclust<-Mclust(elements)} else {
	print("Cluster analysis in progress: using fixed number of clusters...")		
	dataclust<-Mclust(elements, G=1:maxClust)
}

#########################################################
# Computes densities of observations in parameterized MVN mixtures. 
#########################################################

density <- dens(modelName=dataclust$modelName, data = elements, parameters = dataclust$parameters)
summary(density)
#image(interp(X,Y,density), col=1:5, pch = 1:20, cex.lab=1.7,main="Analysis")

#########################################################
#Density Estimation via Model-Based Clustering 
#########################################################

density2 <- densityMclust(data = elements, modelName=dataclust$modelName, parameters = dataclust$parameters)
summary(density2)

#####################################
# Result extraction
#####################################

numProbe<-dataclust$n
numPhase<-dataclust$G
numInd<-dataclust$n
Z<-dataclust$z
phase<-dataclust$classification
Sort.Phase<-rep(0,numInd)

Alloc<-cbind(phase,Z)  ## to find alloacation rates
mean<-dataclust$parameters$mean

Order<-order(mean[1,]) ### sort phases based on M and return index of each one

oSIGMA<- array(0,c(numPar,numPar,numPhase))
oSIGMA[,,1:numPhase]<-dataclust$parameters$variance$sigma[,,Order] #### covaraince of each phase: sorted
sortedMean<-mean[,Order]  ##### Sorted mean

Vol.F<-array(0,numPhase)
Vol.F[1:numPhase]<-dataclust$parameters$pro[Order]*100   ### sorted volume fractions

meanZ<-rep(0,numPhase)
stdM <-matrix(0,numPhase,numPar)

for (j in 1:numPhase)
	{
	for (i in 1:numInd)
		{
		if (phase[i]==Order[j])
		Sort.Phase[i]<-j
		}
	meanZ[j]<-mean(subset(Alloc[,j+1],Alloc[,1]==j))
	
	for(f in 1:numPar){
		stdM[j,f]<-sqrt(oSIGMA[f,f,j])}
	}
	
omeanZ<-meanZ[Order]

Clusto<-cbind(dataset,Sort.Phase,Xm,Ym)
colnames(Clusto) <- c(Par[1:numPar],"Phase","X","Y")
if(csvAsOut==T){
	print("Save clustering file as csv")	
	write.csv(Clusto,paste(rootName,"-clust-all.csv",sep=""))
} else {
	print("Save clustering file as txt")	
	write.table(Clusto,paste(rootName,"-clust-all.txt",sep=""), quote = FALSE, sep = "\t", col.names = NA)
}


##############################################################

for(i in 1:numPhase)
	{
	phaser<-subset(Clusto,Clusto[,numPar+1]==i)
	cov<-oSIGMA[,,i]
	if(csvAsOut==T){	
		write.csv(phaser,paste(rootName,"-clust-",i,".csv",sep=""))
	} else {
		write.table(phaser, paste(rootName,"-clust-",i,".txt",sep=""), quote = FALSE, sep = "\t", col.names = NA)}	
	}

Summary<-sortedMean
for(i in 1:numPar){
	Summary<-rbind(Summary,stdM[,i])}
Summary<-rbind(Summary,Vol.F)

rownames(Summary) <- append(Par[1:numPar],append(paste("sd_",Par[1:numPar],sep=""),"Vol.Fr"))
colnames(Summary) <- colnames(sortedMean,do.NULL=FALSE, prefix="Phase_")

if(csvAsOut==T){
    	SumFile=paste(rootName,"-clust-summary.csv",sep="")
    	write.csv(Summary,file=SumFile)
} else {
    	SumFile=paste(rootName,"-clust-summary.txt",sep="")
    	write.table(Summary, file = SumFile, quote = FALSE, sep = "\t", col.names = NA)}
    	
    	
############################################
# Phase, allocation, Colume fraction plots
############################################

#dev.new(width=dimPlot*2, height=dimPlot)

dataFile<-paste(rootName,"-clust-plots.pdf",sep="")
pdf(dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = F)); 

par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")

for (i in 1:numPar){
	errbar(c(1:numPhase),sortedMean[i,],sortedMean[i,]+stdM[,i],sortedMean[i,]-stdM[,i],cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[i],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[i,])#;   title("A (A)",cex.main=1.8)
}
layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)
barplot(Vol.F,cex.lab=1.5,cex.main=1.8, width = 0.2,col=1:numPhase+1,xlab="Phase",ylab="Volume Fraction [%]",names.arg = c(1:numPhase),cex.names=1.5,cex.axis=1.5, ylim=c(0,80))#;   title("Volume Fractions",cex.main=1.8)

barplot(omeanZ,cex.lab=1.5,cex.main=1.8, width = 0.2,col=1:numPhase+1,xlab="Phase",ylab="Allocation Rate",names.arg = c(1:numPhase),cex.names=1.5,cex.axis=1.5,ylim=c(0,1))#;   title("Uncertainty",cex.main=1.8)

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)
image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)

dev.off()


 
#####################################
# Mapping the phases
#####################################


plotClustMaps <-function(i1,i2){
	plot(y[,i1],y[,i2],col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[i1],ylab=Par[i2],xlim=c(min(y[,i1]),max(y[,i1])),ylim=c(min(y[,i2]),max(y[,i2])))
	legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
	grid(NULL, lwd = 1,equilogs =FALSE)

	#image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")
	image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.4,main="Phase distribution",asp = aspratio)
}

#dev.new(width=dimPlot*2, height=dimPlot)
dataFile<-paste(rootName,"-clust-maps.pdf",sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)
	
for (j in 1:numPar){
	for (i in 1:numPar){
		if(i!=j)
			{plotClustMaps(j,i)}}}

dev.off()


#####################################
# Mapping the clustering plots
#####################################

if(plotClust==T){
	clustFile<-paste(rootName,"-clust-analytics.pdf",sep="")
	pdf(file=clustFile, width=dimPlot, height=dimPlot, onefile=T)
	plot(dataclust, what = "BIC", main = "BIC")
	plot(dataclust, what = "classification", main = "classification")
	plot(dataclust, what = "uncertainty", main = "uncertainty")
	plot(dataclust, what = "density", main = "density")
	dev.off()
	}
	
############################
# Extract data from phases
############################
	
extractDataPhases <-function(i1){
	dataFile<-paste(rootName,paste(paste("-clust-maps-phases_",Par[i1], sep=""),".pdf", sep=""),sep="")
	pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
	R<-y[,i1]
	for (j in 1:numPhase){
		R2 <- R
		Rav = {}
		for(i in 1:length(R)){
			if(Sort.Phase[i]==j){
				Rav = c(Rav, R[i])
				R2[i]=1}
		else {R2[i]=0}}
	
	layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
	par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.0,cex.main=1.3,cex.axis=1.0,cex=1)

	image(interp(X,Y,R,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=paste("Mean HC = ",round(mean(R),2),"\u00b1",round(sd(R),2)))

	image(interp(X,Y,R2,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = 1, main=paste("Phase: ",j," - Average H:C = ",round(mean(Rav),2),"\u00b1",round(sd(Rav),2)), col=1:3)
	}
	dev.off()}
	
extractDataPhases(dataSetForThreshold)
