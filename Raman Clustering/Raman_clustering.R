##########################################################
#
# Cluster analysis of Raman spectral maps
#
# Version 1-20150916a
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

sampleName<- "Dracken-7-tracky_rat_map2_fit2-col-clan"
# sampleName<- "Draken_intensities_map1_fit2-col-clan"

NumPar=6
 Par=c("HC","wG","D1G","D4D5G","DG","D5G")
# Par=c("HC","wG","G","D1","D5","D4")
#Par=c("A","B","C", "D","E","F")

maxClust=6

dimPlot=8
normcoord=T

skimData=F

limClust=F
plotClust=T

############################
# File name management 
############################
rootName=gsub("-clan","",sampleName)
fileName<-paste(sampleName,".txt",sep="")
file<-read.delim(fileName, header = TRUE, sep="\t" )


############################
# Load Libraries 
############################
library(mclust);library(ellipse);library(Hmisc);library(pixmap)
library(matlab); library(akima)
palette=(c("black","red","blue","magenta","green", "yellow"))


############################
# Read the data 
############################
Amatrix<-subset(file, Parameter == Par[1], select = c(Value))
Bmatrix<-subset(file, Parameter == Par[2], select = c(Value))
Cmatrix<-subset(file, Parameter == Par[3] ,select = c(Value))
Dmatrix<-subset(file, Parameter == Par[4],select = c(Value))
Ematrix<-subset(file, Parameter == Par[5],select = c(Value))
Fmatrix<-subset(file, Parameter == Par[6] ,select = c(Value))
#Gmatrix<-subset(file, Parameter == Par[7] ,select = c(Value))
Xmatrix<-subset(file, Parameter =="X" ,select = c(Value))
Ymatrix<-subset(file, Parameter =="Y" ,select = c(Value))

A<-Amatrix[,1]
B<-Bmatrix[,1]
C<-Cmatrix[,1]
D<-Dmatrix[,1]
E<-Ematrix[,1]
F<-Fmatrix[,1]
#G<-Gmatrix[,1]
Xm<-Xmatrix[,1]
Ym<-Ymatrix[,1]


#################################
# Make sure the matrix is numeric
#################################
A[is.na(A)] <- 0
B[is.na(B)] <- 0
C[is.na(C)] <- 0
D[is.na(D)] <- 0
E[is.na(E)] <- 0
F[is.na(F)] <- 0
#G[is.na(G)] <- 0


############################
# Removing spurious data
############################
if(skimData==T){
	n<-1
	while(n==1)
		{	dev.new(width=dimPlot, height=dimPlot)
		plot(A,B, xlab=Par[1], ylab=Par[2], main="Press the mouse right button to refresh, close graph to continue")
		ind<-identify(A,B,tolerance =0.1,labels="M", plot = TRUE)

		A<-A[-ind]
		B<-B[-ind]
		C<-C[-ind]
		D<-D[-ind]
		E<-E[-ind]
		F<-F[-ind]
		Xm<-Xm[-ind]
		Ym<-Ym[-ind]
		dev.off()
		}

	plot(A,B, xlab=Par[1], ylab=Par[2], main="This data will be analyzed")}
	
	
############################
# Plotting plain Raman data
############################

dataFile<-paste(rootName,"-data.pdf",sep="")
pdf(file=dataFile, width=dimPlot, height=dimPlot, onefile=T)
plot(A,B, xlab=Par[1], ylab=Par[2])
plot(A,E, xlab=Par[1], ylab=Par[6])
plot(B,E, xlab=Par[2], ylab=Par[6])
plot(B,C, xlab=Par[2], ylab=Par[3])
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

dataFile<-paste(rootName,"-plain-maps.pdf",sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

image(interp(X,Y,A), main=Par[1])
image(interp(X,Y,B), main=Par[2])

layout(matrix(c(1,2,1,2), 2, 2, byrow = T));  
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

image(interp(X,Y,C), main=Par[3])
image(interp(X,Y,D), main=Par[4])


layout(matrix(c(1,2,1,2), 2, 2, byrow = T));  
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

image(interp(X,Y,E), main=Par[5])
image(interp(X,Y,F), main=Par[6])

#legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase,pch = 1:numPhase%%10+15, cex = 1.5)
#grid(gridx,gridy, lwd = 1,equilogs =FALSE)

dev.off()


############################
# Clustering
############################

dataset<-cbind(A,B,C,D,E,F)

elements<-cbind(A,B,C,D,E,F)

if(limClust==T){
print("Full clustering...")	
dataclust<-Mclust(elements)} else {
print("Fixed phases clustering...")		
dataclust<-Mclust(elements, G=1:maxClust)
}

#####################################
# Density
#####################################

density <- dens(modelName=dataclust$modelName, data = elements, parameters = dataclust$parameters)

#image(interp(X,Y,density), col=1:5, pch = 1:20, cex.lab=1.7,main="Analysis")

#####################################
# Result extraction
#####################################

numProbe<-dataclust$n
numPhase<-dataclust$G
numInd<-dataclust$n
Z<-dataclust$z
phase<-dataclust$classification
Sort.Phase<-rep(0,numInd)
stdA<-rep(0,numPhase)
stdB<-rep(0,numPhase)
stdC<-rep(0,numPhase)
stdD<-rep(0,numPhase)
stdE<-rep(0,numPhase)
stdF<-rep(0,numPhase)

Alloc<-cbind(phase,Z)  ## to find alloacation rates
mean<-dataclust$parameters$mean

Order<-order(mean[1,]) ### sort phases based on M and return index of each one

oSIGMA<- array(0,c(NumPar,NumPar,numPhase))
oSIGMA[,,1:numPhase]<-dataclust$parameters$variance$sigma[,,Order] #### covaraince of each phase: sorted
sortedMean<-mean[,Order]  ##### Sorted mean 

Vol.F<-array(0,numPhase)
Vol.F[1:numPhase]<-dataclust$parameters$pro[Order]*100   ### sorted volume fractions

meanZ<-rep(0,numPhase)

for (j in 1:numPhase)
	{
	for (i in 1:numInd)
		{
		if (phase[i]==Order[j])
		Sort.Phase[i]<-j
		}
	meanZ[j]<-mean(subset(Alloc[,j+1],Alloc[,1]==j))
	stdA[j]<-sqrt(oSIGMA[1,1,j])
	stdB[j]<-sqrt(oSIGMA[2,2,j])
	stdC[j]<-sqrt(oSIGMA[3,3,j])
	stdD[j]<-sqrt(oSIGMA[4,4,j])
	stdE[j]<-sqrt(oSIGMA[5,5,j])
	stdF[j]<-sqrt(oSIGMA[6,6,j])
	}
omeanZ<-meanZ[Order]


Clusto<-cbind(dataset,Sort.Phase,Xm,Ym)
colnames(Clusto) <- c(Par,"Phase","X","Y")
ClusFile<-paste(rootName,"-Clustering-Details-Tot.csv",sep="")
write.csv(Clusto,file=ClusFile)


##############################################################

for(i in 1:numPhase)
	{
	phaser<-subset(Clusto,Clusto[,NumPar+1]==i)
	cov<-oSIGMA[,,i]
	write.csv(phaser,paste(rootName,"-Clustering-Details-",i,".csv",sep=""))
	}

Summary<-rbind(sortedMean,stdA,stdB,stdC,stdD,stdE,stdF,Vol.F)

##write.csv(Summary,paste(sampleName,"Clustering Details.csv",sep=" "),colNames = TRUE,sheet ="Summary",append=TRUE,from=c(4,6))
SumFile<-paste(rootName,"-summary.csv",sep="")
write.csv(Summary,file=SumFile)

#####################################
# Plots - A-B-C-D
#####################################

#dev.new(width=dimPlot*2, height=dimPlot)

dataFile<-paste(rootName,"-clust.pdf",sep="")
pdf(dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,1,2,3,1,1,4,5), 2, 4, byrow = F)); 

par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")

plot(A,B,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[2],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[1:2,1:2,i]
	mu<-sortedMean[1:2,i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}

errbar(c(1:numPhase),sortedMean[1,],sortedMean[1,]+stdA,sortedMean[1,]-stdA,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[1],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[1,])#;   title("A (A)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[2,],sortedMean[2,]+stdB,sortedMean[2,]-stdB,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[2],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[2,])#;   title("B (H)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[3,],sortedMean[3,]+stdC,sortedMean[3,]-stdC,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[3],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[3,])#;   title("A (A)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[4,],sortedMean[4,]+stdD,sortedMean[4,]-stdD,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[4],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[4,])#;   title("B (H)",cex.main=1.8)


#####################################
# Plots - E-F - Vol. Frac. - Alloc
#####################################

#dev.new(width=dimPlot*2, height=dimPlot)
layout(matrix(c(1,1,2,3,1,1,4,5), 2, 4, byrow = F)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5))

par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")

plot(C,D,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[4],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[3:4,3:4,i]
	mu<-sortedMean[3:4,i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}

errbar(c(1:numPhase),sortedMean[5,],sortedMean[5,]+stdE,sortedMean[5,]-stdE,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[5],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[5,])#;   title("A (A)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[6,],sortedMean[6,]+stdF,sortedMean[6,]-stdF,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[6],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[6,])#;   title("B (H)",cex.main=1.8)

barplot(Vol.F,cex.lab=1.5,cex.main=1.8, width = 0.2,col=1:numPhase+1,xlab="Phase",ylab="Volume Fraction [%]",names.arg = c(1:numPhase),cex.names=1.5,cex.axis=1.5, ylim=c(0,80))#;   title("Volume Fractions",cex.main=1.8)

barplot(omeanZ,cex.lab=1.5,cex.main=1.8, width = 0.2,col=1:numPhase+1,xlab="Phase",ylab="Allocation Rate",names.arg = c(1:numPhase),cex.names=1.5,cex.axis=1.5,ylim=c(0,1))#;   title("Uncertainty",cex.main=1.8)

dev.off()
 
#####################################
# Mapping the phases
#####################################

#dev.new(width=dimPlot*2, height=dimPlot)

dataFile<-paste(rootName,"-maps.pdf",sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,B,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[2],xlim=c(min(A),max(A)),ylim=c(min(B),max(B)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

if(normcoord==T){
	X<-(Xm-min(Xm))
	Y<-(Ym-min(Ym))
}else{
		X=Xm;
		Y=Ym;}	

#dev.new(width=7, height=7)
image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")
#legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase,pch = 1:numPhase%%10+15, cex = 1.5)
#grid(gridx,gridy, lwd = 1,equilogs =FALSE)

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,C,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[3],xlim=c(min(A),max(A)),ylim=c(min(C),max(C)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,D,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[4],xlim=c(min(A),max(A)),ylim=c(min(D),max(D)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,E,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[5],xlim=c(min(A),max(A)),ylim=c(min(E),max(E)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,F,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[6],xlim=c(min(A),max(A)),ylim=c(min(F),max(F)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(C,E,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[5],xlim=c(min(C),max(C)),ylim=c(min(E),max(E)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(C,D,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[4],xlim=c(min(C),max(C)),ylim=c(min(D),max(D)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(C,F,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[6],xlim=c(min(C),max(C)),ylim=c(min(F),max(F)))
legend("bottomright", paste("Phase(", 1:numPhase,")"), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")

dev.off()


#####################################
# Mapping the clustering plots
#####################################

if(plotClust==T){
	clustFile<-paste(rootName,"-clust-plots.pdf",sep="")
	pdf(file=clustFile, width=dimPlot, height=dimPlot, onefile=T)
	plot(dataclust, what = "BIC", main = "BIC")
	plot(dataclust, what = "classification", main = "classification")
	plot(dataclust, what = "uncertainty", main = "uncertainty")
	plot(dataclust, what = "density", main = "density")
	dev.off()
	}