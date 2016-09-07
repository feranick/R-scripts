##########################################################
#
# Cluster analysis of Raman spectral maps
#
# Version 2-n7-20160907c
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

##########################################################
# input file is direclty from the data sheet
##########################################################
inputFile<- "regions_1591_X-col.txt"

############################
# Script parameters
############################
maxClust=6

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

#Par=c("HC","wG","D5G","D1G, "D4D5G","DG") # No longer used.

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
    			{Par[i]=as.character(m[1,i+1])}}	

numPar = length(Par)-2
print("numPar: ")
print(numPar)

if(cvsAsIn==F) {
	m=read.table(inputFile, header = FALSE, fill = TRUE)
	y <- matrix(scan(inputFile, n = (nrow(m))*(ncol(m)), what = double(), skip = 1), nrow(m)-1, ncol(m), byrow = TRUE)} else {
		
	m=read.csv(inputFile, header = FALSE, skip = 1)
	#y<-asnumeric(as.character(m))}
	y<-m}


A<-y[,2]
B<-y[,3]
C<-y[,4]
D<-y[,5]
E<-y[,6]
F<-y[,7]
G<-y[,8]

if(labSpec == T) {
	Xm <- y[,10]
	Ym <- -y[,9]
	} else {
		Xm <- y[,9]
		Ym <- y[,10]}
	

############################
# Removing spurious data
############################
if(skimData==T){
	n<-1
	while(n==1)
		{	dev.new(width=dimPlot, height=dimPlot)
		plot(A,D, xlab=Par[1], ylab=Par[4], main="Press the mouse right button to refresh after a selection,\n or to continue if no point is selected")
		ind<-identify(A,D,tolerance = 0.1,labels="M", plot = TRUE)
		
		if(!length(ind)) break
		
		A<-A[-ind]
		B<-B[-ind]
		C<-C[-ind]
		D<-D[-ind]
		E<-E[-ind]
		F<-F[-ind]
		G<-G[-ind]
		Xm<-Xm[-ind]
		Ym<-Ym[-ind]
		dev.off()
		}

	plot(A,D, xlab=Par[1], ylab=Par[4], main="This data will be analyzed")}
	
	#aspratio <- length(unique(Ym))/length(unique(Xm))
	aspratio <- 1	
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

dataFile<-paste(rootName,"-data-maps.pdf",sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

image(interp(X,Y,A,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=Par[1])
image(interp(X,Y,B,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = aspratio, main=Par[2])

image(interp(X,Y,C,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = aspratio, main=Par[3])
image(interp(X,Y,D,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = aspratio, main=Par[4])

image(interp(X,Y,E,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = aspratio, main=Par[5])
image(interp(X,Y,F,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = aspratio, main=Par[6])

image(interp(X,Y,G,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = aspratio, main=Par[7])


#legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase,pch = 1:numPhase%%10+15, cex = 1.5)
#grid(gridx,gridy, lwd = 1,equilogs =FALSE)

dev.off()


############################
# Clustering
############################

dataset<-cbind(A,B,C,D,E,F,G)

elements<-cbind(A,B,C,D,E,F,G)

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
stdA<-rep(0,numPhase)
stdB<-rep(0,numPhase)
stdC<-rep(0,numPhase)
stdD<-rep(0,numPhase)
stdE<-rep(0,numPhase)
stdF<-rep(0,numPhase)
stdG<-rep(0,numPhase)

Alloc<-cbind(phase,Z)  ## to find alloacation rates
mean<-dataclust$parameters$mean

Order<-order(mean[1,]) ### sort phases based on M and return index of each one

oSIGMA<- array(0,c(numPar,numPar,numPhase))
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
	stdG[j]<-sqrt(oSIGMA[7,7,j])
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

Summary<-rbind(sortedMean,stdA,stdB,stdC,stdD,stdE,stdF,stdG,Vol.F)

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

errbar(c(1:numPhase),sortedMean[1,],sortedMean[1,]+stdA,sortedMean[1,]-stdA,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[1],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[1,])#;   title("A (A)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[2,],sortedMean[2,]+stdB,sortedMean[2,]-stdB,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[2],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[2,])#;   title("B (H)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[3,],sortedMean[3,]+stdC,sortedMean[3,]-stdC,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[3],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[3,])#;   title("A (A)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[4,],sortedMean[4,]+stdD,sortedMean[4,]-stdD,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[4],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[4,])#;   title("B (H)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[5,],sortedMean[5,]+stdE,sortedMean[5,]-stdE,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[5],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[5,])#;   title("A (A)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[6,],sortedMean[6,]+stdF,sortedMean[6,]-stdF,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[6],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[6,])#;   title("B (H)",cex.main=1.8)

errbar(c(1:numPhase),sortedMean[7,],sortedMean[7,]+stdG,sortedMean[7,]-stdG,cex.lab=1.5,pch=16,xlab="Phase",ylab=Par[7],col=c(1:numPhase+1),cex=2,cex.axis=1.5,xaxt = "n");   axis(1, at = 1:numPhase,cex.axis=1.5);lines(sortedMean[7,])#;   title("B (H)",cex.main=1.8)

barplot(Vol.F,cex.lab=1.5,cex.main=1.8, width = 0.2,col=1:numPhase+1,xlab="Phase",ylab="Volume Fraction [%]",names.arg = c(1:numPhase),cex.names=1.5,cex.axis=1.5, ylim=c(0,80))#;   title("Volume Fractions",cex.main=1.8)

barplot(omeanZ,cex.lab=1.5,cex.main=1.8, width = 0.2,col=1:numPhase+1,xlab="Phase",ylab="Allocation Rate",names.arg = c(1:numPhase),cex.names=1.5,cex.axis=1.5,ylim=c(0,1))#;   title("Uncertainty",cex.main=1.8)

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)
image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)

#####################################
# Covariance Plots
#####################################

#dev.new(width=dimPlot*2, height=dimPlot)
layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = F)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")

plot(C,D,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[4],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[c(3,4),c(3,4),i]
	mu<-sortedMean[c(3,4),i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)

plot(C,E,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[5],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[c(3,5),c(3,5),i]
	mu<-sortedMean[c(3,5),i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
	
#dev.new(width=dimPlot*2, height=dimPlot)
layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = F)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")


plot(A,C,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[3],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[c(1,3),c(1,3),i]
	mu<-sortedMean[c(1,3),i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)	


plot(A,E,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[5],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[c(1,5),c(1,5),i]
	mu<-sortedMean[c(1,5),i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)


#dev.new(width=dimPlot*2, height=dimPlot)
layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = F)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")

plot(A,D,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[4],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[c(1,4),c(1,4),i]
	mu<-sortedMean[c(1,4),i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase ), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)

plot(A,B,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[2],cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7)
#grid(NULL, lwd = 1,equilogs =FALSE)

for (i in 1:numPhase)
	{
	cov<-oSIGMA[1:2,1:2,i]
	mu<-sortedMean[1:2,i]
	lines(ellipse(cov,centre=mu,level=0.95))
	}
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)


#################################################################
### Debugging only
#################################################################

#layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = F)); 
#par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")

#plot(B,D/C,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[2],ylab="D1/G",cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7,ylim=c(0.5, 3))

#for (i in 1:numPhase) 	{
#	cov<-oSIGMA[c(2,c(4)/c(3)),c(2,c(4)/c(3)),i]
#	mu<-sortedMean[c(2,c(4)/c(3)),i]
#	lines(ellipse(cov,centre=mu,level=0.95)) }

#legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)

#plot(A,D/C,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab="D1/G",cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7,ylim=c(0.5, 3))

#for (i in 1:numPhase)	{
#	cov<-oSIGMA[c(1,c(4)/c(3)),c(1,c(4)/c(3)),i]
#	mu<-sortedMean[c(1,c(4)/c(3)),i]
#	lines(ellipse(cov,centre=mu,level=0.95)) }

#legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)

#layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = F)); 
#par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),family="serif")

#plot(E/C,C,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,ylab=Par[3],xlab="D5/G",cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7,xlim=c(0, 1))

#for (i in 1:numPhase)	{
#	cov<-oSIGMA[c(c(5)/c(3),3),c(c(5)/c(3),3),i]
#	mu<-sortedMean[c(c(5)/c(3),3),i]
#	lines(ellipse(cov,centre=mu,level=0.95))	}

#legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)

#plot(A,E/C,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab="D5/G",cex.lab=1.7,cex.main=2,cex.axis=1.7,cex=1.7,ylim=c(0, 2))

#for (i in 1:numPhase)	{
#	cov<-oSIGMA[c(1,c(5)/c(3)),c(1,c(5)/c(3)),i]
#	mu<-sortedMean[c(1,c(5)/c(3)),i]
#	lines(ellipse(cov,centre=mu,level=0.95))	}

#legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)

#################################################################
#################################################################


dev.off()
 
#####################################
# Mapping the phases
#####################################

#dev.new(width=dimPlot*2, height=dimPlot)

dataFile<-paste(rootName,"-clust-maps.pdf",sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,B,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[2],xlim=c(min(A),max(A)),ylim=c(min(B),max(B)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

if(normcoord==T){
	X<-(Xm-min(Xm))
	Y<-(Ym-min(Ym))
}else{
		X=Xm;
		Y=Ym;}	

#image(interp(X,Y,Sort.Phase), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,main="Analysis")
image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.4,main="Phase distribution",asp = aspratio)

#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,C,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[3],xlim=c(min(A),max(A)),ylim=c(min(C),max(C)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)
#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,D,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[4],xlim=c(min(A),max(A)),ylim=c(min(D),max(D)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)
#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,E,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[5],xlim=c(min(A),max(A)),ylim=c(min(E),max(E)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)
#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,F,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[6],xlim=c(min(A),max(A)),ylim=c(min(F),max(F)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)
#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(A,G,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[1],ylab=Par[7],xlim=c(min(A),max(A)),ylim=c(min(G),max(G)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)
#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(C,E,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[5],xlim=c(min(C),max(C)),ylim=c(min(E),max(E)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)
#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(C,D,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[4],xlim=c(min(C),max(C)),ylim=c(min(D),max(D)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)
#-----

layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

plot(C,F,col=Sort.Phase+1,type="p",pch=Sort.Phase+15,xlab=Par[3],ylab=Par[6],xlim=c(min(C),max(C)),ylim=c(min(F),max(F)))
legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase+1,pch = 1:numPhase%%10+15, cex = 1.5)
grid(NULL, lwd = 1,equilogs =FALSE)

image(interp(X,Y,Sort.Phase,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), col=1:numPhase+1, pch = 1:numPhase%%10+15, cex.lab=1.7,asp = aspratio)

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
	
dataFile<-paste(rootName,paste("-clust-maps-phases.pdf", sep=""),sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)

for (j in 1:numPhase){

A2 <- A
Aav = {}

for(i in 1:length(A)){
	if(Sort.Phase[i]==j){
		Aav = c(Aav, A[i])
		A2[i]=1
			}
	else {A2[i]=0}		
}


layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.0,cex.main=1.3,cex.axis=1.0,cex=1)


image(interp(X,Y,A,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=paste("Mean HC = ",round(mean(A),2),"\u00b1",round(sd(A),2)))

image(interp(X,Y,A2,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = 1, main=paste("Phase: ",j," - Average H:C = ",round(mean(Aav),2),"\u00b1",round(sd(Aav),2)), col=1:3)
}

dev.off()