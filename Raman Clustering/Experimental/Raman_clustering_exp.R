##########################################################
#
# Cluster analysis of Raman spectral maps
#
# Version 3-20160909b-exp
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
maxClust=6

labSpec=T	# Set to true if maps acquired with LabSpec
			# if F, also set csvAsIn = T
			# if T, also set csvAsIn = F
			
cvsAsIn=F	# Set to false is normal txt input 
csvAsOut=F  # Set to false for normal txt output

skimData=T

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

plotData <- function(i1,i2){
	plot(y[,i1],y[,i2], xlab=Par[i1], ylab=Par[i2])}
	
dataFile<-paste(rootName,"-data.pdf",sep="")
pdf(file=dataFile, width=dimPlot, height=dimPlot, onefile=T)

for (i in 2:numPar){
	plotData(1,i)
}
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

plotDataMaps <- function(i1){
	image(interp(X,Y,y[,i1],xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=Par[i1])}

dataFile<-paste(rootName,"-data-maps.pdf",sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

for (i in 1:numPar){
	plotDataMaps(i)
	}

#legend("bottomright", bg = "white", paste("Phase ", 1:numPhase), col=1:numPhase,pch = 1:numPhase%%10+15, cex = 1.5)
#grid(gridx,gridy, lwd = 1,equilogs =FALSE)

dev.off()


############################
# Clustering
############################
dataset <-matrix(NA,length(y[,1]),0)
for (i in 1:numPar) {
	dataset<-cbind(dataset,y[,i])
}
elements<-dataset


#if(limClust==F){
#	print("Cluster analysis in progress: using unlimited number of clusters...")	
#	dataclust<-Mclust(elements)} else {
#	print("Cluster analysis in progress: using fixed number of clusters...")		
#	dataclust<-Mclust(elements, G=1:maxClust)
#}
