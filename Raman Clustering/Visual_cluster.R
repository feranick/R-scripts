##########################################################
#
# Visualize clusters in Raman cluster analysis
#
# Version 2-20150918d
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

##########################################################
# input file is directly used from the clustering analysis
##########################################################
inputFile="Draken_ratios_map1_fit2-col-Clustering-Details-Tot.txt"

############################
# Visualization parameters
############################
hct = 0.2 # H:C threshold
ph = 2 # phase of interest

############################
# Script parameters
############################
dimPlot=8
normcoord=F
skimData=F
limClust=F
plotClust=T
csvAsOut=F  # Set to false is normal txt output

############################
# Load Libraries 
############################
library(Hmisc);library(pixmap)
library(matlab); library(akima)
palette=(c("black","red","blue","magenta","green", "yellow"))

############################
# File load and handling
############################

parValue=gsub(".txt","",inputFile)
if(csvAsOut==TRUE){
    outputFile=paste(parValue,"-clan.csv",sep="")} else {
        
        outputFile=paste(parValue,"-clan.txt",sep="")}

# Get and Set current working directory
(WD <- getwd())
if (!is.null(WD)) setwd(WD)
print(WD)


# Read Matrix From File
m=read.table(inputFile, header = FALSE, fill = TRUE)

Par =matrix(NA, ncol(m)-1, 1)

for(i in 1:ncol(m)-1){
    Par[i]=as.character(m[1,i])
}

m=read.table(inputFile, header = FALSE, fill = TRUE)

y <- matrix(scan(inputFile, n = (nrow(m))*(ncol(m)), what = double(), skip = 1), nrow(m)-1, ncol(m), byrow = TRUE)

A<-y[,2]
B<-y[,3]
C<-y[,4]
D<-y[,5]
E<-y[,6]
F<-y[,7]
G<-y[,8]
Xm<-y[,9]
Ym<-y[,10]


############################
# Removing spurious data
############################
if(skimData==T){
	n<-1
	while(n==1)
		{dev.new(width=dimPlot, height=dimPlot)
		plot(A,B, xlab=Par[1], ylab=Par[2], main="Press the mouse right button to refresh, close graph to continue")
		ind<-identify(A,B,tolerance =0.1,labels="M", plot = TRUE)

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

	plot(A,B, xlab=Par[1], ylab=Par[2], main="This data will be analyzed")}
	

############################
# Extract data from phases
############################
A2 <- noquote(A)

for(i in 1:length(A)){
	if(G[i]==ph){
		if(A[i] <= hct) {A2[i]=1}
			else{A2[i]=2}
			}
	else {A2[i]=0}		
}

############################
# Plotting plain Raman maps
############################


if(normcoord==T){
	X<-(Xm-min(Xm))
	Y<-(Ym-min(Ym))
}else{
		X=Xm;
		Y=Ym;}	


dataFile<-paste(parValue,paste("-phase-",ph,"_hct",hct,".pdf"),sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)


image(interp(X,Y,A,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=Par[1])
image(interp(X,Y,A2,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = 1, main=paste("Phase: ",ph," - H:C threshold = ",hct), col=1:3)

dev.off()

