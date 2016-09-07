##########################################################
#
# Visualize clusters in Raman cluster analysis
#
# Version 3-20160907c
# (three regions H:C)
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

##########################################################
# input file is directly used from the clustering analysis
##########################################################
inputFile="Cluster_matrix-clust-all.txt"

############################
# Visualization parameters
############################
hct = 0.6 # high H:C threshold
hct2 = 0.5 # low H:C threshold

############################
# Script parameters
############################
dimPlot=8
normcoord=F
skimData=F
limClust=F
plotClust=T
csvAsOut=F  # Set to false is normal txt output
legend = T

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
if(csvAsOut==T){
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

numPhases <- length(unique(G))

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
	

#######################################
# If needed, normalize coordinates
#######################################	
	
if(normcoord==T){
	X<-(Xm-min(Xm))
	Y<-(Ym-min(Ym))
}else{
		X=Xm;
		Y=Ym;}	
		
		
############################
# Extract data from phases
############################

dataFile<-paste(parValue,paste("-plot-phases_hct1-",hct2,"_hct2-",hct,".pdf", sep=""),sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)

for (j in 1:numPhases){

A2 <- A
Aav = {}

for(i in 1:length(A)){
	if(G[i]==j){
		Aav = c(Aav, A[i])
		if(A[i] >= hct) {A2[i]=3}
			else{
				if(A[i] >= hct2) {A2[i]=2}
					else {A2[i]=1}
				}}
	else {A2[i]=0}		
	}


layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.0,cex.main=1.3,cex.axis=1.0,cex=1)


image(interp(X,Y,A,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=paste("Mean H:C = ",round(mean(A),2),"\u00b1",round(sd(A),2)))

image(interp(X,Y,A2,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = 1, main=paste("Phase: ",j,"- Average H:C = ",round(mean(Aav),2),"\u00b1",round(sd(Aav),2)), col=1:4)

if(legend == TRUE){
	legend("bottomright", c("Other phases", paste("H:C < ", hct2, sep=""), paste(paste(hct2, " < H:C < ", sep=""), hct, sep=""), paste("H:C > ", hct, sep="")), col=1:4,pch = 1:3%%10+15, cex = 1.2, bg = "white")
	}

}

dev.off()