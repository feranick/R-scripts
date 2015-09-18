##########################################################
#
# Visualize clusters in Raman cluster analysis
#
# Version 2-20150918b
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

sampleName<- "Draken_ratios_map1_fit2-col-Clustering-Details-Tot-clan"

NumPar=7
 Par=c("HC","wG","D1G","D4D5G","DG","D5G", "Phase")
# Par=c("HC","wG","G","D1","D5","D4")

############################
# Visualization parameters
############################
hct = 0.2
ph = 1 # phase of interest

dimPlot=8
normcoord=F
skimData=F
limClust=F
plotClust=T
csvAsOut=F  # Set to false is normal txt output

############################
# File name management 
############################
rootName=gsub("-clan","",sampleName)
fileName<-paste(sampleName,".txt",sep="")
file<-read.delim(fileName, header = TRUE, sep="\t" )


############################
# Load Libraries 
############################
library(Hmisc);library(pixmap)
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
Gmatrix<-subset(file, Parameter == Par[7] ,select = c(Value))
Xmatrix<-subset(file, Parameter =="X" ,select = c(Value))
Ymatrix<-subset(file, Parameter =="Y" ,select = c(Value))

A<-Amatrix[,1]
B<-Bmatrix[,1]
C<-Cmatrix[,1]
D<-Dmatrix[,1]
E<-Ematrix[,1]
F<-Fmatrix[,1]
G<-Gmatrix[,1]
Xm <- Xmatrix[,1]
Ym <- Ymatrix[,1]


#################################
# Make sure the matrix is numeric
#################################
A[is.na(A)] <- 0
B[is.na(B)] <- 0
C[is.na(C)] <- 0
D[is.na(D)] <- 0
E[is.na(E)] <- 0
F[is.na(F)] <- 0
G[is.na(G)] <- 0


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
		G<-G[-ind]
		Xm<-Xm[-ind]
		Ym<-Ym[-ind]
		dev.off()
		}

	plot(A,B, xlab=Par[1], ylab=Par[2], main="This data will be analyzed")}
	
	#aspratio <- length(unique(Ym))/length(unique(Xm))
	aspratio <- 1	
############################
# Plotting plain Raman data
############################

dataFile<-paste(rootName,"-data.pdf",sep="")
pdf(file=dataFile, width=dimPlot, height=dimPlot, onefile=T)
#plot(A,B, xlab=Par[1], ylab=Par[2])
#plot(A,C, xlab=Par[1], ylab=Par[3])
#plot(B,E, xlab=Par[2], ylab=Par[6])
#plot(B,C, xlab=Par[2], ylab=Par[3])
dev.off()



############################
# Extract data from phases
############################
A2 <- A

for(i in 1:length(A)){
	if(G[i]==ph){
		if(A[i]<hct){A2[i]=1}
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


dataFile<-paste(rootName,paste("-phase-",ph,"_hct",hct,".pdf"),sep="")
pdf(file=dataFile, width=dimPlot*2, height=dimPlot, onefile=T)
layout(matrix(c(1,2,1,2), 2, 2, byrow = T)); 
par(mar=c(4,4,4,1),mai=c(0.8,0.8,0.5,0.5),cex.lab=1.3,cex.main=2,cex.axis=1.3,cex=1)

image(interp(X,Y,A,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), xlim = c(min(X), max(X)), ylim = c(min(Y), max(Y)), asp = 1, main=Par[1])
image(interp(X,Y,A2,xo=seq(min(X), max(X), length = length(unique(X))), yo=seq(min(Y), max(Y), length = length(unique(Y)))), asp = aspratio, main=paste("Phase: ",ph," - H:C threshold = ",hct), col=1:3)

dev.off()

