#####################################################################
#
# Autocorrelation matrix for Raman maps 
#
# Version 1.5-20160909h
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
#####################################################################

library('corrplot');library(Hmisc);library(akima); library(fields);library(plotrix);
library(spatstat);

rootName = "Cluster_matrix"
dimPlot = 10
cvsAsIn=T   # Set to false is normal txt input

if(cvsAsIn==T){
	inputFile = paste(rootName,".csv",sep="")
	t = read.csv(inputFile, header = TRUE)
	numCols = length(t)
	labels <- colnames(t[2:(numCols-2)])
} else {		
	inputFile = paste(rootName,".txt",sep="")
	m = read.table(inputFile, header = FALSE, fill = TRUE)
	t <- matrix(scan(inputFile, n = (nrow(m))*(ncol(m)), what = double(), skip = 1), nrow(m)-1, ncol(m), byrow = TRUE)
	numCols = ncol(m)
	labels = matrix(NA, ncol(m)-3, 1)
	for(i in 1:(ncol(m)-3)){
		labels[i]=as.character(m[1,i])}}

dimCols <- ncol(t)-3
Mcorr <- cor(t[,2:(numCols-2)])

pdf(file=paste(rootName,"-CorrMap.pdf",sep=""), width=dimPlot, height=dimPlot, onefile=T)

image(x=seq(dimCols), y=seq(dimCols), z=Mcorr, xlab="", ylab="", axes=FALSE)
axis(3, at = 1:dimCols, labels=labels)
axis(2, at = 1:dimCols, labels=labels)
text(expand.grid(x=seq(dimCols), y=seq(dimCols)), labels=round(c(Mcorr),2))

corrplot(Mcorr, method = "circle", type = "upper") #plot matrix
corrplot(Mcorr, method = "number", type = "upper") #plot matrix

dev.off()