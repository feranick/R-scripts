#####################################################################
#
# Autocorrelation matrix for Raman maps 
#
# Version 1.5-20160909g
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
#####################################################################

library('corrplot');library(Hmisc);library(akima); library(fields);library(plotrix);
library(spatstat);

inputFile = "Cluster_matrix.csv"

rootName=gsub(".csv","",inputFile)
dimPlot = 10;

t = read.csv(inputFile, header = TRUE)
numCols = length(t)

dimCols <- dim(t[,2:numCols-2])[2]
Mcorr <- cor(t[,2:numCols-2])

pdf(file=paste(rootName,"-CorrMap.pdf",sep=""), width=dimPlot, height=dimPlot, onefile=T)

image(x=seq(dimCols), y=seq(dimCols), z=Mcorr, xlab="", ylab="", axes=FALSE)
axis(3, at = 1:dimCols, labels=colnames(t[2:numCols-2]))
axis(2, at = 1:dimCols, labels=colnames(t[2:numCols-2]))
text(expand.grid(x=seq(dimCols), y=seq(dimCols)), labels=round(c(Mcorr),2))

corrplot(Mcorr, method = "circle", type = "upper") #plot matrix
corrplot(Mcorr, method = "number", type = "upper") #plot matrix

dev.off()