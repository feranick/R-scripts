#####################################################################
#
# Autocorrelation matrix for Raman maps 
#
# Version 1.0-20160330
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
#####################################################################


library('corrplot');library(Hmisc);library(akima); library(fields);library(plotrix);
library(spatstat);

inputFile = "1-col.csv"

rootName=gsub(".csv","",inputFile)
dimPlot = 10;
numCols = 6

dimCols <- dim(t[,2:numCols+1])[2]

t = read.csv(inputFile, header = TRUE)

Mcorr <- cor(t[,2:numCols+1])

pdf(file=paste(rootName,"-CorrMap.pdf",sep=""), width=dimPlot, height=dimPlot, onefile=T)

image(x=seq(dimCols), y=seq(dimCols), z=Mcorr, xlab="", ylab="", axes=FALSE)
axis(3, at = 1:dimCols, labels=colnames(t[2:numCols+1]))
axis(2, at = 1:dimCols, labels=colnames(t[2:numCols+1]))
text(expand.grid(x=seq(dimCols), y=seq(dimCols)), labels=round(c(Mcorr),2))

corrplot(Mcorr, method = "circle", type = "upper") #plot matrix
corrplot(Mcorr, method = "number", type = "upper") #plot matrix

dev.off()