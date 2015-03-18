#####################################################################
#
# VisualMap 
#
# Conversion of LabSpec Maps/Matrices into single columns with plots
#
# Version 3-20150318
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
#####################################################################

library(plyr); library(akima); library(Hmisc);library(akima); 
library(fields);library(plotrix);library(spatstat);

csvOut=T  # Set to false is normal txt output
outFile="Collective_map1"
dimPlot= 7  #Dimension of the plots in output (7 default)

XName="X"
YName="Y"

pal = palette()[1:8]

##########################################
# Get list of Files
##########################################

listOfFiles <- list.files(pattern= ".txt")
ft<-as.matrix(listOfFiles)

# Remove extension from file name
for (i in 1:nrow(ft)){
	
	f=as.matrix(gsub(".txt","",ft))
	
}


##########################################
# Read and rescale index for Y coordinates
##########################################

m1=read.table(ft[1], header = FALSE, nrow =1)
x1 <- matrix(scan(ft[1], n = nrow(m1)*ncol(m1)), nrow(m1), ncol(m1), byrow = TRUE)

# Reorder coordinates to match the matrix
for(i in seq(ncol(x1),1,-1)){
        x1[i+1]=x1[i]
        }
x1[1]=0



##########################################
# Read Coordinates
##########################################

# Read X, Y data in matrix
m = read.table(ft[1], header = FALSE, skip =1)
x <- matrix(scan(ft[1], n = nrow(m)*ncol(m), skip=1), nrow(m), ncol(m), byrow = TRUE)

X = matrix(NA , (nrow(x)), 1)
Y = matrix(NA , (ncol(x)-1), 1)

k=1
for(i in 1:nrow(x)){
    for(j in 1:ncol(x)-1){
        
        X[i]=x[i,1]
        Y[j]=x1[j+1]
        k=k+1
    }}


# Format coordinates for file
Cm=matrix(NA , (nrow(x))*(ncol(x)-1), 2)

# Transfer data into new matrix
k=1
for(i in 1:nrow(x)){
    for(j in 2:ncol(x)){
        
        Cm[k,2]=-x[i,1]
        Cm[k,1]=x1[j]
        k=k+1
    }}


##########################################
# Read Actual Data from multiple maps
##########################################

data<-matrix(NA, (nrow(x))*(ncol(x)-1), nrow(ft))


for(p in 1:nrow(ft)){
	
	m = read.table(ft[p], header = FALSE, skip =1)
	x <- matrix(scan(ft[p], n = nrow(m)*ncol(m), skip=1), nrow(m), ncol(m), byrow = TRUE)
	numProg <- fdata<-matrix(NA, (nrow(x))*(ncol(x)-1), 1)

k=1
for(i in 1:nrow(x)){
    for(j in 2:ncol(x)){
        
        data[k,p]=x[i,j]
        numProg[k]=k
        k=k+1
    }}	
    
    }

fdata<-matrix(NA, (nrow(x))*(ncol(x)-1), nrow(ft)+3)

# combine data into one matrix    
fdata<-cbind(numProg,data,Cm)  

colnames(fdata) <- c("",f,XName,YName)

##########################################
# Save to File
##########################################


#write.table(a, file = outputFile, col.names = T, row.names = F)
if(csvOut==TRUE){
    write.csv(fdata, file = paste(outFile,"-data.csv",sep=""), row.names = F)} else {
        write.table(fdata, file = paste(outFile,"-data.txt",sep=""), quote = FALSE, sep = "\t", col.names = T, row.names = F)}
        
##########################################
# Save plots of maps
##########################################        
        
        

pdf(file=paste(outFile,"-maps.pdf",sep=""), width=dimPlot, height=dimPlot, onefile=T)


for(p in 1:nrow(ft)){
	
	interpol = interp(Cm[,1],Cm[,2],data[,p])
	image.plot(interpol, cex.lab=1.7,main=f[p], legend.args=list( text="",cex=1.0, side=3, line=1), zlim=c(min(data[,p]),max(data[,p])),xlab="[um]",ylab="[um]")
	
		image.plot(interpol, cex.lab=1.7,main=f[p], legend.args=list( text="",cex=1.0, side=3, line=1), zlim=c(min(data[,p]),max(data[,p])),xlab="[um]",ylab="[um]", col = pal)
	
	interpol2 = as.im(interpol)
	plot(interpol2, main=f[p])
}
dev.off()
