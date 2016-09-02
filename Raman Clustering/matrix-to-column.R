##########################################################
#
# Conversion of LabSpec Maps/Matrices into single columns
# for cluster analysis
#
# Version 4-20160902a
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

fullOutput = "Cluster_matrix.txt"
fullOut=T # Set to true for full matrix ready for clustering

XName="X"
YName="Y"
csvOut=T  # Set to false is normal txt output

##########################################
# Get list of Files
##########################################

listOfFiles <- list.files(pattern= ".txt")
inputFile<-as.matrix(listOfFiles)

##########################################
# Remove extension from file name
##########################################

for (i in 1:nrow(inputFile)){
	parValue=as.matrix(gsub(".txt","",inputFile))
	if(csvOut==TRUE){
    		outputFile=paste(parValue,"-col.csv",sep="")} else {    
    		outputFile=paste(parValue,"-col.txt",sep="")}	
}

##########################################
# Get and Set current working directory
##########################################

(WD <- getwd())
if (!is.null(WD)) setwd(WD)
	print(WD)

##############################################
# Read and initialize Full Matrix (if needed)
##############################################

if(fullOut==TRUE){
	m = read.table(inputFile[1], header = FALSE, skip =1)
	x <- matrix(scan(inputFile[1], n = nrow(m)*ncol(m), skip=1), nrow(m), ncol(m), byrow = TRUE)
	fullMatrix <- matrix(NA , (nrow(x))*(ncol(x)-1), nrow(inputFile)+3)

	#Initialize column names for fullMatrix
	colnames(fullMatrix) <- vector(mode="character", length=nrow(inputFile)+3)
	colnames(fullMatrix)[1] <- ""
	for(q in 1:nrow(inputFile)){
		colnames(fullMatrix)[q+1] <- parValue[q]	}
	colnames(fullMatrix)[nrow(inputFile)+2] <- XName
	colnames(fullMatrix)[nrow(inputFile)+3] <- YName
	print(colnames(fullMatrix))
	}

##############################################
# Process and save each matrix
##############################################

for(p in 1:nrow(inputFile)){
	# Read Y coordinates from File
	m1=read.table(inputFile[p], header = FALSE, nrow =1)
	x1 <- matrix(scan(inputFile[p], n = nrow(m1)*ncol(m1)), nrow(m1), ncol(m1), byrow = TRUE)

	# Reorder coordinates to match the matrix
	for(i in seq(ncol(x1),1,-1)){
        x1[i+1]=x1[i]
        }
	x1[1]=0

	# Read data in matrix
	m = read.table(inputFile[p], header = FALSE, skip =1)
	x <- matrix(scan(inputFile[p], n = nrow(m)*ncol(m), skip=1), nrow(m), ncol(m), byrow = TRUE)

	# Define new matrix
	a=matrix(NA , (nrow(x))*(ncol(x)-1), 4)

	# Define header for new matrix
	colnames(a) <- c("",parValue[p],XName,YName)
	print(colnames(a))

	# Transfer data into new matrix
	k=1
	for(i in 1:nrow(x)){
    		for(j in 2:ncol(x)){
        a[k,1]=k
        a[k,2]=x[i,j]
        a[k,3]=x[i,1]
        a[k,4]=x1[j]
        
        # Process data full matrix
        if(fullOut==TRUE){
        		fullMatrix[k,p+1]=x[i,j]
        		if(p==1){
        			fullMatrix[k,1]=k
        			fullMatrix[k,nrow(inputFile)+2]=x[i,1]
        			fullMatrix[k,nrow(inputFile)+3]=x1[j]
        		}}
        k=k+1
    }}

	# Save new file
	#write.table(a, file = outputFile, col.names = T, row.names = F)
	if(csvOut==TRUE){
    		write.csv(a, file = outputFile[p], row.names = F)} else {
        write.table(a, file = outputFile[p], quote = FALSE, sep = "\t", col.names = T, row.names = F)}
}

##############################################
# Save full matrix (if needed)
##############################################

if(fullOut==TRUE){
	print("Full matrix for clustering saved in: Cluster_matrix.txt")
	write.table(fullMatrix, file = fullOutput, quote = FALSE, sep = "\t", col.names = T, row.names = F)}