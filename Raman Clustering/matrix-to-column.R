##########################################################
#
# Conversion of LabSpec Maps/Matrices into single columns
#
# Version 3-20150929a
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

XName="X"
YName="Y"
csvOut=T  # Set to false is normal txt output

##########################################
# Get list of Files
##########################################

listOfFiles <- list.files(pattern= ".txt")
inputFile<-as.matrix(listOfFiles)

# Remove extension from file name
for (i in 1:nrow(inputFile)){
	
	parValue=as.matrix(gsub(".txt","",inputFile))
	
if(csvOut==TRUE){
    outputFile=paste(parValue,"-col.csv",sep="")} else {    
    outputFile=paste(parValue,"-col.txt",sep="")}
	
}

# Get and Set current working directory
(WD <- getwd())
if (!is.null(WD)) setwd(WD)
print(WD)

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
        k=k+1
    }}

	# Save new file
	#write.table(a, file = outputFile, col.names = T, row.names = F)
	if(csvOut==TRUE){
    		write.csv(a, file = outputFile[p], row.names = F)} else {
        write.table(a, file = outputFile[p], quote = FALSE, sep = "\t", col.names = T, row.names = F)}
}