##########################################################
#
# Conversion of tables with variables into required MClust formatting
# 
# CSV input version
#  
# Version 20130513a
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

inputFile="MHQ.csv"

csvOut=F  # Set to false is normal txt output

# Headers and Definitions
parValue=gsub(".csv","",inputFile)
if(csvOut==TRUE){
    outputFile=paste(parValue,"-clan.csv",sep="")} else {
        
        outputFile=paste(parValue,"-clan.txt",sep="")}

# Get and Set current working directory
(WD <- getwd())
if (!is.null(WD)) setwd(WD)
print(WD)


# Read Matrix From File
m=read.table(inputFile, header = FALSE, fill = TRUE, sep=",")

c =matrix(NA, ncol(m), 1)

for(i in 1:ncol(m)){
    c[i]=as.character(m[1,i])
}

m=read.table(inputFile, header = FALSE, fill = TRUE, sep=",")

 y <- matrix(scan(inputFile, n = (nrow(m))*(ncol(m)), what = character(0), sep=","), nrow(m), ncol(m), byrow = TRUE)

d=matrix(NA , (nrow(y)-1)*(ncol(y)-1), 2)
colnames(d) <- c("Parameter", "Value")
print(colnames(d))

k=1
for(i in 2:nrow(y)){
    for(j in 2:ncol(y)){
        
        d[k,1]=c[j]
         d[k,2]=y[i,j]
         k=k+1
        
        
    }}

# Save new file
#write.table(a, file = outputFile, col.names = T, row.names = F)
if(csvOut==TRUE){
    write.csv(d, file = outputFile, row.names = F)} else {
    write.table(d, file = outputFile, quote = FALSE, sep = "\t", col.names = T, row.names = F)}