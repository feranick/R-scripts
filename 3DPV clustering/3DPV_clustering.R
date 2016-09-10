##########################################################
#
# Group data from several files into one with filename
#
# Version 1-20160909f
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
##########################################################

library(plyr)
listOfFiles <- list.files(pattern= ".txt")
ft<-as.matrix(listOfFiles)

# Remove extension from file name
for (i in 1:nrow(ft)){
	f=gsub(".txt","",ft)}

# Method 1
#d <- ldply(listOfFiles, read.table)
#str(d)

# Method 2
#e <- do.call(rbind, lapply(listOfFiles, read.table))  # set in rows
e <- do.call(cbind, lapply(listOfFiles, read.table))   # set in columns

am<-as.matrix(e)
b=matrix("" , nrow(am)+1,nrow(f) )
    for (i in 1:nrow(f)){
		b[1,i]=f[i,]
		for (j in 1:nrow(am)){
			b[j+1,i]=am[j,2*i]}}

tb<-as.matrix(t(b))

#write.table(b,file="1b.csv", col.names = FALSE, row.names = F, sep=",")
#write.table(tb,file="1tb.csv", col.names = FALSE, row.names = F, sep=",")

# Method 1: through direct definition of matrix
#elements=matrix("", nrow(tb), ncol(tb)-1)
#for (k in 2:ncol(tb)) {
#	elements[,k-1]=tb[,k]
#	}

# Method 2: thorough cbind
elements<-NULL
for (k in 2:ncol(tb)){
	elements<-cbind(elements, tb[,k])}
colnames(elements)<-seq(1,ncol(tb)-1)
class(elements)<-"numeric"
elements
print(sprintf("Number of triangles: %d", ncol(tb)-1))
#write.table(elements,file="clus.csv", col.names = FALSE, row.names = F, sep=",")

library(mclust)
dclust<-Mclust(elements)
plot(dclust)
