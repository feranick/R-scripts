##########################################################
#
# Group data from several files into one with filename
#
# Version 1-20130421a
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
	
	f=gsub(".txt","",ft)
	
}

# Method 1
#d <- ldply(listOfFiles, read.table)
#str(d)

# Method 2
#e <- do.call(rbind, lapply(listOfFiles, read.table))  # set in rows
e <- do.call(cbind, lapply(listOfFiles, read.table))   # set in columns

am<-as.matrix(e)

b=matrix("" , nrow(am)+1,nrow(f)+1 )

for (i in seq(1,nrow(f),1)){
	b[1,i+1]=f[i,]
	for (j in 1:nrow(am)){
		b[j+1,1]=am[j,1]
		b[j+1,i+1]=am[j,2*i]
	}
	
	
}

write.table(b,file="1b.csv", col.names = FALSE, row.names = F, sep=",")





