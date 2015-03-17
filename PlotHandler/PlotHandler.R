#####################################################################
#
# Plot Handler
#
# Version 1-20130526a
#
# Nicola Ferralis - ferralis@mit.edu
#
# Released under GNU Public License v.3
#
#####################################################################
library(Hmisc)

inputFile="test.txt"

# Read X, Y data in matrix
t = read.table(inputFile, header = FALSE)

# Transform into (x,y) matrix
m <- matrix(scan(inputFile, n = nrow(t)*ncol(t)), nrow(t), ncol(t), byrow = TRUE)

# Transform into vectors
x=as.vector(t[,1])
y=as.vector(t[,2])

# Perform fit and summarize results
fit<-nls(y~a+b*x^2, data.frame(x,y), trace=T)
#fit<-nls(y ~ a+b*x^2, data.frame(x,y), start=list(a=0, b=1), trace=T)
summary(fit)

# Extract sigma y
sigma2=0
for(i in 1:nrow(t)){
	sigma2=((y[i]-predict(fit)[i])^2)/(nrow(t)-2)
	}
sigma=sqrt(sigma2)
sigma

# Plot data with error bars 
errbar(x,y,y+sigma,y-sigma)
lines(x, predict(fit))
title(formula(fit))