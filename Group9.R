## This R code is just a basic overview of how to use the softImpute package. SoftImpute fits 
## a low-rank matrix approximation to a matrix with missing values via nuclear-norm
## regularization. The algorithm fills in the missing values with the current guess, 
## and then solves the optimization problem on the complete matrix using a soft-thresholded SVD. 
## An SVD object is returned, with components "u", "d", and "v". 
## For large matrices there is a special sparse-matrix class named``Incomplete'' that
## efficiently handles all computations.
## Authors: Trevor Hastie, Rahul Mazumder.

install.packages("softImpute")
install.packages("Matrix")

##Dimensions
set.seed(101)
n=200
p=100
J=50
np=n*p


## Generate a matrix (with dimensions from above) with Gaussian random entries:
missfrac=0.3
x=matrix(rnorm(n*J),n,J)%*%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
ix=seq(np)

## Take a sample of the specified size (np*missfrac) from the elements of ix withouth replacement:
imiss=sample(ix,np*missfrac,replace=FALSE)
xna=x
xna[imiss]=NA   ##the indices of the matrix given by imiss are filled with NA

### Matrix in the class "Incomplete"
xnaC=as(xna,"Incomplete")

### Use softImpute with "svd" algorithm 
fit1=softImpute(xnaC,rank=50,lambda=30,type="svd")
##Impute the missing values using complete(), which returns the full matrix
ximp=complete(xna,fit1) 

### We can standardize the matrix with biScale function
xnas=biScale(xna)  
fit2=softImpute(xnas,rank=50,lambda=10)
ximp=complete(xna,fit2)


