#ultimate goal is to check if JL really preserve approximate distance

#calculate the possible minimum dimensionality given n data points and e tolerance
dimension <- function(size,eps)
{
  for(i in seq(eps)){
    return(floor((4 * log(size) / ((eps ^ 2) / 2 - (eps ^ 3) / 3))))
  }
}

#compare different sample size
dimension(6, 0.1)
dimension(1000, 0.1)
#compare different tolerance
e.v <- c(0.1, 0.5)
dimension(10000, e.v)


#create gaussian matrix with given 
#sample size(sample size after JL reduction is row size), 
#sample dimension(sample dimension is the column size) and tolerance
gaussian.matrix <- function(rows, cols, eps){
  rows          <- dimension(rows, eps)
  random.matrix <- matrix(rnorm(rows * cols, mean = 0, sd = 1), rows, cols)
  return(random.matrix)
}

#create a random gaussian matrix which can reduce dimension using JL
##mat <- gaussian.matrix(10000, 10000, 0.1)
#see if the dimension is indeed reduced
#compare with former result
##print(nrow(mat))
#check if the norm is indeed 1
##norm(1/sqrt(dimension(ncol(mat),0.1))*mat[,10],type = "2")

#apply random gaussian matrix to data, given data as random matrix here
nrow <- 1000
ncol <- 1000

created.matrix <- matrix(runif(nrow * ncol, min = 0, max = 1000),
                         nrow = nrow, ncol = ncol)

#apply the gaussian matrix to data and reduce its dimension
#get data with same size and smaller dimension in the end

projection <- function(datamat, eps){
  result <- (1 / sqrt(dimension(ncol(datamat), eps))) *
    gaussian.matrix(ncol(datamat), nrow(datamat), eps) %*% datamat
  return(result)
}

#try projection function on created.matrix
##projection(created.matrix,0.1)

#view the matrix after multiplying with gaussian matrix
##new.matrix <- projection(created.matrix,0.1)
#check if the dimension is indeed reduced and 
#check if the sample size is still the same
##print(nrow(new.matrix))
##print(ncol(new.matrix))



#now check if the approximate distance is indeed preserved
#input is the data as matrix
#if the dimension reducton does not preserve distance
#the return value will give at which data point the distance is not preserved
#if the dimension reduction preserves distance, the function will return true
#either way the function will give a histogram of the ratio of the distances before and after reduction
inspect <- function(datamat, eps){
  newmat    <- projection(datamat, eps)
  dist.old  <- dist(t(datamat)) ^ 2
  dist.new  <- dist(t(newmat)) ^ 2
  denom.hist<- as.matrix(dist.new / dist.old - 1)
  denom     <- as.matrix(abs(dist.new / dist.old - 1))
  n         <- NROW(denom)
  a         <- 0
  for (i in 1:(n - 1)){
    for (j in (i + 1):n){
      if (denom[i,j] > eps){
        a <- a + 1
        return(list(c(i, j, denom[i,j])), hist(denom.hist, xlim = c(min(denom[i,j],-eps), max(denom[i,j], eps)), freq = F))}
    }
  }
if (a == 0){
  return(list(("true"), hist(denom.hist, xlim = c(min(denom[i,j],-eps), max(denom[i,j], eps)), freq = F)))}
}

inspect(created.matrix,0.3)

