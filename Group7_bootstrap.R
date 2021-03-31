library(monomvn)
library(glmnet)
library(boot)
library(mvtnorm)

data(diabetes)

## Reformat the data
rawx <- matrix(diabetes$x,nrow=442) #extract x from diabetes
name <- attributes(diabetes$x)$dimnames[[2]]
colnames(rawx) <- name #rename x
dim(rawx) #442 rows, 10 columns
rawy <- diabetes$y #extract y from diabetes
rawdata <- cbind(rawx,rawy) 

data <- scale(rawdata)
x <- data[,1:10]
y <- data[,11]
dmatrix <- model.matrix(y~x) #design matrix

set.seed(10)

### Non-parametric bootstrap

## Define a function to do CV within bootstrap samples
cv.fn <- function(data,index){ #index is the bootstrap sample index
  x <- data[index,1:10]
  y <- data[index,11]
  cv.out <- cv.glmnet(x,y,alpha=1,nfolds=10) #alpha=1, lasso
  bestlam <- cv.out$lambda.min #the best lambda chosen by CV
  lasso.mod <- glmnet(x,y,alpha=1,lambda=bestlam)
  coef <- as.vector(coef(lasso.mod))
  return(coef)
}

#bootstrap
np.boot <- boot(data,cv.fn,R=1000)
np.t0 <- np.boot$t0 # the estimator from original data set
np.t <- np.boot$t # estimators from all the bootstrap samples
dim(np.t) # 1000, 11. Every row corresponds to a bootstrap sample
colnames(np.t) <- c("intercept","beta1","beta2","beta3","beta4",
                 "beta5","beta6","beta7","beta8","beta9","beta10")
boxplot(np.t[,2:11],horizontal=T, main="Non-parametric Bootstrap Samples",xlab="Coefficients")
abline(v=0)

#the probability for betas to be zero
np.beta.zero <- np.t==0 #np.beta.zero is a 0-1 matrix where 
                        #1 represents coefficients with value 0
np.zero.prob <- colSums(np.beta.zero)/1000
barplot(np.zero.prob[2:11],horiz=T, main="Non-parametric Bootstrap Probability of 0",xlab="Probability")



## Parametric Bootstrap

## Get parameters for estimating F 
## We will sample y from y|beat, sigma is from N(X*beta, sigma^2I_N×N)
## And we need estimators for beta and sigma
lm.mod <- lm(y~x) 
coef <- lm.mod$coefficients #least square estimator
sigma <- summary(lm.mod)$sigma #residual standard error
mean <- dmatrix %*% coef # mean is X*beta
cov <- sigma*diag(length(y)) #covariance matrix
mle <- list(mean=mean,cov=cov) 

##Define a function to sample from the parametric distribution 
rg <- function(data, mle){ #mle is the parameters we get from original data
  mean <- mle$mean
  cov <- mle$cov
  y_new <- rmvnorm(1,mean,cov) #multivariate normal distribution
  data[,11] <- y_new
  return(data)
}

##Define a function to do CV within bootstrap samples
cv.fn2 <- function(data){
  x <- data[,1:10]
  y <- data[,11]
  cv.out <- cv.glmnet(x,y,alpha=1,nfolds=10)
  bestlam <- cv.out$lambda.min
  lasso.mod <- glmnet(x,y,alpha=1,lambda=bestlam)
  coef <- as.vector(coef(lasso.mod))
  return(coef)
}

para.boot <- boot(data,cv.fn2,1000,sim="parametric",ran.gen=rg,mle=mle) 
para.t0 <- para.boot$t0 # the estimator from original data set
para.t <- para.boot$t  # estimators from all the bootstrap samples
dim(para.t) #1000,11. Every row corresponds to a bootstrap sample
colnames(para.t) <- c("intercept","beta1","beta2","beta3","beta4","beta5",
                 "beta6","beta7","beta8","beta9","beta10")
boxplot(para.t[,2:11],horizontal=T, main="Parametric Bootstrap Samples", xlab="Coefficients")
abline(v=0)

#the probability for betas to be zero
para.beta.zero <- para.t==0 #para.beta.zero is a 0-1 matrix where 
                            #1 represents coefficients with value 0
para.zero.prob <- colSums(para.beta.zero)/1000
barplot(para.zero.prob[2:11],horiz=T, main="Parametric Boostrap Probability of 0", xlab="Probability")

