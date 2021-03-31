## Simulate a Dataset
set.seed(0)

simulate_X <- function(n, p) {
  ## Set up X
  X <- matrix(rnorm(n * p), n, p)
  
  ## Scale X: let sum(X[,j]^2)=n
  X <- apply(X, 2, function(X) {
  		return((X - mean(X)) / sqrt(var(X) * (n - 1) / n))
    })
  return(X)
}

X <- simulate_X(n = 10000, p = 10)

## Check sum(X[, j]^2)

apply(X ^ 2, 2, sum)

## Simulate Y

simulate_Y <- function(X, n, p) {
  b <- (1:p) / 1.2
  Y <- X %*% b + rnorm(n)
  ## Center Y
  Y <- Y - mean(Y)
  return(Y)
}

Y <- simulate_Y(X, n = 10000, p = 10)

## Apply Coordinate Descent to calculate coefficents of parameters

Coordinate_Descent <- function(X, Y, iteration, lambda, p, n) {
  b = rep(0, p) 										# Initialize beta
  
  for (i in 1:iteration) {
    for (j in 1:p) {
      r_partial <- Y - X[, -j] %*% as.matrix(b[-j])		# Calculate partial residual
      xr 		<- sum(r_partial * X[, j]) / n  		# Calculate the inner product between r_partial and X[,j]
      x_squared <- sum((X[, j]) ^ 2) / n  				# Calculate the inner product between X^T and X

      b[j]		<- (abs(xr) - lambda) / x_squared
      b[j]		<- sign(xr) * ifelse(b[j] > 0, b[j], 0) # Soft thresholding      
    }    
  }  
  return(b)
}


## Calculate coefficients of parameters using glmnet package at lambda = 0.05

library(glmnet)

glm_lasso <-
  glmnet(
    X,
    Y,
    alpha = 1,
    lambda = 0.05,
    intercept = F,
    standardize = F,
    thresh = 1.0e-15
  )

beta_CD = Coordinate_Descent(X, Y, iteration = 1000, 
                             lambda = 0.05, n = 10000, p = 10)

max(abs(coef(glm_lasso)[-1] - beta_CD))


## Check the difference between the coefficients obtained by glmnet package and
##	 the coefficients obtained by manually implementing coordinate descent

## Case 1: Vary the number of iterations; Check the maximinum absolute difference

iteration_value <- c(1:10, 20, 50, 100, 250)

beta_diff_max = c()
beta_diff_sum = c()

for (i in iteration_value) {
  beta_CD <- Coordinate_Descent(X,
                                Y,
                                iteration = i,
                                lambda = 0.05,
                                p = 10,
                                n = 10000)
  beta_diff 	<- beta_CD - coef(glm_lasso)[-1]
  beta_diff_max	<- c(beta_diff_max, max(abs(beta_diff)))
  beta_diff_sum	<- c(beta_diff_sum, sum(abs(beta_diff)))
}

beta_diff_max
beta_diff_sum

## Plot the result 
par(mfrow = c(2, 1))
plot(iteration_value,
     beta_diff_max,
     col  = "red",
     main = "Maximum absolute difference vs. Number of iteration")
plot(iteration_value,
     beta_diff_sum,
     col  = "green",
     main = "Sum of absolute difference vs. Number of iteration")

## Case 2: For different values of n, check how many iterations needed for the estimation of coefficients to converge

## convergence threshold: 1e-15

## Fix lambda = 0.05, p = 10, iteration = 1000

number_n = c(10, 50, 100, 500, 1000, 10000)

## Modify the function to obtain how many iterations are needed for convergence

Coordinate_Descent1 <- function(X, Y, iteration, lambda, p, n) {
  beta_curr	<- rep(0, p) # Initialize beta for the current iteration
  beta_prev	<- rep(0, p) # Initialize beta before the iteration
  change 	  <- matrix(nrow = iteration, ncol = p)
  for (i in 1:iteration) {
    for (j in 1:p) {
      r_partial 	<- Y - X[,-j] %*% as.matrix(beta_curr[-j])  
      xr			<- sum(r_partial * X[, j]) / n  
      x_squared 	<- sum((X[, j]) ^ 2) / n  
      beta_curr[j]	<- (abs(xr) - lambda) / x_squared
      beta_curr[j]	<- sign(xr) * ifelse(beta_curr[j] > 0, beta_curr[j], 0)     
    }
    
    change[i, ]	<- beta_curr - beta_prev	# Calculate the diffence in coefficients after one iteration

    if (max(abs(change[i, ])) < 1e-15) {
      break  								# If the difference is smaller than the threshold, break the loop 
    }
    beta_prev <- beta_curr					# Update beta before next iteration
  }
  
  return(i)
}

iteration_conv <- c()

for (k in number_n){
  X_new		  <- simulate_X(n = k, p = 10)
  Y_new	  	<- simulate_Y(X = X_new, n = k, p = 10)
  num_iter	<- Coordinate_Descent1(X = X_new, Y = Y_new, 
                                 	iteration = 2000, lambda = 0.05, 
                                 	p = 10, n = k)
  iteration_conv <- c(iteration_conv, num_iter)
}

## Display the result
matrix(cbind(number_n, iteration_conv), 
       nrow = length(number_n), ncol = 2, 
       dimnames = list(c(1:length(number_n)), c("Number of n","Iteration")))


## Case 3: Vary the number of predictors p; Check how many iterations are needed for convergence

#  number_p <- c(5, 10, 20, 30, 50, 60)
number_p <- c(5, 10, 20, 30, 40)  
  ## Fix n = 10000, max_iteration = 1000, lambda = 0.05
  
  iteration_conv <- c()

  for (p in number_p){
    X_new	 	  <- simulate_X(n = 10000, p)
    Y_new		  	<- simulate_Y(X = X_new, n = 10000,p)
    num_iter		<- Coordinate_Descent1(X = X_new, Y = Y_new, 
                                    	iteration = 1000, lambda = 0.05, 
                                    	p, n = 10000)
    iteration_conv	<- c(iteration_conv, num_iter)
  }
  
  matrix(cbind(number_p, iteration_conv), nrow = length(number_p), 
         ncol = 2, dimnames = list(c(1:length(number_p)), c("Number of p","Iteration")))




#################################################################################
## Compare Lasso and  LAR solution on diabetes dataset

# install.packages("lars")
library(lars)
## load diabetes data from 
data(diabetes)
colnames(diabetes$x)		#442 obs with 10 columns
colnames(diabetes$x2)		#422 obs with 64 columns: include interaction terms

## Rescale X: sum(x[,j])^2) equals to 1

diabetes_X <- apply(diabetes[, 'x'], 2, 
                    function(X) return((X - mean(X)) / (sqrt(var(X) * (length(X) - 1)))))
apply(diabetes_X ^ 2, 2, sum)

## center Y 
diabetes_Y <- diabetes[, "y"] - mean(diabetes[, "y"])

## Fit Least Angle Regression 

lar.fit <- lars(diabetes_X, diabetes_Y, type = "lar", 
          		 intercept = TRUE, normalize = FALSE)

lasso.fit <- glmnet(diabetes_X, diabetes_Y, intercept = TRUE, alpha = 1,
              	standardize = FALSE, thresh = 1e-16)

## Plot the solution path from LAR and Lasso

par(mfrow = (c(1, 2)))
plot(lar.fit)
plot(lasso.fit, main = "LASSO")


## Get lambda values from LAR

lambda_lar <- lar.fit$lambda
#lasso.fit$lambda

## Compare solution for corresponding lambda value 
lasso.fit1 <-
  glmnet(diabetes_X, diabetes_Y, intercept = TRUE, alpha = 1,
         standardize = FALSE, thres = 1e-18, 
         lambda = lambda_lar/nrow(diabetes_X))
  
## Calculate the difference in estimation of parameters

max_diff <- c()

for (i in 1:length(lambda_lar)){
  		diff 		  <- coef(lar.fit)[i, ] - coef(lasso.fit1)[-1, i] 
  		max_diff	<- c(max_diff, max(abs(diff)))
}

max_diff

## Compare LAR with OLS

OLS.fit <- lm(diabetes_Y~diabetes_X)

abs(coef(OLS.fit)[-1] - coef(lar.fit)[11, ])

# --------------------------------------------------------------------

# Screening rules

# --------------------------------------------------------------------

## From Trevor Hastie, Robert Tibshirani, and 
## Martin Wainwright (2015), Statistical learning with sparsity: 
## the Lasso and generalizations, p127: 
## "The first variable to enter the model has largest 
## absolute inner-product lambda_max=max_j |x_j^T*y|,
## which also defines the entry value for lambda."

## attaching necessary libraries and dataset, preparing said dataset
library(glmnet)
library(ISLR)
data(Hitters)
Hitters <- na.omit(Hitters)

## define x and y, as we cannot use y~x notation in glmnet
x <- model.matrix(Salary~., Hitters)[, -1]
y <- Hitters$Salary

## standardisation
x <- apply(x, 2, function(k) (k - mean(k)) / 
             sqrt(sum((k - mean(k)) ^ 2) / length(k)))
## not necessary to get correct results
y <- (y - mean(y)) / sd(y) 

## grid for lambda values - refining step-by-step to find lambda_max
grid <- 10 ^ seq(from = 4, to = -3, length = 100)  
grid <- 10 ^ seq(from = 0, to = -1, length = 100)  
grid <- seq(from = 0.57, to = 0.55, length = 100)  
grid <- seq(from = 0.5657, to = 0.5659, length = 100)  
grid <- seq(from = 0.5658859, to = 0.5658879, length = 100)  

## fitting lasso 
lasso.mod <- glmnet(x, y, alpha = 1, lambda = grid, thresh = 1e-12)

## smallest lambda with all zero coefficients is at lambda = 0.5659 
## (i.e. this should be lambda_max)
lasso.mod 
## extracting relevant columns
lasso.df <- as.data.frame(do.call(cbind, lasso.mod[c("df", "lambda")]))
## extracting relevant rows
lasso.df[54:59, ]

## calculating the inner product |y^tx_j| for each column x_j
absinprod <- (y %*% x) / length(y) 

## according to the book it should be lambda_max=max_j |x_j^T*y|
max(absinprod) 



## ------------------------------------------------------------------------------------------------
## trying the next inequality 
## -----------------------------------------------------------------------------------------------

## trying to check |x_j^T*(y-y_lambda^hat)|=lambda for j in 
##   the active set and < for j not in the active set
## from Trevor Hastie, Robert Tibshirani, and 
## Martin Wainwright (2015), Statistical learning with sparsity: 
## the Lasso and generalizations, p127
grid <- 10 ^ seq(from = 4, to = -3, length = 100)  
lasso.mod <- glmnet(x, y, alpha = 1, lambda = grid, thresh = 1e-12)
lasso.mod # to choose a lambda
(lam <- lasso.mod$lambda[84]) # this is the chosen lambda

## predict to get y.hat
lasso.pred <- predict(lasso.mod, s = lam, newx = x) 
y.hat <- lasso.pred[, 1]
## predict to find active set
lasso.coef <- predict(lasso.mod, s = lam, type = "coefficients", newx = x) 

## calculating inner product 
absinprod2 <- apply(x, 2, function(k) abs(sum((y - y.hat) * k)) / length(y))

## lasso.coef - to see which predictors are in the active set
## absinprod2 - absolute inner product
## lam - chosen lambda
cbind(lasso.coef, round(c(0, absinprod2), 5), round(lam, 5))



