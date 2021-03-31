library(ggplot2)
library(reshape2)
library(glmnet)
library(expm)
library(glasso)

## Graphical Lasso (Algorithm 9.1 in the book) for a range of different values for lambda
#using the non-discrete variables of the diabetes dataset:
#https://web.stanford.edu/~hastie/StatLearnSparsity/data.html

predictors <- diabetes[,c(3,4,6,8,9)]
predictors <- scale(predictors)
S <- cov(predictors)
p <- 5


#repeat the algorithm for different values of lambda
x=seq(from=0, to=1.2, by=0.02)
for (lambda in x)
{
  
  iter <- 10
  # 1. Initialize W==S:
  W = S + lambda*diag(p)  
  # 2. Repeat until convergence:
  
  
  
  for (i in 1:iter) 
  {
    for  (j in 1:p)                  
    {
      ## partition the matrices:
      ## P.1: all but j-th row and column:
      S11 <- S[-j,-j]
      W11 <- W[-j,-j]
      #Theta11 <- Theta[-j,-j]
      ## P.2: j-th row and column:
      s12 <- S[-j,j]
      w12 <- W[-j,j]
      #theta12 <- Theta[-j,j]
      s22 <- S[j,j]
      w22 <- W[j,j]
      #theta22 <- Theta[j,j]
      ## optimization (using package glmnet):
      fit <- glmnet(-sqrtm(W11), -solve(sqrtm(W11))%*%s12, family = "gaussian", alpha = 1, 
                    ## devide by (p-1):  
                    lambda =  lambda/(p-1),
                    standardize = FALSE, intercept = FALSE, thresh = 1e-5)
      beta <- coef(fit)[-1]
      ## update:
      W[-j,j]  <- W11%*%beta
      W[j,-j] <- W11%*%beta
      
    }
  }
  Theta <- solve(W)
  # for each lambda save the entries of the upper triangle of the estimated precision matrix 
  if (exists("coeff") == FALSE) 
  {
    coeff <- c(log(lambda), Theta[1,c(2,3,4,5)],Theta[2,c(3,4,5)],Theta[3,c(4,5)],Theta[4,5])
  } 
  else
  {
    coeff <- rbind(coeff, c(log(lambda), Theta[1,c(2,3,4,5)],Theta[2,c(3,4,5)],Theta[3,c(4,5)],Theta[4,5]))                
  }
}


## plot the coefficients vs. log-lambda
colnames(coeff) <- c("logl", "BMI-BP", "BMI-S2", "BMI-S4","BMI-S5","BP-S2","BP-S4","BP-S5","S2-S4","S2-S5","S4-S5") 
rownames(coeff) <- c()
mltd <- melt(coeff[, -1], id = "logl")
mltd <- cbind(coeff[, 1], mltd)

ggplot(mltd, aes(x = mltd[, 1], y = value, colour = Var2)) + 
  geom_line() + 
  ylab(label="Coefficients vs log-lambda") + 
  xlab("log-lambda")




