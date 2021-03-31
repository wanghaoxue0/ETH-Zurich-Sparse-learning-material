#Solving the Graph Selection Problem for a multivariate Gaussian distribution
#using the non-discrete variables of the diabetes dataset:
#https://web.stanford.edu/~hastie/StatLearnSparsity/data.html

#Grraphical Lasso (Algorithm 9.1 in the book)


#######using the glasso package#########
library(glasso)
library(igraph)
library(glmnet)
library(expm)


predictors <- diabetes[,c(3,4,6,8,9)]

#scale to get zero mean distribution
predictors <- scale(predictors)
N <-nrow(predictors)
p <-ncol(predictors)

#Initial matrix for glasso algorithm
S <- cov(predictors)

#regularization parameter (page 252 in the book)
lambda <- 2*sqrt(log(p)/N)

Theta <- glasso(S,lambda)
colnames(Theta$wi) <- c("BMI", "BP", "S2", "S4","S5")
Theta$wi

#plot the underlying graph of the distribution
G <- graph_from_adjacency_matrix(Theta$wi, mode="min",weighted=TRUE)
G <- simplify(G)
plot(G)


#####manual graphical lasso######
lambda <- 0.4
iter <- 20
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
    
    ## P.2: j-th row and column:
    s12 <- S[-j,j]
    w12 <- W[-j,j]
    
    s22 <- S[j,j]
    w22 <- W[j,j]
    
    ## optimization (using package glmnet):
    fit <- glmnet(-sqrtm(W11), -solve(sqrtm(W11))%*%s12, family = "gaussian", alpha = 1, 
                  ## devide by (p-1):  
                  lambda =  lambda/(p-1),
                  standardize = FALSE, intercept = FALSE, thresh = 1e-12)
    beta <- coef(fit)[-1]
    ## update:
    W[-j,j]  <- W11%*%beta
    W[j,-j] <- W11%*%beta
  }
}  

## plot the underlying graph
Theta <- solve(W)
G1 <- graph_from_adjacency_matrix(Theta, mode="min",weighted=TRUE)
G1 <- simplify(G1)
plot(G1)


#Neighborhood based (Meinshausen,Bühlmann) for gaussian graphical model (Algorithm 9.2 in the book)
lambda <- 2*sqrt(log(p)/N)
Theta2 <- glasso(S,lambda,approx=TRUE)
colnames(Theta2$wi) <- c("BMI", "BP", "S2", "S4","S5")
Theta2$wi

#plot the underlying graph of the distribution
G2 <- graph_from_adjacency_matrix(Theta2$wi,mode="min",weighted=TRUE)
G2 <- simplify(G2)
plot(G2)

