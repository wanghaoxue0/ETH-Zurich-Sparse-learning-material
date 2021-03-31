library(glmnet)

set.seed(1) 
n <- 600
p <- 1

Y <- matrix(rnorm(n * p), ncol = p)      # simulate observations
Y <- Y + rep(c(0, 1, 3), each = n/3)     # add changepoints
Y <- Y - mean(Y)                         # standardize for glmnet
Y <- Y / sqrt(mean(Y^2))

X  <- matrix(0, n, n)                     # define design matrix
X[lower.tri(X, diag = T)] <- 1
cm <- colMeans(X) 
X  <- t(t(X) - cm) 

lower <- 0.1; upper <- 1; stepsize <- 0.05;               # specify choices of regularization parameters
nfits <- (upper - lower)/stepsize + 1

RES           <- matrix(NA, ncol = 601, nrow = nfits)
RES_predicted <- matrix(NA, ncol = 600, nrow = nfits)

for(j in 1 : nrow(RES)){                                  # perform fit for each lambda
  fit <- glmnet(x = X, y = Y, alpha = 1,
               lambda = seq(lower, upper, stepsize)[j], standardize = F,
               penalty.factor = c(0, rep(1, n - 1)),      # note the 0, the first summand in the penalty is left unpenalized
               intercept = F, thresh = 10^-8)
  
  RES[j,]           <- as.numeric(coef(fit))
  RES_predicted[j,] <- predict(fit, newx = X)
}

RES_theta <- t(X %*% t(RES[, -1]))

plot(Y)
for(i in 1 : nrow(RES)){
   lines(RES_theta[i, ], col = i)
}
abline(v = c(200, 400), col = 2)                          # indicate the true changepoints

#######
# piecewise linear using trend-filtering
#######

D_n_1 		<- matrix(0, n - 1, n)
diag(D_n_1)	<- -1

for(i in 1 : (n - 1)){
   D_n_1[i, i + 1] <- 1
}

D			<- D_n_1[1:(n - 2), 1:(n - 1)] %*% diag(1, n - 1) %*% D_n_1

X	<- round(solve(rbind(c(1, rep(0, n - 1)), c(-1, 1, rep(0, n - 2)), D)))
cm  <- colMeans(X) 
X   <- t(t(X) - cm)

# takes a bit longer
lower <- 0.1; upper <- 1; stepsize <- 0.1;
nfits <- (upper - lower)/stepsize + 1

RES           <- matrix(NA, ncol = 601, nrow = nfits)
RES_predicted <- matrix(NA, ncol = 600, nrow = nfits)

for(j in 1:nrow(RES)){
  fit <- glmnet(x = X, y = Y, alpha = 1,
               lambda = seq(lower, upper, stepsize)[j], standardize = F,
               penalty.factor = c(0, 0, rep(1, n - 2)),   # note the difference when "0, 0" is put
                                                          # to "1, 1" for the first two entries
               intercept = F, thresh = 10^-8)
  RES[j,] <- as.numeric(coef(fit))
  RES_predicted[j,] <- predict(fit, newx = X)
}

RES_theta <- t(X %*% t(RES[, -1]))


plot(Y)
for(i in 1:nrow(RES)){
   lines(RES_theta[i, ], col = i)
}
abline(v = c(200, 400), col = 2)
