set.seed(0)
options(warn = -1)

# generate simple model with only intercept
model_intercept <- list(
  x = matrix(rep(1, 1000), ncol = 1),
  beta = list(
    '1' = 1,
    '2' = -1,
    '3' = 0
  ),
  alpha = c(150, 700),
  sigma = 0.1
)

# help function for later
r <- function(alpha, j){
        alpha[j + 1] - alpha[j]
}

# function to simulate y values from a model
SimulateModel <- function(model){

  n     <- nrow(model$x)
  p     <- ncol(model$x)

  alpha <- c(1, model$alpha, n + 1)
  k     <- length(alpha) - 1
  y     <- rep(NA, n)

  for (j in 1:k){
    x     <- model$x[alpha[j] : (alpha[j + 1] - 1), ]
    beta  <- as.matrix(model$beta[[as.character(j)]])
    y[alpha[j] : (alpha[j + 1] - 1)] <- x %*% beta + model$sigma * rnorm(r(alpha, j))
  }
  y
}

# generate y values for simple model
model_intercept$y <- SimulateModel(model_intercept)

# function to plot model
plot_cp <- function(y, alpha, main = 'Time Series with change points', ylab = 'y', plot_max = FALSE){
  plot(y, main = main, type='n', ylab = ylab, xlab = 'time')
  lines(y)
  abline(v = alpha, col='red')
  if (plot_max){
    points(which.max(y), max(y), cex = 2, pch = 4)
  }
}

# plot simple model
plot_cp(model_intercept$y, model_intercept$alpha)

# simple loss function
MSE_mean <- function(x, y, n_obs, ...){
                 sum((y - mean(y))^2) / n_obs
}

# show fit in plot & calculate MSE
abline(h = mean(model_intercept$y), col = 'green')
MSE_mean(model_intercept$x, model_intercept$y, 1000)

# function to calculate decrease in Loss if x is split at i
SplitLoss <- function(x, y, FUN, ...){
  n_obs <- nrow(x)

  SplitLossFun <- function(i){
    FUN(x[1:i, ], y[1:i], n_obs, ...) + FUN(x[(i + 1) : n_obs, ], y[(i + 1) : n_obs], n_obs, ...)
  }

  SSE_0 <- FUN(x, y, n_obs, ...)
  SSE   <- apply(matrix(2 : (n_obs-2)), 1, SplitLossFun)
  c(0, 0, SSE_0 - SSE, 0)
}

# evaluate above function
model_intercept$SplitLoss <- SplitLoss(model_intercept$x, model_intercept$y, MSE_mean)

#plot
x1 = 200
dev.off()
par(mfrow = c(2,1))
plot_cp(model_intercept$y, model_intercept$alpha)
lines(c(x1 + 1, 1000), c(mean(model_intercept$y[(x1 + 1) : 1000]), mean(model_intercept$y[(x1 + 1) : 1000])), col = 'green')
lines(c(1, x1), c(mean(model_intercept$y[1 : x1]), mean(model_intercept$y[1 : x1])), col = 'green')
abline(v = x1, col = 'green', lty = 'dashed')

plot_cp(model_intercept$SplitLoss, model_intercept$alpha, main = 'SplitLoss', ylab = 'delta Loss')
lines(c(x1,x1), c(-1, model_intercept$SplitLoss[x1]), col = 'green', lty = 'dashed')


# change model to include more error
model_intercept$sigma     <- 3
model_intercept$y         <- SimulateModel(model_intercept)
model_intercept$SplitLoss <- SplitLoss(model_intercept$x, model_intercept$y, MSE_mean)

#plot again
dev.off()
par(mfrow = c(2, 1))
plot_cp(model_intercept$y, model_intercept$alpha)
plot_cp(model_intercept$SplitLoss, model_intercept$alpha, main='SplitLoss', plot_max = T)


# now high dimensional model
p <- 400
n <- 400

model_hd <- list(
  x       = matrix(rnorm(n*p), ncol=p ),
  beta    = list(
      '1' = c(rnorm(10), rep(0, p - 10)),
      '2' = c(rep(0,5), rnorm(10), rep(0, p - 15)),
      '3' = c(rnorm(5), rep(0,5), rnorm(5), rep(0,p - 15))
  ),
  alpha = c(floor(n / 4), floor(2 * n / 3)),
  sigma = 0.1
)

dev.off()
par(mfrow=c(3, 1))
plot(1:15, model_hd$beta$'1'[1:15], main = 'beta1',  pch = 16, cex = 2); abline(v = c(5.5, 10.5))
plot(1:15, model_hd$beta$'2'[1:15], main = 'beta2',  pch = 16, cex = 2); abline(v = c(5.5, 10.5))
plot(1:15, model_hd$beta$'3'[1:15], main = 'beta3',  pch = 16, cex = 2); abline(v = c(5.5, 10.5))


library(glmnet)

#function for high dimensional loss via LASSO / glmnet
MSE_hd <- function(x, y, n_obs, ...){
  fit <- glmnet(x, y, lambda = sqrt(n_obs / nrow(x)) * list(...)$lambda)
  sum((y - x%*%fit$beta)^2) / n_obs
}

#get y values
model_hd$y <- SimulateModel(model_hd)

#get sensible lambda
lambda <- cv.glmnet(model_hd$x, model_hd$y)$lambda.min

model_hd$SplitLoss <- SplitLoss(model_hd$x, model_hd$y, MSE_hd, lambda = lambda)

dev.off()
plot_cp(model_hd$SplitLoss, model_hd$alpha, plot_max = T)

# what happens if we choose lambda incorrectly??
model_hd$SplitLoss_x5  <- SplitLoss(model_hd$x, model_hd$y, MSE_hd, lambda = 5   * lambda)
model_hd$SplitLoss_d5  <- SplitLoss(model_hd$x, model_hd$y, MSE_hd, lambda = 1/5 * lambda)

dev.off()
par(mfrow=c(3,1))
plot_cp(model_hd$SplitLoss_x5, model_hd$alpha, main='Time Series with change points, lambda high',    plot_max = T)
plot_cp(model_hd$SplitLoss,    model_hd$alpha, main='Time Series with change points, lambda optimal', plot_max = T)
plot_cp(model_hd$SplitLoss_d5, model_hd$alpha, main='Time Series with change points, lambda low',     plot_max = T)
