#library to compute gradient and Hessian
library("pracma")

#we use the following built-in dataset as an example
data("anscombe")

#RSS is the objective function (in terms of the coefficient)
fun <- function(v){
  return(sum((anscombe$y1-v[1]-v[2]*anscombe$y2)^2))
}

#calculating the real optimum by built-in linear regression
lin_mod <- lm(y1 ~ y2, anscombe)
opt = lin_mod[["coefficients"]]

#function of the gradient descent iteration
#arguments: v0 - starting point, tol - tolerance of error
grad_desc <- function(v0, tol = 10^-3){
  v = v0
  i = 1
  err = vector("numeric")
  
  #INITIAL STEP
  #calculating step size with bactracking rule
  s = 1
  while (fun(v - s*grad(fun, v)) > fun(v) - 0.5*s*dot(grad(fun, v), grad(fun, v))){
    s = 0.8*s
  }
  #gradient descent step
  v_prev = v
  v = v - s*grad(fun, v)
  #calculating error
  err[i] = abs(fun(v) - fun(v_prev))
  
  #ITERATING UNTIL ERROR IS UNDER THE REQUIRED TOLERANCE
  while (err[i] > tol){
    #calculating stepsize with bactracking rule
    s = 1
    while (fun(v - s*grad(fun, v)) > fun(v) - 0.5*s*dot(grad(fun, v), grad(fun, v))){
      s = 0.8*s
    }
    #gradient descent step
    v_prev = v
    v = v - s*grad(fun, v)
    #calculating error
    i = i + 1
    err[i] = abs(fun(v) - fun(v_prev))
  }
  return(list(v, err))
}

#function of the Newton's method
#arguments: v0 - starting point, tol - tolerance of error
newton <- function(v0, tol = 10^-3){
  v = v0
  i = 1
  err = vector("numeric")
  
  #INITIAL STEP
  #newton step
  delta = -(solve(hessian(fun, v)) %*% grad(fun, v))
  #calculating stepsize with bactracking rule
  s = 1
  while (fun(v + s*delta) > fun(v) + 0.5*s*dot(grad(fun, v), delta)){
    s = 0.8*s
  }
  v_prev = v
  v = v + s*delta
  #calculating error
  err[i] = abs(fun(v) - fun(v_prev))
  
  #ITERATING UNTIL ERROR IS UNDER THE REQUIRED TOLERANCE
  while (err[i] > tol){
    #newton step
    delta = -(solve(hessian(fun, v)) %*% grad(fun, v))
    #calculating stepsize with bactracking rule
    s = 1
    while (fun(v + s*delta) > fun(v) + 0.5*s*dot(grad(fun, v), delta)){
      s = 0.8*s
    }
    v_prev = v
    v = v + s*delta
    #calculating error
    i = i + 1
    err[i] = abs(fun(v) - fun(v_prev))
    
  }
  return(list(v, err))
}

#set the starting point of the iteration
v0 = c(12897, -3335)

result1 = grad_desc(v0, 10^-3)
result2 = newton(v0, 10^-3)

#plotting the logarithm of the error against the number of iterations
min1 = min(log(result1[[2]]))
min2 = min(log(result2[[2]]))
max1 = max(log(result1[[2]]))
max2 = max(log(result2[[2]]))
plot(seq(length(result1[[2]])), log(result1[[2]]), xlim = c(0, max(length(result1[[2]]), length(result2[[2]]))),
     type = "l", xlab = "number of iterations", ylab = "log(error)", main = "convergence rate",
     ylim = c(min(min1, min2), max(max1, max2)))
lines(seq(length(result2[[2]])), log(result2[[2]]), pch = 2, col = 2)
legend("topright", c("grad descent", "newton's method"), col = c(1, 2), lty = c(1, 1))