
### Ridge Regression and the Lasso
# perform ridge regression and the lasso in order to predict Salary on the Hitters data

library(glmnet)      # function glmnet() can be used to fit ridge and lasso models

library(ISLR)
data(Hitters)
str(Hitters)
?Hitters
Hitters=na.omit(Hitters)      # remove missing values

# define x and y, as we cannot use y~x notation in glmnet
# x has to be a matrix and y a vector
x=model.matrix(Salary~.,Hitters)[,-1]          
y=Hitters$Salary

# Ridge Regression
grid=10^seq(from = 10, to = -2, length=100)  # grid for lambda values
ridge.mod=glmnet(x,y,alpha=0,lambda=grid)    # alpha argument determines what type of model is fit
                                             # alpha=0 corresponds to ridge
                     # associated with each value of lambda is a vector of ridge 
                     # regression coefficients
dim(coef(ridge.mod)) # each row corresponds to a predictor (plus intercept), 
                     # each column to a value of lambda


### for a large value of lambda, we exptect the coefficient estimates to be much smaller,
### in terms of l_2-norm, as compared to when a small value of lambda is used

ridge.mod$lambda[50] # =grid[50]: lambda = 11'498
coef(ridge.mod)[,50] # corresponding coefficient estimates
sqrt(sum(coef(ridge.mod)[-1,50]^2)) # corresponding l_2-norm

ridge.mod$lambda[60] # =grid[60]: lambda = 705
coef(ridge.mod)[,60]
sqrt(sum(coef(ridge.mod)[-1,60]^2))


# create training and test set  
set.seed(1)                          #set random seed s.t. result can be reproduced
train=sample(1:nrow(x), nrow(x)/2)   # indices for training data
test=(-train)                        # indices we *leave out* for test data
y.test=y[test]
x.test=x[test,]

# conduct ridge regression on training data
ridge.mod = glmnet(x[train,], y[train], alpha=0, lambda=grid, thresh=1e-12)

# choose lambda via cross validation (use built-in cross validation function cv.glmnet())
?cv.glmnet   # conducts 10 fold cross validation by default
set.seed(1)
cv.out=cv.glmnet(x[train,],y[train],alpha=0)
par(mfrow=c(1,1))
plot(cv.out)
bestlamda=cv.out$lambda.min
bestlamda  # 212

# evaluate test MSE for lambda=212
ridge.pred=predict(ridge.mod,s=bestlamda,newx=x[test,])
mean((ridge.pred-y.test)^2)       #96016

# refit on the full data 
out=glmnet(x,y,alpha=0)   # use default values for lambda
predict(out,type="coefficients",s=bestlamda)[1:20,]  # plug in lambda chosen by CV
# as expected, none of the coefficients are zero (ridge regression
# does not perform variable selection!)



### lasso

# By using alpha=1 we can do similar things for the lasso:
# Conduct lasso regression on training data
lasso.mod = glmnet(x[train,], y[train], alpha=1, lambda=grid, thresh=1e-12)
plot(lasso.mod)        # we can see from the coefficient plot that depending on lambda,
                       # some of the coefficients will be exactly equal to zero

cv.out = cv.glmnet(x[train,],y[train],alpha=1)     #perform CV on training set
par(mfrow=c(1,1))
plot(cv.out)
bestlam=cv.out$lambda.min
bestlam  # 35.3

# evaluate test MSE for lambda=35.3
lasso.pred=predict(lasso.mod,s=bestlam,newx=x[test,])
mean((lasso.pred-y.test)^2)         # 102'442
# test MSE is very similar to the test MSE of the ridge regression (96016)


# refit on the full data
out=glmnet(x,y,alpha=1)   # use default values for lambda
predict(out,type="coefficients",s=bestlam)[1:20,]  # plug in lambda chosen by CV

# lasso has substantial advantage over ridge regression because
# the resulting coefficients are sparse

