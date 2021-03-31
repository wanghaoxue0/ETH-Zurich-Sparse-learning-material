install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)
data(QuickStartExample)
fit = glmnet(x, y)
#default 
plot(fit)
#by default does lasso
cvfit = cv.glmnet(x, y)
plot(cvfit)
#It includes the cross-validation curve (red dotted line), and upper and lower standard deviation curves along the lamda sequence (error bars). Two selected lamda¡¯s are indicated by the vertical dotted lines
#lambda.min is the value of lamda that gives minimum mean cross-validated error. The other lamda saved is lambda.1se, which gives the most regularized model such that error is within one standard error of the minimum

fit2= glmnet(x, y, alpha = 0.2, weights = c(rep(1,50),rep(2,50)), nlambda = 20)
plot(fit2)
#As an example, we set alpha = 0.2 (more like a ridge regression), and give double weights to the latter half of the observations.

cv1=cv.glmnet(x,y,alpha=1)
cv.5=cv.glmnet(x,y,alpha=.5)
cv0=cv.glmnet(x,y,alpha=0)
par(mfrow=c(2,2))
plot(cv1);plot(cv.5);plot(cv0)

par(mfrow=c(1,1))
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name)
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))
#we see that lasso (alpha=1) does about the best here for small lamda. We also see that the range of lambdas used differs with alpha.

#Coefficient upper and lower bounds
#Suppose we want to fit our model,but limit the coefficients to be bigger than -0.7 and less than 0.5. This is easily achieved via the upper.limits and lower.limits arguments:
tfit=glmnet(x,y,lower=-.7,upper=.5)
plot(tfit)

#panalty factors 
#This argument allows users to apply separate penalty factors to each coefficient. Its default is 1 for each parameter, but other values can be specified. In particular, any variable with penalty.factor equal to zero is not penalized at all!
p.fac = rep(1, 20)
p.fac[c(5, 10, 15)] = 0
pfit = glmnet(x, y, penalty.factor = p.fac)
plot(pfit, label = TRUE)
#we need to kept the 5th, 10th and 15th variable all the time
#We see from the labels that the three variables with 0 penalty factors always stay in the model, while the others follow typical regularization paths and shrunken to 0 eventually.

#now lets see a different family!!!
#Binomial family allows logistic regression 
data(BinomialExample)
#the response variable y should be either a factor with two levels, or a two-column matrix of counts or proportions.
bfit = glmnet(x, y, family = "binomial")
predict(bfit, newx = x[1:5,], type = "class", s = c(0.05, 0.01))
#we make prediction of the class labels at lamda= 0.05, 0.01.
#class¡± produces the class label corresponding to the maximum probability
#¡°response¡± gives the fitted probabilities

predict(bfit, newx = x[1:5,], type = "response", s = c(0.05, 0.01))

cvfitb = cv.glmnet(x, y, family = "binomial", type.measure = "class")
plot(cvfitb)
