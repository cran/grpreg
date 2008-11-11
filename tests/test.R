library(grpreg)
X <- matrix(rnorm(500),ncol=10)
b <- rnorm(10)
y <- rnorm(X%*%b)
Data <- list(X=X,y=y,family="gaussian",group=c(1,1,1,2,2,2,3,3,3,3))
tol <- .01

coef <- lm(Data$y~Data$X)$coef
beta <- grpreg(Data,lambda=0,penalty="gLasso",eps=.0001,max.iter=500)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")
beta <- grpreg(Data,lambda=0,penalty="gBridge",eps=.0001,max.iter=500)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")
beta <- grpreg(Data,lambda=0,penalty="gMCP",eps=.0001,max.iter=500)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")

Data.std <- grpreg(Data,penalty="gMCP")$Data

coef <- lm(Data.std$y~0+Data.std$X)$coef
beta <- grpreg(Data,lambda=0,penalty="gLasso",eps=.0001,max.iter=500)$beta.std
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")
beta <- grpreg(Data,lambda=0,penalty="gBridge",eps=.0001,max.iter=500)$beta.std
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")
beta <- grpreg(Data,lambda=0,penalty="gMCP",eps=.0001,max.iter=500)$beta.std
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")

X <- matrix(rnorm(1000),ncol=10)
b <- rnorm(10)
y <- 1*(rnorm(X%*%b)>0)
Data <- list(X=X,y=y,family="binomial",group=c(1,1,1,1,2,2,3,3,3,3))
tol <- .01

coef <- glm(Data$y~Data$X,family="binomial")$coef
beta <- grpreg(Data,lambda=0,penalty="gLasso",eps=.0001,max.iter=500)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")
beta <- grpreg(Data,lambda=0,penalty="gBridge",eps=.0001,max.iter=500)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")
beta <- grpreg(Data,lambda=0,penalty="gMCP",eps=.0001,max.iter=500)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")

Data.std <- grpreg(Data,penalty="gMCP")$Data

coef <- glm(Data.std$y~0+Data.std$X,family="binomial")$coef
beta <- grpreg(Data,lambda=0,penalty="gLasso",eps=.0001,max.iter=500)$beta.std
if (max(abs(coef - beta)) > tol) stop("binomial check failed")
beta <- grpreg(Data,lambda=0,penalty="gBridge",eps=.0001,max.iter=500)$beta.std
if (max(abs(coef - beta)) > tol) stop("binomial check failed")
beta <- grpreg(Data,lambda=0,penalty="gMCP",eps=.0001,max.iter=500)$beta.std
if (max(abs(coef - beta)) > tol) stop("binomial check failed")

