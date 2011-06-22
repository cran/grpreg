setupLambda <- function(X,y,group,family,J,K,penalty,lambda.max,lambda.min,nlambda,gamma)
  {
    if (missing(lambda.max))
      {
        fit <- glm(y~0+X[,group==0],family=family)
        if (family=="gaussian")
          {
            r <- fit$residuals
            w <- rep(1,length(y))
          }
        if (family=="binomial")
          {
            eta <- X[,group==0]%*%matrix(fit$coef,ncol=1)
            pi. <- exp(eta)/(1+exp(eta))
            w <- sqrt(pi.*(1-pi.))
            r = (y - pi.)/w;
          }
        lambda.max <- .C("determineMax",double(1),as.double(X[,group!=0]),as.double(r),as.double(w),as.integer(group),family,as.integer(length(y)),as.integer(sum(group!=0)),as.integer(J),as.integer(K),penalty,gamma)[[1]] + 1e-6
      }
    if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
    else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
    if (penalty=="gBridge") lambda <- rev(lambda)
    return(lambda)
  }
