formatPar <- function(Data,par)
  {
    if (par$penalty=="gLasso")
      {
        if (is.null(par$gamma)) par$penpars <- par$delta <- .0005
        else par$penpars <- par$delta
        if (par$delta<0) stop("Delta must be a positive number")
      }
    if (par$penalty=="gBridge")
      {
        if (!is.null(par$gamma)) par$penpars <- par$gamma
        else par$penpars <- par$gamma <- 1/2
      }
    if (par$penalty=="gMCP")
      {
        if (!is.null(par$a)) par$penpars <- par$a
        else if (Data$family=="gaussian") par$penpars <- par$a <- 3
        else if (Data$family=="binomial") par$penpars <- par$a <- 30
      }
    if (is.null(par$lambda))
      {
        if (is.null(par$lambda.max))
          {
            if (sum(Data$group==0)>0) fit <- glm(Data$y~0+Data$X[,Data$group==0],family=Data$family)
            else fit <- NULL
            if (Data$family=="gaussian")
              {
                if (is.null(fit)) r <- Data$y
                else r <- fit$residuals
                w <- rep(1,Data$n)
              }
            if (Data$family=="binomial")
              {
                if (is.null(fit)) eta <- rep(0,Data$n)
                else eta <- Data$X[,Data$group==0]%*%matrix(fit$coef,ncol=1)
                pi. <- exp(eta)/(1+exp(eta))
                w <- sqrt(pi.*(1-pi.))
                r = (Data$y - pi.)/w;
              }
            par$lambda.max <- .C("determineMax",double(1),as.double(Data$X[,Data$group!=0]),as.double(r),as.double(w),as.integer(Data$group),Data$family,as.integer(Data$n),as.integer(Data$p),as.integer(Data$J),as.integer(Data$K),par$penalty,as.double(par$penpars))[[1]]
          }
        if (is.null(par$lambda.min))
          {
            if (Data$n > Data$p) par$lambda.min <- .001*par$lambda.max
            else par$lambda.min <- .05*par$lambda.max
            par$lambda <- exp(seq(log(par$lambda.max),log(par$lambda.min),len=par$n.lambda))
          }
        else if (par$lambda.min==0)
          {
            if (Data$n > Data$p) par$lambda.min <- .001*par$lambda.max
            else par$lambda.min <- .05*par$lambda.max
            par$lambda <- c(exp(seq(log(par$lambda.max),log(par$lambda.min),len=par$n.lambda-1)),0)
          }
        else par$lambda <- exp(seq(log(par$lambda.max),log(par$lambda.min*par$lambda.max),len=par$n.lambda))
      }
    else par$n.lambda <- length(par$lambda)
    if (par$penalty=="gBridge") par$lambda <- rev(par$lambda)
    par$lambda2 <- par$lambda*par$lambda2
    return(par)
  }
