criteria.grpreg <- function(fit,method="BIC",df)
  {
    n <- length(fit$lambda)
    L <- numeric(n)
    for (i in 1:n) L[i] <- calcL(fit$Data,fit$beta[,i])

    if (method=="AIC") return(2*L + 2*df)
    if (method=="BIC") return(2*L + log(fit$Data$n)*df)
    if (method=="GCV") return(2*L/(1-df/fit$Data$n)^2)
    stop(paste("method",method,"not recognized"))
  }
