criteria.grpreg <- function(X,y,beta,family,method,df)
  {
    l <- ncol(beta)
    n <- length(y)
    L <- numeric(l)
    for (i in 1:l) L[i] <- calcL(X,y,beta[,i],family)

    if (method=="AIC") return(2*L + 2*df)
    if (method=="BIC") return(2*L + log(n)*df)
    if (method=="GCV") return(2*L/(1-df/n)^2)
    stop(paste("method",method,"not recognized"))
  }
