calcL <- function(Data,beta)
  {
    if (Data$family=="gaussian") return(apply((Data$y - Data$X %*% beta)^2,2,sum)/2)
    if (Data$family=="binomial")
      {
        eta <- Data$X %*% beta
        pi. <- exp(eta)/(1+exp(eta))
        return(-sum(Data$y*log(pi.)+(1-Data$y)*log(1-pi.)))
      }
  }
