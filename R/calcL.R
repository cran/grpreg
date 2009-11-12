calcL <- function(X,y,beta,family)
  {
    if (family=="gaussian") return(apply((y - X %*% beta)^2,2,sum)/2)
    if (family=="binomial")
      {
        eta <- X %*% beta
        pi. <- exp(eta)/(1+exp(eta))
        return(-sum(y*log(pi.)+(1-y)*log(1-pi.)))
      }
  }
