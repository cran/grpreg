unstandardize <- function(beta,Data)
  {
    val <- beta
    val[1,] <- beta[1,]-apply(Data$meanx*beta[-1,,drop=F]/Data$normx,2,sum)
    val[-1,] <- beta[-1,]/Data$normx
    return(val)
  }
