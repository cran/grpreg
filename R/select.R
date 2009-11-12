select.grpreg <- function(obj, X, y, criterion="BIC", df.method="default", ...)
  {
    if (df.method=="default") df <- obj$df
    else if (df.method=="active") df <- apply(obj$beta!=0,2,sum)
    else stop(paste("df.method",df.method,"is not an allowed choice\n"))
    IC <- criteria.grpreg(cbind(1,X),y,obj$beta,obj$family,criterion,df)
    i <- which.min(IC)
    if (min(obj$lambda) == obj$lambda[i]) warning(paste("minimum lambda selected for",obj$penalty))
    else if ((max(obj$lambda) == obj$lambda[i]) & obj$penalty=="gBridge") warning("maximum lambda selected")
    return(list(beta=obj$beta[,i],
                lambda=obj$lambda[i],
                df=df[i],
                IC=IC))
  }
select <- function(obj,...) UseMethod("select")
