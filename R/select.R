select.grpreg <- function(obj,method="BIC",df="default",...)
  {
    if (df=="default") df <- obj$df
    else if (df=="active") df <- apply(obj$beta!=0,2,sum)
    else stop(paste("df",df,"is not an allowed choice"))
    IC <- criteria.grpreg(obj,method,df)
    i <- which.min(IC)
    if (min(obj$lambda) == obj$lambda[i]) warning(paste("minimum lambda selected for",obj$penalty))
    else if ((max(obj$lambda) == obj$lambda[i]) & obj$penalty=="gBridge") warning("maximum lambda selected")
    return(list(beta=obj$beta[,i],
                lambda=obj$lambda[i],
                df=df[i],
                IC=IC))
  }
select <- function(obj,...) UseMethod("select")
