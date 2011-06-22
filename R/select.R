select.grpreg <- function(obj, criterion=c("BIC","AIC","GCV"), df.method=c("default","active"), ...)
  {
    criterion <- match.arg(criterion)
    df.method <- match.arg(df.method)
    ll <- logLik(obj,df.method=df.method)
    df <- as.numeric(attr(ll,"df"))
    
    if (criterion=="AIC") IC <- AIC(ll)
    if (criterion=="BIC") IC <- BIC(ll)
    if (criterion=="GCV") IC <- (1/obj$n) * (-2) * as.numeric(ll) / (1-df/obj$n)^2
    i <- which.min(IC)

    if (min(obj$lambda) == obj$lambda[i]) warning(paste("minimum lambda selected for",obj$penalty))
    else if ((max(obj$lambda) == obj$lambda[i]) & obj$penalty=="gBridge") warning("maximum lambda selected")
    return(list(beta=obj$beta[,i],
                lambda=obj$lambda[i],
                df=df[i],
                IC=IC))
  }
select <- function(obj,...) UseMethod("select")
