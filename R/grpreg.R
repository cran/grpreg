## User-level function for fitting group penalization models
## This function:
##   1. Formats the data and options
##   2. For a vector of regulation parameters:
##      a. Obtains initial values by calling getInitial()
##      b. Fits the model by calling .C("gpFit")
##   3. Returns a list containing:
##      a. The formatted data and options
##      b. A matrix of fitted coefficients (ncol = length(lambda))
##         i. Standardized
##         ii. Transformed back to original scale
##      c. The number of iterations until convergence for each lambda value

grpreg <- function(Data,penalty,lambda=NULL,n.lambda=100,lambda.min=NULL,lambda.max=NULL,lambda2=.001,eps=.005,max.iter=100,verbose=FALSE,monitor=NULL,warn.conv=TRUE,...)
  {
    ## Format
    par <- list(penalty=penalty,lambda=lambda,n.lambda=n.lambda,lambda.min=lambda.min,lambda.max=lambda.max,lambda2=lambda2,eps=eps,max.iter=max.iter,verbose=verbose,monitor=monitor,...)
    Data <- formatData(Data)
    par <- formatPar(Data,par)

    path <- .C("gpPathFit",double(Data$p*par$n.lambda),integer(par$n.lambda),double(par$n.lambda),as.double(Data$X),as.double(Data$y),as.integer(Data$group),Data$family,as.integer(Data$n),as.integer(Data$p),as.integer(Data$J),as.integer(Data$K),par$penalty,as.double(par$lambda),as.integer(par$n.lambda),as.double(par$eps),as.integer(par$max.iter),as.integer(par$verbose),as.integer(par$monitor-1),as.integer(length(par$monitor)),as.double(par$penpars),as.double(par$lambda2),as.integer(warn.conv))
    beta <- matrix(path[[1]],nrow=Data$p,byrow=T,dimnames=list(colnames(Data$X),round(par$lambda,digits=3)))
    val <- list(beta=unstandardize(beta,Data),
                beta.std=beta,
                lambda=par$lambda,
                lambda2=par$lambda2[1]/par$lambda[1],
                penalty=par$penalty,
                df=path[[3]],
                iter=path[[2]],
                par=par[c("n.lambda","lambda.min","lambda.max","eps","delta","max.iter")],
                Data=Data)
    class(val) <- "grpreg"
    return(val)
  }
