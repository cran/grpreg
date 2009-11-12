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

grpreg <- function(X, y, group=1:ncol(X), family="gaussian", penalty="gMCP", lambda, n.lambda=100, lambda.min=ifelse(n>p,.001,.05), lambda.max, lambda2=.001, eps=.005, max.iter=1000, delta=1e-8, gamma=.5, a=ifelse(family=="gaussian",3,30), verbose=FALSE, monitor=NULL, warn.conv=TRUE)
  {
    ## Setup/error check
    if (length(group)!=ncol(X)) stop("group does not match X")
    if (is.null(colnames(X))) colnames(X) <- paste("V",1:ncol(X),sep="")
    if (delta <= 0) stop("Delta must be a positive number")
    J <- max(group)
    K <- as.numeric(table(group))
    if (!identical(as.integer(sort(unique(group))),as.integer(1:J))) stop("Groups must be consecutively numbered 1,2,3,...")
    n <- nrow(X)
    meanx <- apply(X,2,mean)
    normx <- sqrt(apply((t(X)-meanx)^2,1,sum))/sqrt(n)
    if (any(normx < 0.0001)) stop("X contains columns which are numerically constant.  If you are attempting to specify an intercept, please remove these columns; an intercept is included automatically.")
    X <- scale(X,meanx,normx)
    X <- cbind(1,X)
    group <- c(0,group)
    colnames(X)[1] <- "(Intercept)"
    p <- ncol(X)
    if (missing(lambda)) lambda <- setupLambda(X,y,group,family,J,K,penalty,lambda.max,lambda.min,n.lambda,a,gamma)
    l <- length(lambda)

    path <- .C("gpPathFit",double(p*l),integer(l),double(l),as.double(X),as.double(y),as.integer(group),family,as.integer(n),as.integer(p),as.integer(J),as.integer(K),penalty,as.double(lambda),as.integer(l),as.double(eps),as.integer(max.iter),as.integer(verbose),as.integer(monitor-1),as.integer(length(monitor)),as.double(delta),as.double(gamma),as.double(a),as.double(lambda2*lambda),as.integer(warn.conv))
    beta <- matrix(path[[1]],nrow=p,byrow=T,dimnames=list(colnames(X),round(lambda,digits=4)))

    val <- list(beta=unstandardize(beta,meanx,normx),
                family=family,
                group=group,
                lambda=lambda,
                lambda2=lambda2,
                penalty=penalty,
                df=path[[3]],
                iter=path[[2]])
    class(val) <- "grpreg"
    return(val)
  }
