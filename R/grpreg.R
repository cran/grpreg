grpreg <- function(X, y, group=1:ncol(X), family=c("gaussian","binomial"), penalty=c("gMCP","gBridge","gLasso"), nlambda=100, lambda, lambda.min=ifelse(n>p,.001,.05), lambda.max, alpha=.999, eps=.005, max.iter=1000, delta=1e-8, gamma=switch(penalty,gMCP=3,gBridge=0.5), verbose=FALSE, warn.conv=TRUE)
  {
    ## Check for errors
    family <- match.arg(family)
    if (length(group)!=ncol(X)) stop("group does not match X")
    if (is.null(colnames(X))) colnames(X) <- paste("V",1:ncol(X),sep="")
    if (delta <= 0) stop("Delta must be a positive number")
    J <- max(group)
    K <- as.numeric(table(group))
    if (!identical(as.integer(sort(unique(group))),as.integer(1:J))) stop("Groups must be consecutively numbered 1,2,3,...")

    ## Scale X
    n <- nrow(X)
    meanx <- apply(X,2,mean)
    normx <- sqrt(apply((t(X)-meanx)^2,1,sum))/sqrt(n)
    if (any(normx < 0.0001)) stop("X contains columns which are numerically constant.  If you are attempting to specify an intercept, please remove these columns; an intercept is included automatically.")
    XX <- cbind(1,scale(X,meanx,normx))
    group <- c(0,group)
    colnames(XX)[1] <- "(Intercept)"
    p <- ncol(XX)

    ## Setup lambda
    if (missing(lambda)) lambda <- setupLambda(XX,y,group,family,J,K,penalty,lambda.max,lambda.min,nlambda,gamma)
    l <- length(lambda)

    ## Fit
    path <- .C("gpPathFit",double(p*l),integer(l),double(l),as.double(XX),as.double(y),as.integer(group),family,as.integer(n),as.integer(p),as.integer(J),as.integer(K),penalty,as.double(lambda),as.integer(l),as.double(eps),as.integer(max.iter),as.integer(verbose),as.double(delta),as.double(gamma),as.double(alpha*lambda),as.integer(warn.conv))
    beta <- matrix(path[[1]],nrow=p,byrow=T,dimnames=list(colnames(XX),round(lambda,digits=4)))
    beta <- unstandardize(beta,meanx,normx)

    val <- list(beta=beta,
                family=family,
                group=group,
                lambda=lambda,
                alpha=alpha,
                loss = calcL(cbind(1,X),y,beta,family),
                n = length(y),
                penalty=penalty,
                df=path[[3]],
                iter=path[[2]])
    class(val) <- "grpreg"
    return(val)
  }
