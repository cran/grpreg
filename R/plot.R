plot.grpreg <- function(x,col=NULL,alpha=FALSE,color=TRUE,type="l",lty=1,pch=NULL,xlab=expression(lambda),ylab=expression(hat(beta)),...)
  {
    if (is.logical(alpha))
      {
        if (alpha)
          {
            if (color) alpha <- exp(1-log(sum(x$Data$group!=0),14))
            else  alpha <- exp(1-log(sum(x$Data$group!=0),7))
          }
        else alpha <- 1
      }
    if (!is.numeric(alpha)) stop("alpha must be either logical or numeric")
    if (color)
      {
        if (is.null(col))
          {
            group <- x$Data$group[x$Data$group!=0]
            n.g <- max(x$Data$group)
            col <- rep(hcl(seq(0,360,len=(n.g+1)),c=1000,l=100,alpha=alpha)[1:n.g],table(group))
          }
      }
    else
      {
        col <- rgb(0,0,0,alpha=alpha)
      }
    if (is.null(pch))
      {
        p <- nrow(x$beta)
        pch <- 1:p
      }
    matplot(x$lambda,t(x$beta[x$Data$group!=0,]),type=type,xlab=xlab,ylab=ylab,pch=pch,col=col,lty=lty,...)
  }
