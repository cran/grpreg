plot.grpreg <- function(x, color=TRUE, alpha=1, type="l", pch=1:p, lty=1, xlab=expression(lambda), ylab=expression(hat(beta)), legend.loc, ...)
  {
    zeros <- which(apply(abs(x$beta),1,sum)==0)
    ind <- -1*c(1,zeros)
    beta <- x$beta[ind,,drop=FALSE]
    p <- nrow(beta)
    g1 <- as.factor(x$group[ind])
    g2 <- as.numeric(g1)
    n.g2 <- max(g2)
    if (color) col <- hsv(seq(0,1,len=(n.g2+1)),alpha=alpha)[g2]
    else col <- rgb(0,0,0,alpha=alpha)
    l <- x$lambda
    matplot(l,t(beta),type=type,xlab=xlab,ylab=ylab,pch=pch,col=col,lty=lty,xlim=rev(range(l)),...)
    abline(h=0)    
    if(!missing(legend.loc))
      {
        if (length(col)==1) col <- rep(col,p)
        legend(legend.loc,lty=lty,legend=g1[!duplicated(g2)],col=col[!duplicated(g2)])
      }
  }
