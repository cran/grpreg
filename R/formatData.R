formatData <- function(Data)
  {
    if (is.null(Data$X)) stop("Must supply Data$X")
    if (is.null(Data$y)) stop("Must supply Data$y")
    if (is.null(Data$family)) Data$family <- "gaussian"
    if (!is.null(Data$group))
      {
        if (length(Data$group)!=ncol(Data$X)) stop("Data$group does not match Data$X")
      }
    else Data$group <- 1:ncol(Data$X)
    if (is.null(colnames(Data$X))) colnames(Data$X) <- paste("V",1:ncol(Data$X),sep="")
    Data$J <- max(Data$group)
    Data$K <- as.numeric(table(Data$group))
    if (!identical(as.integer(sort(unique(Data$group))),as.integer(1:Data$J))) stop("Groups must be consecutively numbered 1,2,3,...")
    
    Data$n <- nrow(Data$X)
    Data$meanx <- apply(Data$X,2,mean)
    Data$normx <- sqrt(apply((t(Data$X)-Data$meanx)^2,1,sum))/sqrt(Data$n)
    if (any(Data$normx < 0.0001)) stop("X contains columns which are numerically constant.  If you are attempting to specify an intercept, please remove these columns; an intercept is included automatically.")
    Data$X <- scale(Data$X,Data$meanx,Data$normx)
    Data$X <- cbind(1,Data$X)
    Data$group <- c(0,Data$group)
    colnames(Data$X)[1] <- "(Intercept)"
    Data$p <- ncol(Data$X)

    return(Data)
  }
