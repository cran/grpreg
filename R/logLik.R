logLik.grpreg <- function(object, df.method=c("default","active"), REML=FALSE, ...) {
  df.method <- match.arg(df.method)
  n <- as.integer(object$n)
  df <- if (df.method=="active") apply(coef(object)!=0, 2, sum) else object$df
  if (object$family=="gaussian") {
    rdf <- if (REML) n-df else n
    RSS <- object$loss
    l <- -n/2 * (log(2*pi) + log(RSS) - log(rdf)) - rdf/2
    df <- df + 1
  } else if (object$family=='poisson') {
    y <- object$y
    ind <- y != 0
    l <- -object$loss/2 + sum(y[ind]*log(y[ind])) - sum(y) - sum(lfactorial(y))
  } else {
    l <- -object$loss/2
  }
  structure(l, df=df, nobs=n, class='logLik')
}
logLik.grpsurv <- function(object, df.method=c("default","active"), ...) {
  df.method <- match.arg(df.method)
  n <- as.integer(object$n)
  df <- if (df.method=="active") apply(coef(object)!=0, 2, sum) else object$df
  structure(-object$loss/2, df=df, nobs=n, class='logLik')
}
