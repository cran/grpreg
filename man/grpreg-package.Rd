\name{grpreg-package}
\alias{grpreg-package}
\docType{package}
\title{Regularization paths for regression models with grouped covariates}
\description{
  This package fits regularization paths for linear or logistic
  regression models penalized by the group lasso, group bridge, or group
  MCP methods.  The algorithm is based on the idea of a locally
  approximated coordinate descent, and is stable and very fast.
  }
\details{
\tabular{ll}{
Package: \tab grpreg\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-11-11\cr
License: \tab GPL-2\cr
}
Accepts a list \code{Data} containing the response, design matrix,
family and covariate groupings, and produces the regularization path
over a grid of values for the tuning parameter \code{lambda}.  Also
provides methods for selecting the optimal point along the path using a
variety of information criteria and for plotting the paths.
}
\references{Breheny, P. and Huang, J. (2008) \emph{Penalized Methods for
    Bi-level variable selection}.  Tech report No. 393, Department of
  Statistics and Actuarial Science, University of Iowa.\url{http://www.stat.uiowa.edu/techrep/tr393.pdf}}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\keyword{models}
\keyword{regression}
\keyword{package}
\examples{
data(birthwt.grpreg)
Data.gaussian <- list(y=birthwt.grpreg$bwt,
                      X=as.matrix(birthwt.grpreg[,c(-1,-2)]),
                      family="gaussian",
                      group=c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8))
Data.binomial <- list(y=birthwt.grpreg$low,
                      X=as.matrix(birthwt.grpreg[,c(-1,-2)]),
                      family="binomial",
                      group=c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8))

fit1.gLasso <- grpreg(Data.gaussian,"gLasso")
fit1.gBridge <- grpreg(Data.gaussian,"gBridge",lambda.max=0.08)
fit1.gMCP <- grpreg(Data.gaussian,"gMCP")

fit2.gLasso <- grpreg(Data.binomial,"gLasso")
fit2.gBridge <- grpreg(Data.binomial,"gBridge",lambda.max=0.06)
fit2.gMCP <- grpreg(Data.binomial,"gMCP")

plot(fit1.gMCP)
select(fit2.gLasso)
}
