\name{grpreg}
\alias{grpreg}
\title{fit a group penalized regression path}
\description{Fit paths for group lasso, group bridge, or group MCP at
a grid of values of the regularization parameter lambda. Fits linear and
logistic regression models.}
\usage{
grpreg(X, y, group=1:ncol(X), family="gaussian", penalty="gMCP", lambda,
n.lambda=100, lambda.min=ifelse(n>p,.001,.05), lambda.max, lambda2=.001,
eps=.005, max.iter=1000, delta=1e-8, gamma=.5,
a=ifelse(family=="gaussian",3,30), verbose=FALSE, monitor=NULL,
warn.conv=TRUE)
}
\arguments{
  \item{X}{The design matrix, without an intercept.  \code{grpreg}
  standardizes the data and includes an intercept by default.}
  \item{y}{The response vector.}
  \item{group}{A vector of consecutive integers describing the grouping
    of the coefficients (see example below).}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{penalty}{The penalty to be applied to the model.  One of
    "gLasso" for group lasso, "gBridge" for group bridge, or "gMCP" for
    group MCP.}
  \item{lambda}{A user supplied sequence of \code{lambda} values.
  Typically, this is left unspecified, and the function automatically
  computes a vector of lambda values that ranges uniformly on the log 
  scale from \code{lambda.min} to \code{lambda.max}.}
  \item{n.lambda}{The number of \code{lambda} values.  Default is 100.}
  \item{lambda.min}{The smallest value for \code{lambda}, as a fraction
  of \code{lambda.max}.  Default is .001 if the number of observations
  is larger than the number of covariates and .05 otherwise.} 
  \item{lambda.max}{The largest value for \code{lambda}.  Default is to
  use the smallest value for which all penalized coefficients are 0.
  This can be calculated for group lasso and group MCP, but not for
  group bridge.  If \code{lambda.max} is left unspecified for group
  bridge, \code{grpreg} will make a rough guess.}
  \item{lambda2}{By default, a small L2 penalty is included alongside
  the group penalty.  \code{lamdba2} controls the magnitude of this
  penalty, as a fraction of \code{lambda}.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
  relative change in any coefficient is less than \code{eps}.  Default
  is \code{.005}.  See details.}
  \item{max.iter}{Maximum number of iterations.  Default is 1000.  See
    details.}
  \item{delta}{The amount by which the group lasso penalty is bounded
  away from 0.  Smaller values of \code{delta} produce more accurate
  results, but converge more slowly.  See example.}
  \item{gamma}{Tuning parameter of the bridge penalty; defaults to 1/2}
  \item{a}{Tuning parameter of the group MCP penalty - defaults to 3 for
  linear regression and 30 for logistic regression;}
  \item{warn.conv}{Should the function give a warning if it fails to
  converge?  Default is TRUE.  See details.} 
  \item{verbose}{Get a running update on what the algorithm is doing.
  Default is FALSE.} 
  \item{monitor}{Monitor the iterations of a vector of covariates.  If
  set to a numeric vector, for example \code{c(3,5)}, the algorithm will
  display iterates of the third and fifth covariates as it progresses.}
}
\details{
  The sequence of models indexed by \code{lambda} is fit using a locally
  approximated coordinate descent algorithm.  For logistic regression
  models, some care is taken to avoid model saturation; the algorithm
  may exit early in this setting.  The objective function is defined to
  be
  \deqn{1/(2*n)RSS + penalty}
  for \code{"gaussian"} and
  \deqn{-1/nobs loglik + \lambda*penalty}
  for \code{"binomial"}, where the likelihood is from a traditional
  generalized linear model for the log-odds of an event.\cr
  \cr
  This algorithm is stable and generally converges quite rapidly to
  values close to the solution.  However, especially when p is large
  compared with n, \code{grpreg} may fail to converge at low values of
  \code{lambda}, where models are nonidentifiable or nearly singular,
  Often, this is not the region of the coefficient path that is most
  interesting.  The default behavior warning the user when convergence
  criteria are not met may be distracting in these cases, and can be
  modified with \code{warn.conv} (convergence can always be checked
  later by inspecting the value of \code{iter}).\cr
  \cr
  If models are not converging, increasing \code{max.iter} may not be
  the most efficient way to correct this problem.  Consider increasing
  \code{eps}, increasing \code{n.lambda}, or increasing
  \code{lambda.min} in addition to increasing \code{max.iter}.
}
\value{
  An object with S3 class \code{"grpreg"} containing:
  \item{beta}{The fitted matrix of coefficients.  The number of rows is
    equal to the number of coefficients, and the number of columns is
    equal to \code{n.lambda}.}
  \item{family}{Same as above.}
  \item{lambda}{The sequence of \code{lambda} values in the path.}
  \item{lambda2}{Same as above.}
  \item{penalty}{Same as above.}
  \item{df}{A vector of length \code{n.lambda} containing estimates of
    effective number of model parameters all the points along the
    regularization path.  For details on how this is calculated, see
    reference.}
  \item{iter}{A vector of length \code{n.lambda} containing the number
    of iterations until convergence at each value of \code{lambda}.}
}
\references{Breheny, P. and Huang, J. (2009) Penalized methods for
  bi-level variable selection.  \emph{Statistics and its interface},
  \bold{2}: 369--380.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{plot} and \code{select} methods.}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="gLasso")
plot(fit)
fit <- grpreg(X,y,group,penalty="gBridge",lambda.max=.08)
plot(fit)
fit <- grpreg(X,y,group,penalty="gMCP")
plot(fit)
select(fit,X,y)

# Role of delta
# Note that fit2 cuts down on the # of iterations,
# but creates artifacts near 0
fit1 <- grpreg(X,y,group,penalty="gLasso")
fit2 <- grpreg(X,y,group,penalty="gLasso",delta=1e-4)
par(mfrow=c(2,1))
plot(fit1)
plot(fit2)
dev.off()
plot(fit1$lambda,fit1$iter,type="l")
lines(fit2$lambda,fit2$iter,lty=2)

# Logistic regression
y <- birthwt.grpreg$low
fit <- grpreg(X,y,group,penalty="gLasso")
plot(fit)
fit <- grpreg(X,y,group,penalty="gBridge",lambda.max=.03)
plot(fit)
fit <- grpreg(X,y,group,penalty="gMCP")
plot(fit)
select(fit,X,y)
}
\keyword{models}
\keyword{regression}
