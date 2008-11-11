\name{grpreg}
\alias{grpreg}
\title{fit a group penalized regression path}
\description{Fit paths for group lasso, group bridge, or group MCP at
a grid of values of the regularization parameter lambda. Fits linear and
logistic regression models.}
\usage{
grpreg(Data,penalty,lambda=NULL,n.lambda=100,lambda.min=NULL,lambda.max=NULL,lambda2=.001,eps=.005,max.iter=100,verbose=FALSE,monitor=NULL,warn.conv=TRUE,...)
}

\arguments{
  \item{Data}{A list containing: \code{y}, the response vector;
    \code{X}, the design matrix; \code{family}, either "gaussian" or
    "binomial", depending on the response; and \code{group}, a vector
    of consecutive integers describing the grouping of the coefficients
    (see example below). The design matrix should not contain an
    intercept; \code{grpreg} standardizes the data and includes an
    intercept by default.}
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
    have the algorithm calculate the smallest value for which all
    penalized coefficients are 0.  This can be done exactly for group
    lasso and group MCP, but only guessed at for group bridge.}
  \item{lambda2}{By default, a small L2 penalty is included alongside
    the group penalty.  \code{lamdba2} controls the magnitude of this
    penalty, as a fraction of \code{lambda}.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{.005}.  See details.}
  \item{max.iter}{Maximum number of iterations.  Default is 100.  See
    details.}
  \item{warn.conv}{Should the function give a warning if it fails to
    converge?  Default is TRUE.  See details.}
  \item{verbose}{Get a running update on what the algorithm is doing.
    Default is FALSE.}
  \item{monitor}{Monitor the iterations of a vector of covariates.  If
    set to a numeric vector, for example \code{c(3,5)}, the algorithm
    will display iterates of the third and fifth covariates as it progresses.}
  \item{...}{Other parameters specific to one of the penalized method.
    These include \code{delta}, the amount by which the group lasso
    penalty is bounded away from 0 - defaults to \code{0.0005};
    \code{a}, the tuning parameter of the group MCP penalty - defaults
    to 3 for linear regression and 30 for logistic regression; and
    \code{gamma}, the tuning parameter of the bridge penalty- defaults
    to 1/2.}
}
\details{
  The sequence of models indexed by \code{lambda} is fit using a locally
  approximated coordinate descent algorithm.  For logistic regression
  models, some care is taken to avoid overfitting to unstable models.
  The algorithm may exit early and give unstable results in this
  setting.  The objective function is defined to be \deqn{1/(2*n)RSS +
    penalty} for \code{"gaussian"} and \deqn{-1/nobs loglik +
    \lambda*penalty} for \code{"binomial"}, where the likelihood is from
  a traditional generalized linear model for the log-odds of an
  event.\cr
  \cr
  This algorithm is stable and rapidly converges to values close to the
  solution.  However, it also displays a linear rate of convergence
  (Newton-Raphson algorithms, in contrast, display quadratic rates of
  convergence), meaning that it may take a large number of iterations to
  reach accuracy out to many decimal places.  Furthermore, some areas of
  the regularization path may contain models that are nonidentifiable or
  nearly singular.  Thus, the algorithm may fail to satisfy
  convergence criteria at certain points while yielding accurate
  solutions over the region of interest.  The default behavior warning
  the user when convergence criteria are not met is often distracting,
  and can be modified with \code{warn.conv} (convergence can always be
  checked later by inspecting the value of \code{iter}).  If models are
  not converging, consider increasing \code{eps}, increasing
  \code{n.lambda}, or increasing \code{lambda.min} before increasing
  \code{max.iter}.}
\value{
  An object with S3 class \code{"grpreg"} containing:
  \item{beta}{The fitted matrix of coefficients.  The number of rows is
    equal to the number of coefficients, and the number of columns is
    equal to \code{n.lambda}.}
  \item{beta.std}{Same as \code{beta}, only on the standardized scale.}
  \item{lambda}{The sequence of \code{lambda} values in the path.}
  \item{penalty}{The user-supplied value of \code{penalty}.}
  \item{df}{A vector of length \code{n.lambda} containing estimates of
    effective number of model parameters all the points along the
    regularization path.  For details on how this is calculated, see
    reference.}
  \item{iter}{A vector of length \code{n.lambda} containing the number
    of iterations until convergence at each value of \code{lambda}.}
  \item{par}{Processed list of parameters used by \code{grpreg}
    during fitting.}
  \item{Data}{Processed list containing data used by \code{grpreg}
    during fitting.}}
\references{Breheny, P. and Huang, J. (2008) \emph{Penalized Methods for
    Bi-level variable selection}.  Tech report No. 393, Department of
  Statistics and Actuarial Science, University of Iowa.\url{http://www.stat.uiowa.edu/techrep/tr393.pdf}}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{plot} and \code{select} methods.}
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
## An example of an unimportant failure to converge.
## Note that the plot looks fine, and that the algorithm fails to
## converge for only one value of lambda.
fit2.gBridge <- grpreg(Data.binomial,"gBridge",lambda.max=0.06)
fit2.gBridge$iter
plot(fit2.gBridge)
fit2.gMCP <- grpreg(Data.binomial,"gMCP")

select(fit2.gLasso)
}
\keyword{models}
\keyword{regression}
