\name{select.grpreg}
\alias{select}
\alias{select.grpreg}
\title{select an point along a grpreg path}
\description{Selects a point along the regularization path of a fitted
  grpreg object according to the AIC, BIC, or GCV criteria.}
\usage{
select(obj,\dots)
\method{select}{grpreg}(obj, X, y, criterion="BIC", df.method="default",
...)
}
\arguments{
  \item{obj}{A fitted grpreg object.}
  \item{X,y}{The data used to fit the model.}
  \item{criterion}{The criterion by which to select the regularization
    parameter.  One of \code{"AIC"}, \code{"BIC"}, or \code{"GCV"};
    default is \code{"BIC"}.} 
  \item{df.method}{How should effective model parameters be calculated?
    One of: \code{"active"}, which counts the number of nonzero
    coefficients; or \code{"default"}, which uses the calculated
    \code{df} returned by \code{grpreg}.  Default is \code{"default"}.}
  \item{\dots}{For S3 method compatibility.}}
\details{
  The criteria are defined as follows, where \code{L} is the loss function
  (usually, \code{-loglik}) and \code{n} is the sample size: \cr
  \deqn{AIC = 2*L + 2*df} \cr
  \deqn{BIC = 2*L + log(n)*df} \cr
  \deqn{GCV= 2*L/((1-df/n)^2)}}
\value{
  A list containing:
  \item{lambda}{The selected value of the regularization parameter,
    \code{lambda}.}
  \item{beta}{The vector of coefficients at the chosen value of
    \code{lambda}.}
  \item{df}{The effective number of model parameters at the chosen value
    of \code{lambda}.}
  \item{IC}{A vector of the calculated model selection criteria for each
    point on the regularization path.}
  }
\references{Breheny, P. and Huang, J. (2009) Penalized methods for
  bi-level variable selection.  \emph{Statistics and its interface},
  \bold{2}: 369--380.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{grpreg}}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="gMCP")
select(fit,X,y)
select(fit,X,y,crit="AIC",df="active")
plot(fit$lambda,select(fit,X,y)$IC,xlim=rev(range(fit$lambda)))
}
\keyword{models}
\keyword{regression}
