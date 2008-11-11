\name{select.grpreg}
\alias{select}
\alias{select.grpreg}
\title{select an point along a grpreg path}
\description{Selects a point along the regularization path of a fitted
  grpreg object according to the AIC, BIC, or GCV criteria.}
\usage{
select(obj,\dots)
\method{select}{grpreg}(obj,method="BIC",df="default",\dots)
}
\arguments{
  \item{obj}{A fitted grpreg object.}
  \item{method}{The criteria by which to select the regularization
    parameter.  One of \code{"AIC"}, \code{"BIC"}, or \code{"GCV"}; default is \code{"BIC"}.}
  \item{df}{How should effective model parameters be calculated?  One
    of: \code{"active"}, which counts the number of nonzero
    coefficients; or \code{"default"}, which uses the calculated
    \code{df} returned by \code{grpreg}.  Default is
    \code{"default"}.}
  \item{\dots}{For S3 method compatibility.}}
\details{
  The criteria are defined as follows, where \deqn{L} is the loss function
  (usually, \code{-loglik}) and \deqn{n} is the sample size: \cr
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
\references{Breheny, P. and Huang, J. (2008) \emph{Penalized Methods for
    Bi-level variable selection}.  Tech report No. 393, Department of
  Statistics and Actuarial Science, University of Iowa.\url{http://www.stat.uiowa.edu/techrep/tr393.pdf}}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{grpreg}}
\examples{
data(birthwt.grpreg)
Data <- list(y=birthwt.grpreg$low,
             X=as.matrix(birthwt.grpreg[,c(-1,-2)]),
             family="binomial",
             group=c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8))

fit <- grpreg(Data,"gMCP")

select(fit)
select(fit,method="AIC",df="active")
plot(fit$lambda,select(fit)$IC,type="l")
}
\keyword{models}
\keyword{regression}

 
