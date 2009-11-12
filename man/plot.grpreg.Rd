\name{plot.grpreg}
\alias{plot.grpreg}
\title{plot coefficients from a "grpreg" object}
\description{Produces a plot of the coefficient paths for a fitted \code{"grpreg"} object.}
\usage{
\method{plot}{grpreg}(x, color=TRUE, alpha=1, type="l", pch=1:p, lty=1,
xlab=expression(lambda), ylab=expression(hat(beta)), legend.loc, ...) 
}
\arguments{
  \item{x}{Fitted \code{"grpreg"} model.}
  \item{color}{TRUE for color plots, FALSE for black and white.
    Defaults to TRUE; overridden by \code{col} if both are specified.}
  \item{alpha}{Controls alpha-blending.  Default is alpha=1.}
  \item{type}{Type of plot.  Default is a line plot.  See \code{matplot}
    for details.}
  \item{lty}{Line type.  Defaults to solid.  See See \code{matplot}
    for details.}
  \item{pch}{Plotting character.  Defaults to no plotted characters.
    See \code{matplot} for details.}
  \item{xlab}{Horizontal axis label.}
  \item{ylab}{Vertical axis label.}
  \item{legend.loc}{Where should the legend go?  If left unspecified, no
    legend is drawn.  See \code{legend} for details.}
  \item{\dots}{Other graphical parameters to \code{matplot}}
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
fit <- grpreg(X,y,group,penalty="gMCP",lambda.min=.5)
plot(fit)
plot(fit,color=FALSE)
plot(fit,color=TRUE,legend.loc="topleft")
}
\keyword{models}
\keyword{regression}
