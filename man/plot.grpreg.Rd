\name{plot.grpreg}
\alias{plot.grpreg}
\title{plot coefficients from a "grpreg" object}
\description{Produces a plot of the coefficient paths for a fitted \code{"grpreg"} object.}
\usage{
\method{plot}{grpreg}(x,col=NULL,alpha=FALSE,color=TRUE,type="l",lty=1,pch=NULL,xlab=expression(lambda),ylab=expression(hat(beta)),...)
}
\arguments{
  \item{x}{Fitted \code{"grpreg"} model.}
  \item{col}{A vector of colors with length equal to the number of
    penalized coefficients in the model.  By default, coefficients will
    be colored by group, with colors equally spaced in the HSV color space.}
  \item{alpha}{Controls alpha-blending.  Can be: FALSE, which turns off
    alpha-blending; TRUE, causing the function to automatically generate
    a reasonable level of alpha-blending depending on the number of
    coefficients; or a numeric vector of user-supplied alpha values.
    The usual recycling rules apply.}
  \item{color}{TRUE for color plots, FALSE for black and white.
    Defaults to TRUE; overrides \code{col} if both are specified.}
  \item{type}{Type of plot.  Default is a line plot.  See \code{matplot}
    for details.}
  \item{lty}{Line type.  Defaults to solid.  See See \code{matplot}
    for details.}
  \item{pch}{Plotting character.  Defaults to no plotted characters.See
    \code{matplot} for details.}
  \item{xlab}{Horizontal axis label.}
  \item{ylab}{Vertical axis label.}
  \item{\dots}{Other graphical parameters to \code{matplot}}
}
\references{Breheny, P. and Huang, J. (2008) \emph{Penalized Methods for
    Bi-level variable selection}.  Tech report No. 393, Department of
  Statistics and Actuarial Science, University of Iowa.\url{http://www.stat.uiowa.edu/techrep/tr393.pdf}}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{grpreg}}
\examples{
data(birthwt.grpreg)
Data.gaussian <- list(y=birthwt.grpreg$bwt,
                      X=as.matrix(birthwt.grpreg[,c(-1,-2)]),
                      family="gaussian",
                      group=c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8))

fit <- grpreg(Data.gaussian,"gMCP")

plot(fit)
plot(fit,color=FALSE,lty=rep(1:8,table(Data.gaussian$group)))
plot(fit,col=rep(c("red","blue","black","green","yellow","brown","purple","gray"),table(Data.gaussian$group)))
}
\keyword{models}
\keyword{regression}
