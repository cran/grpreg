\name{grpreg-internal}
\title{Internal grpreg functions}
\alias{calcL}
\alias{formatData}
\alias{criteria.grpreg}
\alias{formatPar}
\alias{unstandardize}
\description{Internal grpreg functions}
\usage{
calcL(Data,beta)
formatData(Data)
formatPar(Data,par)
criteria.grpreg(fit,method="BIC",df)
unstandardize(beta,Data)
}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\details{These are not intended for use by users. \code{calcL}
  calculates the loss function of a model given a formatted Data object
  and a vector of coefficients beta.  \code{formatData} formats the
  user-supplied data in the form required internally; \code{formatPar}
  does the same with control parameters.  \code{criteria.grpreg}
  calculates the AIC, BIC, and GCV information criteria.
  \code{unstandardize} returns standardized coefficients to their
  original scale.
  }
\keyword{internal}
