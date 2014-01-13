\name{plotClassification}
\alias{plotClassification}
\alias{plotClassification,MixModel-method}

\title{
  Plot classification obtained from a mixture model.
}

\description{
  This method visualizes the assignment of data points to the mixture
  components of the given mixture model. The components are plotted on
  the y-axis and the data on the x-axis. Data points are plotted in the
  color of the respective mixture component.
}

\usage{
plotClassification(object, method, ...)
}

\arguments{
  \item{object}{
    An object of \code{\link{MixModel-class}}.
  }
  \item{method}{
    Depending on the type of the mixture model (ML, Bayes),
    different approaches to obtain a classification are available. Also
    the default approach may vary.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{plot}}.
  }
}

\details{
  If \code{method} is given, it must be a valid option for method
  \code{classification}. E.g., if \code{\link{bayesMixModel}} was
  used to create the mixture model, valid options are "maxDens",
  "median" and "mode".

  Arguments "col" and "pch" can be given to specify the color and
  the shape of the points plotted. Their length must equal to the
  number of components.
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\seealso{
  \code{\link{MixModel-class}}
  \code{\link{listClassificationMethods}}
}

\examples{
z <- c(rnorm(100, 0, 10), rnorm(100, 0, 2), rexp(100, 1/2), -rexp(100, 1/2))
mm <- mlMixModel(z, normNull=1:2, expNeg=3, expPos=4,
   sdNormNullInit=c(1, 2), rateExpNegInit=c(1/2), rateExpPosInit=c(1/2),
   pi=c(1/4, 1/4, 1/4, 1/4), maxIter=50, tol=0.01)
plotClassification(mm)
}
