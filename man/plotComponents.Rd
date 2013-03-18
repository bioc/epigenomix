\name{plotComponents}
\alias{plotComponents}
\alias{plotComponents,MixModel-method}

\title{
  Plots the mixture density together with the densities of all single components.
}

\description{
  This function plots the mixture pdf, the estimated data pdf and
  the weighted pdfs of all components of the given mixture model.
  The plot is useful to assess the fit of the model.
}

\usage{
plotComponents(object, density = FALSE, ...)
}

\arguments{
  \item{object}{
    A \code{\link{MixModel-class}} object to be plotted.
  }
  \item{density}{
    A logical indicating whether the data distribution should be plotted
    as histogram (FALSE) or as density (TRUE) estimated using kernel
    density estimation.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{plot}}.
  }
}

\details{
  If the argument "col" is given, the first color is used for the mixture
  pdf. The following colors (2 to n+1) are used for the n mixture components'
  pdfs. If \code{density} is set to \code{TRUE}, a further color (n+2) must be
  given, that is used for the data pdf.
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\seealso{
  \code{\link{MixModel-class}}
}

\examples{
z <- c(rnorm(100, 0, 1), rnorm(100, 0, 2), rexp(100, 1/2), -rexp(100, 1/2))
mm <- mlMixModel(z, normNull=1:2, expNeg=3, expPos=4,
   sdNormNullInit=c(1, 2), rateExpNegInit=c(1/2), rateExpPosInit=c(1/2),
   pi=c(1/4, 1/4, 1/4, 1/4), maxIter=50, tol=0.01)
plotComponents(mm)
}
