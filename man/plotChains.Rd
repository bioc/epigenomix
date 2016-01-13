\name{plotChains}
\alias{plotChains}
\alias{plotChains,MixModelBayes-method}

\title{
  Produces trace plots for a Bayesian mixture model
}

\description{
  This function method draws trace plots for a Bayesian mixture model,
  e.g. visualizes the course of the Markov Chains. Inspection of the
  Markov Chains is important to determine convergence of the chains,
  which is necessary for sensible results.
}

\usage{
plotChains(object, chain, component, itb = 1, thin = 1, cols, ...)
}

\arguments{
  \item{object}{
    An object of \code{\link{MixModelBayes-class}}
  }
  \item{chain}{
    A character of length one giving the name of the paramter, which
    chain should be plotted. Can be omitted, if \code{component} is
    given. Then, all parameters of the given components are plotted.    
  }
  \item{component}{
    An integer specifying the components, which parameter chains should
    be plotted. Can be omitted, if \code{chain} is given. Then, all
    trace plots are generated for all components having the parameter
    specified via argument \code{chain}.
  }
  \item{itb}{
    Number of iterations used for burn-in. The burn-in is relative to
    the output of \code{\link{bayesMixModel}}, e.g., any burn-in
    specified here is added to the burn-in that was specified when
    calling \code{\link{bayesMixModel}}.
  }
  \item{thin}{
    Thinning value for the iterations after burn-in. The thinning is
    relative to the output of \code{\link{bayesMixModel}}, e.g., any
    thinning specified here multiplies by the thinning that was
    specified in \code{\link{bayesMixModel}}.    
  }
  \item{cols}{
    Number of columns to be used in the plot. Optional, if omitted, the
    number of columns and rows are choosen be the method itself.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{plot}}.
  }
}

\details{
  The number of iterations necessary until a Markov chain reaches
  stationarity depends on the specific model and data. For any inference
  based on Markov Chain Monte Carlo methods, it is therefore necessary
  to inspect the convergence of Markov Chains. One way to do this is
  visual inspection of trace plots using this method.

  If argument \code{main} is passed to this method, it should have as
  many elements as chains are plotted. Otherwise, vector \code{main} is
  reapted.
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
  Martin Schaefer (martin.schaefer@udo.edu)
}

\seealso{
\code{\link{bayesMixModel}}, \code{\link{MixModelBayes-class}}
}

\examples{
z <- c(rnorm(1000, 0, 3), rnorm(1000, 0, 5), rexp(1000, 5), -rexp(1000, 5))
mm <- bayesMixModel(z, normNull=1:2, expNeg=3, expPos=4,
  sdNormNullInit=c(1, 2), rateExpNegInit=8, rateExpPosInit=8,
  shapeNorm0=c(1, 1), scaleNorm0=c(1, 1),
  shapeExpNeg0=c(1, 1), scaleExpNeg0=c(1, 1),
  shapeExpPos0=c(1, 1), scaleExpPos0=c(1, 1),
  sdAlpha=1, itb=200, nmc=1000, thin=10)
plotChains(mm, chain="pi")
plotChains(mm, component=c(2,3))
}
