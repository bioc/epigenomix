\name{mlMixModel}
\alias{mlMixModel}
\alias{mlMixModel,numeric-method}

\title{
  Fits a mixture model using the maximum likelihood principle.
}

\description{
  This method calculates the maximum likelihood estimations of a mixture
  model using the expectation-maximization (EM) algorithm. The mixture
  model may consists of normal components (with a fixed expectation of 0)
  and exponential components, which may be mirrored in order to model
  negative values.
}

\usage{
mlMixModel(z, normNull = c(), expNeg = c(), expPos = c(), sdNormNullInit = c(), rateExpNegInit = c(), rateExpPosInit = c(), piInit = c(), maxIter = 500, tol = 0.001)
}

\arguments{
  \item{z}{
    Observed values.
  }
  \item{normNull}{
    Indices of the normal components (that have \eqn{\mu=0}{mu=0}).
  }
  \item{expNeg}{
    Indices of the mirrored exponential components.
  }
  \item{expPos}{
    Indices of the exponential components.
  }
  \item{sdNormNullInit}{
    Initial standard deviations of the normal components.
  }
  \item{rateExpNegInit}{
    Initial rates ("lambda") of the exponential components.
  }
  \item{rateExpPosInit}{
    Initial rates ("lambda") of the exponential components.
  }
  \item{piInit}{
    Initial weights of the components.
  }
  \item{maxIter}{
    Maximum number of iterations.
  }
  \item{tol}{
    Threshold for convergence. The minimum log likelihood gain between two
    iterations that must be achieved to continue.
  }
}

\details{
  The EM algorithm is known to converge slowly in some cases and local
  maxima may avoid finding the optimal solution. Users should try
  different initial values and different convergence criteria.

  The components' indices do not influence the result, but may influence
  the order in which components are listed or ploted by downstream
  methods. Indices must be successive integers from 1 to n.
}

\value{
  An object of \code{\link{MixModelML-class}} storing results, data,
  initial values and information about the convergence.
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\seealso{
  \code{\link{MixModelML-class}}
}

\examples{
z <- c(rnorm(1000, 0, 0.5), rnorm(1000, 0, 1))
mm <- mlMixModel(z, normNull=1:2, sdNormNullInit=c(0.1, 0.2),
   pi=c(1/2, 1/2), maxIter=500, tol=0.001)
mm

z <- c(rnorm(1000, 0, 3), rnorm(1000, 0, 5), rexp(1000, 5), -rexp(1000, 5))
mm <- mlMixModel(z, normNull=1:2, expNeg=3, expPos=4,
   sdNormNullInit=c(1, 2), rateExpNegInit=8, rateExpPosInit=8,
   pi=c(1/4, 1/4, 1/4, 1/4), maxIter=500, tol=0.001)
mm
}
