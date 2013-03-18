\name{MixModelML-class}
\Rdversion{1.1}
\docType{class}
\alias{MixModelML-class}
\alias{convergence}
\alias{convergence,MixModelML-method}

\title{Class \code{"MixModelML"}}

\description{
  This class stores a mixture model fitted by a maximum likelihood approach.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MixModelML", ...)}.
Usually, objects are created by \code{\link{mlMixModel}}.
}

\section{Slots}{
  \describe{
    \item{\code{convergence}:}{Object of class \code{"list"} storing
      information about the convergence of the EM algorithm.}
    \item{\code{mmData}:}{Object of class \code{"numeric"} storing the data.}
    \item{\code{configuration}:}{Object of class \code{"list"} storing 
        configuration. See notes for details.}
    \item{\code{results}:}{Object of class \code{"list"} storing
      results. See notes for details.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{MixModel}"}, directly.
}

\section{Methods}{
  \describe{
    \item{convergence}{\code{signature(object = "MixModelML")}: Access
      to the convergence information.}
  }
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\note{
  In addition to the content described in \code{\linkS4class{MixModel}},
  the following elements are present:
  Slot \code{configuration}:
  \enumerate{
    \item{convergence}{A list storing the maximum number of allowed
      iterations. And delta log likelihood limit, that is interpreted as
      convergence, if the delta log likelihood falls below that limit.}
  }
  Slot \code{results} is exactly like in \code{\linkS4class{MixModel}}.
  Slot \code{convergence}:
  \enumerate{
    \item{iterations}{Number of iterations ran.}
    \item{deltaLogLik}{Delta of log likelihood observed in the last
      iteration.}
    \item{logLik}{Log likelihood of the model fit.}
  }
}

\seealso{
\code{\link{mlMixModel}}
\code{\linkS4class{MixModel}}
}

\examples{
showClass("MixModelML")
}

\keyword{classes}
