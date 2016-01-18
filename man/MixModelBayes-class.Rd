\name{MixModelBayes-class}
\Rdversion{1.1}
\docType{class}
\alias{MixModelBayes-class}
\alias{chains}
\alias{chains,MixModelBayes-method}
\alias{acceptanceRate}
\alias{acceptanceRate,MixModelBayes-method}

\title{Class \code{"MixModelBayes"}}
\description{
  This class stores a Bayesian mixture model fitted by MCMC methods.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MixModelBayes", ...)}. 
}
\section{Slots}{
  \describe{
    \item{\code{chains}:}{Object of class \code{"list"} storing the course of the Markov chains for each parameter.}
    \item{\code{mmData}:}{Object of class \code{"numeric"} storing the data.}
    \item{\code{configuration}:}{Object of class \code{"list"} storing 
        configuration. See notes for details.}
    \item{\code{results}:}{Object of class \code{"list"} storing results. See notes for details.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{MixModel}"}, directly.
}

\section{Methods}{
  \describe{
    \item{chains}{\code{signature(object = "MixModelBayes")}: Gives
      access to the \code{chains} slot of the object.}
     \item{acceptanceRate}{\code{signature(object = "MixModelBayes")}: Gives
       the acceptance rate for the parameter of the Dirichlet distribution.
       Acceptance rates between 0.3 and 0.7 are usually desired. Values
       not smaller than 0.1 (not larger than 0.9) might still be
       acceptable. The acceptance rate is only meaningful if the option
       \code{weightsPrior} was set to the Finite-dimensional Dirichlet
       prior.
     }
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
    \item{inits}{As in \code{\linkS4class{MixModel}}.}
    \item{priors}{A list specifying the prior distributions for the
      parameters of the components and the parameter of the Dirichlet
      process.}
    \item{chain}{A list with the technical specifications for the Markov
      Chains.}
  }
  Slot \code{results} is exactly like in \code{\linkS4class{MixModel}}.
  Slot \code{chains}:
  \enumerate{
    \item{components}{A list giving the values for the parameters of the
      components in each iteration after burn-in and application of thinning.}
    \item{pi}{A matrix giving the values for the weights pi of the
      components in each iteration after burn-in and application of thinning.}
    \item{dirichletParameter}{A vector giving the values for dirichlet 
      Parameter in each iteration after burn-in and application of thinning.}
    \item{classification}{A matrix giving the number of genes classified to
      each components in each iteration after burn-in and application of thinning.}
  }
}


\seealso{
  \code{\link{bayesMixModel}}
\code{\linkS4class{MixModel}}
}

\examples{
showClass("MixModelBayes")
}

\keyword{classes}
