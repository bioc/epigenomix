\name{MixModel-class}
\Rdversion{1.1}
\docType{class}
\alias{MixModel-class}
\alias{classification}
\alias{classification,MixModel,character-method}
\alias{classification,MixModel,missing-method}
\alias{components}
\alias{components,MixModel-method}
\alias{mmData}
\alias{mmData,MixModel-method}
\alias{dim,MixModel-method}
\alias{length,MixModel-method}
\alias{listClassificationMethods}
\alias{listClassificationMethods,MixModel-method}
\alias{show,MixModel-method}
\alias{summary,MixModel-method}
\alias{weights,MixModel-method}

\title{Class \code{"MixModel"}}

\description{
  This class stores a fitted mixture model.
}

\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{mmData}:}{Object of class \code{"numeric"} storing the data.}
    \item{\code{configuration}:}{Object of class \code{"list"} storing
      configuration. See notes for details.}
    \item{\code{results}:}{Object of class \code{"list"} storing
      results. See notes for details.}
  }
}

\section{Methods}{
  \describe{
    \item{classification}{\code{signature(object = "MixModel", method =
	"character")}: Assess classification results.}
    \item{classification}{\code{signature(object = "MixModel", method =
	"missing")}: Assess classification results.}
    \item{components}{\code{signature(object = "MixModel")}: Assess
      mixture components.}
    \item{mmData}{\code{signature(object = "MixModel")}: Assess data.}
    \item{dim}{\code{signature(x = "MixModel")}: Assess dimension,
      i.e. numer of data points and number of components.}
    \item{length}{\code{signature(x = "MixModel")}: Number of data points.}
    \item{listClassificationMethods}{\code{signature(object =
	"MixModel")}: List available classification methods.}
    \item{show}{\code{signature(object = "MixModel")}: Print an object
      of \code{MixModel} on screen. }
    \item{summary}{\code{signature(object = "MixModel")}: Returns a list
      of data frames summarizing the parameter estimations for each component.}
    \item{weights}{\code{signature(object = "MixModel")}: Asses the
      components weights.}
  }
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\note{
  Slots \code{configuration} and \code{results} are lists with named elements. The following
  elements make up the minimum set of element that must be present. Depending on the method
  that was used to fit the mixture model, more elements may be present.

  Slot \code{configuration} has at least one element.
  \enumerate{
    \item{inits}{A list with at least two elements: \code{component} and \code{pi}. 
                 \code{components} contains a list of objects of
                 \code{\link{MixtureComponent-class}} storing the inital parameters
                 of the mixture components. \code{pi} is a vector storing
                 the initial components' weights.}
  }
  Slot \code{results} has at least three elements.
  \enumerate{
    \item{components}{A list of objects of \code{\link{MixtureComponent-class}} storing
                      the fitted mixture components.}
    \item{pi}{A numeric vector holding the estimated components' weights.}
    \item{classification}{A list of numeric vectors of the same length as \code{data}
                          storing the classification results.}
  }
}

\seealso{
  \code{\link{mlMixModel}}
  \code{\link{bayesMixModel}}
  \linkS4class{MixModelML}
  \linkS4class{MixModelBayes}
}

\examples{
showClass("MixModel")
}

\keyword{classes}
