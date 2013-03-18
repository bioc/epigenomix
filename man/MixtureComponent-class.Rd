\name{MixtureComponent-class}
\Rdversion{1.1}
\docType{class}
\alias{MixtureComponent-class}
\alias{show,MixtureComponent-method}

\title{Class \code{"MixtureComponent"}}

\description{
  A class representing a mixture component.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MixtureComponent", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} giving the
      name or type of the mixture component.}
    \item{\code{parameters}:}{Object of class \code{"list"} storing the
      parameters of corresponding distribution.}
    \item{\code{pdf}:}{Object of class \code{"function"} giving the pdf
      of the mixture component.}
    \item{\code{color}:}{Object of class \code{"character"} giving the
      color of the component that is used by plotting methods.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "MixtureComponent")}: A method
      plotting a summary of the component on screen.}
  }
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\note{
  The element in \code{parameters} should be named by the argument names
  of \code{pdf} such that this call works:
  \code{do.call(object@pdf, c(list(x=data), object@parameters))}
}

\seealso{
  \code{\linkS4class{MixModel}}
}

\examples{
showClass("MixtureComponent")
}

\keyword{classes}
