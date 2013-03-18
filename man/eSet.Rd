\name{eSet}
\alias{eSet}
\docType{data}
\title{
  Example gene expression data set.
}
\description{
  The \code{\linkS4class{ExpressionSet}} stores 2 replicates for each of
  two different conditions. Data were obtained from Affymetrix
  MouseGene 1.0 ST arrays.
}

\usage{data(eSet)}

\format{
  An object of class \code{\linkS4class{ExpressionSet}}.
}

\details{
  The example data contains a subset of 200 probesets located on
  chromosome 1. Data were RMA normalized.
}

%%\source{
%%
%%}

\examples{
data(eSet)
eSet
pData(eSet)
}

\keyword{datasets}
