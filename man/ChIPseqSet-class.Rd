\name{ChIPseqSet-class}
\Rdversion{1.1}
\docType{class}
\alias{ChIPseqSet-class}
\alias{ChIPseqSet}
\alias{ChIPseqSet,matrix,GRanges-method}
\alias{ChIPseqSet,matrix,GRangesList-method}
\alias{chipVals<-}
\alias{chipVals<-,ChIPseqSet,matrix-method}
\alias{chipVals}
\alias{chipVals,ChIPseqSet-method}

\title{Class \code{"ChIPseqSet"}}
\description{
  A class for storing count data obtained from ChIP-seq experiments by
  counting the number of reads lying within regions. The class extends
  \code{\linkS4class{SummarizedExperiment}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{ChIPseqSet(chipVals=countDataMatrix,
                   rowData=genomicRegions,
                   colData=DataFrame(row.names=colnames(countDataMatrix)),
                   ...)}.
  However, one will most likely create a \code{ChIPseqSet} object by
  calling \code{\link{summarizeReads}}.
}

\section{Slots}{
  \describe{
    \item{\code{exptData}:}{An optional \code{SimpleList} of arbitrary
      content describing the overall experiment.}
    \item{\code{rowData}:}{Object of class \code{"GRanges"} or \code{"GRangesList"}
      containing the genomic regions where the reads were counted.}
    \item{\code{colData}:}{Object of class \code{"DataFrame"}
      containing information on variable values of the samples. Some methods
      require the total library size of each sample stored in a column
      titled totalCounts.}
    \item{\code{assays}:}{Object of class \code{SimpleList} of a
      matrix, named \code{chipVals} containing the read counts 
      per genomic region.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{SummarizedExperiment}"}, directly.
}
\section{Methods}{
  \describe{
    \item{chipVals}{\code{signature(object = "ChIPseqSet")}: Returns the
      matrix with read counts.}
    \item{chipVals<-}{\code{signature(object = "ChIPseqSet", value = "matrix")}:
      Sets the matrix with read counts.}
  }
}

\author{
Hans-Ulrich Klein (hklein@broadinstitute.org)
}

\seealso{
\code{\link{summarizeReads}}, \code{\link{normalizeChIP}}
}

\examples{
showClass("ChIPseqSet")
}

\keyword{classes}

