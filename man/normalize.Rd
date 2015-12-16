\name{normalize}
\alias{normalize}
\alias{normalize,ChIPseqSet-method}
\alias{normalize,ExpressionSet-method}

\title{Normalization of ChIP-seq and other count data}
\description{This function implements some methods for between-sample normalization
  of count data. Although these methods were developed for RNA-seq data,
  they are also useful for ChIP-seq data normalization after reads
  were counted within regions or bins. Some methods may also be applied
  to count data after within-sample normalization (e.g. TPM or RPKM
  values).}
\usage{
\S4method{normalize}{ChIPseqSet}(object, method, isLogScale = FALSE, trim = 0.3, totalCounts)
\S4method{normalize}{ExpressionSet}(object, method, isLogScale = FALSE, trim = 0.3, totalCounts)
}

\arguments{
  \item{object}{
    An object of class \code{\linkS4class{ChIPseqSet}} or
    \code{\linkS4class{ExpressionSet}} that contains the raw data.
  }
  \item{method}{
    Normalization method, either "scale", "scaleMedianRegion",
    "quantile" or "tmm".
  }
  \item{isLogScale}{
    Indicates whether the raw data in \code{object} is already
    logarithmized. Default value is \code{FALSE}. Logarithmized
    data will be returned on the log scale, non logarithmized data
    will remain on its original scale.
  }
  \item{trim}{
    Only used if \code{method} is "tmm". Indicates the fraction of
    data points that should be trimmed before calculating the mean.
    Default value is 0.3.
  }
  \item{totalCounts}{
    Only used if \code{method} is "scale". A vector giving the total
    number of reads for each sample. The Vector's length must equal the
    number of samples in \code{object}. Default values are
    the sums over all features for each sample
    (i.e. colsums of \code{object}).
  }
}

\details{
 The following normalization methods are implemented:
  \enumerate{
    \item{scale}{Samples are scaled by a factor such that all samples
      have the same number \eqn{N} of reads after normalization, where 
      \eqn{N} is the median number of reads observed accross all samples. If
      the argument \code{totalCounts} is missing, the total numbers of
      reads are calculated from the given data. Otherwise, the values
      in \code{totalCounts} are used.}
    \item{scaleMedianRegion}{The scaling factor \eqn{s_j} for the
      \eqn{j}-th sample is defined as
      \deqn{s_j = median_i \frac{k_{ij}}{\prod_{v=1}^m k_{iv}}.}
      \eqn{k_{ij}} is the value of region \eqn{i} in sample
      \eqn{j}. See Anders and Huber (2010) for details.}
    \item{quantile}{Quantile normalization is applied to the ChIP-seq
      values such that each sample has the same cdf after
      normalization.}
    \item{tmm}{The trimmed mean M-value (tmm) normalization was
      proposed by Robinson and Oshlack (2010). Here, the logarithm
      of the scaling factor for sample \eqn{i} is calculated as
      the trimmed mean of
      \deqn{\log(k_{i,j}/m_{j}).}
      Variable \eqn{m_{j}} denotes the geometric mean of region \eqn{j}.
      Argument \code{trim} is set to 0.3 as default value, so that
      the smallest 15\% and the largest 15\% of the log ratios are
      trimmed before calculating the mean.}
  }
}

\value{
  An object of the same class as the input \code{object} with
  the normalized data.
}

\references{
  Anders and Huber. Differential expression analysis for sequence count
  data. Genome Biol. 2010;11(10):R106.\\
  Robinson and Oshlack. A scaling normalization method for differential
  expression analysis of RNA-seq data. Genome Biol. 2010;11(3):R25
}

\author{
  Hans-Ulrich Klein (hklein@broadinstitute.org)
}

\examples{
  set.seed(1234)
  chip <- matrix(c(rpois(20, lambda=10), rpois(20, lambda=20)), nrow=20,
                 dimnames=list(paste("feature", 1:20, sep=""), c("sample1", "sample2")))
  rowRanges <- GRanges(IRanges(start=1:20, end=1:20),
                     seqnames=c(rep("1", 20)))
  names(rowRanges) = rownames(chip)
  cSet <- ChIPseqSet(chipVals=chip, rowRanges=rowRanges)

  tmmSet <- normalize(cSet, method="tmm", trim=0.3)
  mean(log(chipVals(tmmSet))[, 1], trim=0.3) -
      mean(log(chipVals(tmmSet))[, 2], trim=0.3) < 0.01

  quantSet <- normalize(cSet, method="quantile")
  all(quantile(chipVals(quantSet)[, 1]) == quantile(chipVals(quantSet)[, 2]))
}

\keyword{normalization}
