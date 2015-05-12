\name{integrateData}
\alias{integrateData}
\alias{integrateData,ExpressionSet,ChIPseqSet,character,character-method}
\alias{integrateData,ExpressionSet,ChIPseqSet,character,missing-method}
\alias{integrateData,ExpressionSetIllumina,ChIPseqSet,character,character-method}
\alias{integrateData,ExpressionSetIllumina,ChIPseqSet,character,missing-method}


\title{
Calculates a normalized correlation score from ChIP-seq and microarray
gene expression data.
}
\description{
This function calculates the product of the standardized differences
between two conditions in ChIP-seq data and the respective standardized
differences in gene expression data. A score close to zero means that there
are no (large) differences in at least one of the two data sets. If the
score is positive, equally directed differences exist in both data
sets. In case of a negative score, differences have unequal signs in the
two data sets.
}

\usage{
integrateData(expr, chipseq, factor, reference)
}

\arguments{
\item{expr}{
  An \code{\link{ExpressionSet}} holding the gene expression data.
}
\item{chipseq}{
  A \code{\linkS4class{ChIPseqSet}} holding the ChIP-seq data.
}
\item{factor}{
  A \code{character} giving the name of the factor that describes the
  conditions to be compared. The factor must be present in the pheno data
  slot of the objects expr and chipseq. Further, the factor must have
  exactly two levels and the level names must be the same in both objects.
}
\item{reference}{
  Optionally, the name of the factor level that should be used as
  reference. If missing, the first level of \code{factor} in the object
  \code{expr} is used.
}
}

\details{
  Let A and B denote the gene expression value of one probe in the group
  of interest and in the reference group defined by the argument
  \code{reference}. And let X and Y be the ChIP-seq values assigned to
  that probe. This functions returnes for each probe 
  \deqn{Z = (A-B)/\sigma_{ge} \times (X-Y)/\sigma_{chip},}
  where \eqn{\sigma_{ge}} is the standard deviation estimated from all
  observed difference in the gene expression data and \eqn{\sigma_{chip}}
  the standard deviation in the ChIP-seq data.

  If there is more than one sample in any group and data set, the
  average of the replicates is calcuated first and than plugged into
  the formula above.

  Not all features in \code{expr} must also be in \code{chipseq} and
  vice versa. Features present in only one of the two data types are
  omitted.
}

\value{
  A matrix with five columns. The first 4 columns store the (average)
  expression values and the (average) ChIP-seq values for each of the two
  conditions. The fith columns store the correlation score. The row names
  equal common feature names of \code{expr} and \code{chipseq}.
}

\author{
Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\seealso{
\code{\link{summarizeReads}}
\code{\link{normalizeChIP}}
}

\examples{
ge <- matrix(c(5,12,5,11,11,10,12,11), nrow=2)
row.names(ge) <- c("100_at", "200_at")
colnames(ge) <- c("c1", "c2", "t1", "t2")
geDf <- data.frame(status=c("control", "control", "treated", "treated"),
  row.names=colnames(ge))
eSet <- ExpressionSet(ge, phenoData=AnnotatedDataFrame(geDf))

chip <- matrix(c(10,20,20,22), nrow=2)
row.names(chip) <- c("100_at", "200_at")
colnames(chip) <- c("c", "t")
rowRanges <- GRanges(IRanges(start=c(10,50), end=c(20,60)), seqnames=c("1","1"))
names(rowRanges) = c("100_at", "200_at")
chipDf <- DataFrame(status=factor(c("control", "treated")),
  totalCount=c(100, 100),
  row.names=colnames(chip))
cSet <- ChIPseqSet(chipVals=chip, rowRanges=rowRanges, colData=chipDf)

integrateData(eSet, cSet, factor="status", reference="control")
}
