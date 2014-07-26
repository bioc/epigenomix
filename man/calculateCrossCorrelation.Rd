\name{calculateCrossCorrelation}
\alias{calculateCrossCorrelation}
\alias{calculateCrossCorrelation,GRanges-method}


\title{
Calculate the cross correlation for a given GRanges object.
}
\description{
This method calculates the cross correlation, i.e. the pearson
correlation between the coverages of the positive and negative strand
from a DNA sequencing experiment. The cross correlation can be used as a
quality measure in ChIP-seq experiments (Kharchenko et al. 2008). Cross
correlation can also be used to estimate the fragment size by deteriming
the shift (given in base pairs) that maximizes the cross correlation. 
}

\usage{
calculateCrossCorrelation(object, shift=c(200,250,300), bin=10, mode="none", minReads=10000, chrs=NA, mc.cores=1)
}

\arguments{
  \item{object}{
  An \code{\link{GRanges}} object containing the aligned reads.
}
  \item{shift}{
  The number of bases that the negative strand is shifted towards its
  three prime end. This can be a vector, if the correlation should be
  calculated for different shifts.
}
  \item{bin}{
  If bin is larger than one, the coverage is calculated for bins of size
  \code{bin} and not for each single base. This speeds up calculations
  and might be beneficial in cases of low coverage. Note that shifting
  is performed after binning, so that the shift(s) should be multiple of
  \code{bin} (otherwise, it is rounded).
}
  \item{mode}{
  \code{mode} defines how bases (or bins) without reads are
  handled. \code{both} means that only bases covered on both strands
  are included when calculating the correlation. \code{one} means, that
  the base has to be covered on at least one strand whereas all bases
  are used independent of their coverage if \code{mode} is \code{none}.
}
  \item{minReads}{
  If not at least \code{minReads} are mapped to a chromosome, the
  chromosome is omitted.
}
  \item{chrs}{
  A \code{character} vector with the chromosomes that should be included
  into the calculation. \code{NA} means all chromosomes.
}
  \item{mc.cores}{
  Number of cores to be used.
}
}

\details{
  Only the 5 prime start positions of the reads are used for calculating
  the coverage. So, after removing duplicates in a single end sequencing
  experiment, the coverage can not be larger than one, if the bin size is
  set to one. (In this setting, mode \code{both} makes no sense.)
  If bin is larger than one, the coverage within a bin is aggregated.
  Then, the correlation is calculated for each shift. A shift
  (given in basepairs) should be multiple of the bin size
  (given in basepairs, too). If not, the binnend coverage is shifted by
  round(shift/bin) elements.

  The different modes define whether regions without coverage or with
  only one covered strand should used. The original implementation in
  the package \"spp\" does not make use of regions without
  coverage. However, this seems to be a loss of information, since no
  coverage has also a biological meaning in a ChIP-seq experiment. If
  the fragment size is approximately 500bp, setting \code{shift=seq(200, 800,
  10)}, \code{bin=10} and \code{mode="none"} should be a good setting.

  After the cross correlation was calculated for each chromosome,
  the weighted mean correlation across all chromosomes is
  calculated.  The weight for a specific chromosome equals the fraction
  of all reads that were aligned to that chromosome.
}

\value{
  A numeric vector with the cross correlation for each shift. The names
  of the vector corresponds to the shifts.
}

\references{
Kharchenko PV, Tolstorukov MY and Park PJ. Design and analysis of
ChIP-seq experiments for DNA-binding proteins. Nat Biotechnol 2008, 26(12):1351-9

Landt SG et al., ChIP-seq guidelines and practices of the ENCODE and
modENCODE consortia. Genome Res. 2012, 22(9):1813-31
}

\author{
Hans-Ulrich Klein (hklein@broadinstitute.org)
}

\seealso{
\code{\link{GRanges-class}}
}

\examples{
triangularKernel <- function(x, pos, h)  {
  res <- 1 - (abs(x - pos) / h)
  res[res < 0] <- 0
  return(res)
}
covPos <- round(triangularKernel(1:100, 60, 50) * 100)
covNeg <- round(triangularKernel(1:100, 65, 50) * 100)

reads <- GRanges(IRanges(start=c(rep(seq_along(covPos), covPos), rep(seq_along(covNeg), covNeg) - 9),
                         width=10),
                         strand=c(rep("+", sum(covPos)), rep("-", sum(covNeg))),
                         seqnames=rep("1", sum(covPos)+sum(covNeg)))

calculateCrossCorrelation(reads, shift=c(0,10), bin=1, mode="none", minReads=1)
cor(covPos, covNeg)
cor(covPos[1:(length(covPos)-10)], covNeg[11:length(covNeg)])

calculateCrossCorrelation(reads, shift=c(0,10), bin=1, mode="one", minReads=1)
cor(covPos[covPos != 0 | covNeg != 0], covNeg[covPos != 0 | covNeg != 0])

calculateCrossCorrelation(reads, shift=c(0,10), bin=1, mode="both", minReads=1)
cor(covPos[covPos != 0 & covNeg != 0], covNeg[covPos != 0 & covNeg != 0])


covPos2 <- round(triangularKernel(1:100, 60, 50) * 50)
covNeg2 <- round(triangularKernel(1:100, 68, 50) * 50)
reads2 <- GRanges(IRanges(start=c(rep(seq_along(covPos2), covPos2), rep(seq_along(covNeg2), covNeg2) - 9),
                          width=10),
                          strand=c(rep("+", sum(covPos2)), rep("-", sum(covNeg2))),
                          seqnames=rep("2", sum(covPos2)+sum(covNeg2)))
seqlevels(reads2) <- c("1", "2")
allReads <- c(reads, reads2)

calculateCrossCorrelation(allReads, shift=5, minReads=1, bin=1, mode="none")
cor1 <- cor(covPos[1:(length(covPos)-5)], covNeg[6:length(covNeg)])
cor2 <- cor(covPos2[1:(length(covPos2)-5)], covNeg2[6:length(covNeg2)])
cor1 * (sum(c(covPos, covNeg))/length(allReads)) +
  cor2 * (sum(c(covPos2, covNeg2))/length(allReads))
}

\keyword{ cross correlation }
