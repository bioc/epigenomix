\name{calculateCrossCorrelation}
\alias{calculateCrossCorrelation}


\title{
Calculate the cross correlation for a given GRanges object.
}
\description{
This method calculates the cross correlation, i.e. the pearson
correlation between the coverages of the positive and negative strand
from a DNA sequencing experiment. The cross correlation can be used as a
quality measure in ChIP-seq experiments (Kharchenko et al. 2008). Cross correlation
can also be used to estimate the fragment size by deteriming the shift
(given in base pairs) that maximizes the cross correlation. 
}

\usage{
calculateCrossCorrelation(object, srange = c(100, 500), bin = 10, minReads = 100000, mc.cores = 1)
}

\arguments{
  \item{object}{
  An \code{\link{GRanges}} object containing the aligned reads.
}
  \item{srange}{
  The number of bases that the negative strand is shifted towards its
  three prime end. This can be a vectore, if the correlation should be
  calculated for different shifts.
}
  \item{bin}{
  A \code{character} vector with the chromosomes that should be included
  into the calculation. \code{NA} means all chromosomes. 
}
  \item{minReads}{
  If not at least \code{minReads} are mapped to a chromosome, the
  chromosome is omitted.
}
  \item{mc.cores}{
  Number of cores to be used.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kharchenko PV, Tolstorukov MY and Park PJ. Design and analysis of
ChIP-seq experiments for DNA-binding proteins. Nat Biotechnol 2008, 26(12):1351-9 
}
\author{
Hans-Ulrich Klein (hklein@broadinstitute.org)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
triangularKernel <- function(x, pos, h)  {
  res <- 1 - (abs(x - pos) / h)
  res[res < 0] <- 0
  return(res)
}
covPos <- round(triangularKernel(1:100, 60, 50) * 100)
covNeg <- round(triangularKernel(1:100, 65, 50) * 100)

reads <- GRanges(IRanges(start=c(rep(seq_along(covPos), covPos), rep(seq_along(covNeg), covNeg) - 
9),
                         width=10),
                         strand=c(rep("+", sum(covPos)), rep("-", sum(covNeg))),
                         seqnames=rep("1", sum(covPos)+sum(covNeg)))

calculateCrossCorrelation(reads, shift=c(0,10), bin=1, minReads=1)
cor(covPos, covNeg)

# zeros are added at the end
calculateCrossCorrelation(reads, shift=c(4,6), bin=1, minReads=1)
cor(covPos, c(covNeg[6:length(covNeg)], rep(0, 5)))

# leading zeros are added at both sides
which(covNeg != 0)[1]
calculateCrossCorrelation(reads, shift=c(5,20), bin=10, minReads=1)
cor(c(rep(0, 5), covPos), c(covNeg[16:length(covNeg)], rep(0, 20)))

#another chromosome
covPos2 <- round(triangularKernel(1:100, 60, 50) * 50)
covNeg2 <- round(triangularKernel(1:100, 68, 50) * 50)

#reads2 <- GRanges(IRanges(start=c(rep(seq_along(covPos2), covPos2), rep(seq_along(covNeg2), 
#covNeg2) - 9),
#                         width=10),
#                         strand=c(rep("+", sum(covPos2)), rep("-", sum(covNeg2))),
#                         seqnames=rep("2", sum(covPos2)+sum(covNeg2)))
#seqlevels(reads2) <- c("1", "2")
#allReads <- c(reads, reads2)

#calculateCrossCorrelation(allReads, shift=5)
#cor1 <- cor(covPos, c(covNeg[6:length(covNeg)], rep(0, 5)))
#cor2 <- cor(covPos2, c(covNeg2[6:length(covNeg2)], rep(0, 5)))
#cor1 * (sum(c(covPos, covNeg))/length(allReads)) + cor2 * (sum(c(covPos2, 
#covNeg2))/length(allReads))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
