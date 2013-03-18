\name{mappedReads}
\alias{mappedReads}
\docType{data}
\title{
  Mapped reads obtained from a anti-histone ChIP-seq experiment.
}

\description{
  The \code{\linkS4class{GRangesList}} contains two elements:
  "CEBPA_WT_1" and "CEBPA_KO_1". Both list elements are
  \code{\linkS4class{GRanges}} objects storing mapped reads
  from anti-H3K4me3 ChIP-seq experiments. The first sample was a
  wild-type mouse cell line. The second sample was obtained from
  the same cell line after CEPBA knock-out.
}

\usage{data(mappedReads)}
\format{
  A \code{\linkS4class{GRangesList}} with two
  \code{\linkS4class{GRanges}}.
}

\details{
  Duplicated reads and reads mapping to more than one genomic location
  were removed. Reads were extended to the estimated DNA fragment size
  of 200bp towards the 3 prime end. Further, only reads lying within
  certain regions of chromomse 1 were kept to reduce storage space.
}

%%\source{
%% 
%%}

\examples{
data(mappedReads)
names(mappedReads)
mappedReads[[1]]
}

\keyword{datasets}
