\name{fpkm}
\alias{fpkm}
\docType{data}
\title{
  Example RNA-seq data set.
}
\description{
  The \code{data.frame} stores transcription values obtained
  from the Cufflinks software for two samples (CEBPA_WT and
  CEBPA_KO). Transcription values are given in fragments per
  kilobase of transcripts per million fragments (FPKM).
}

\usage{data(fpkm)}

\format{
  An object of class \code{data.frame}.
}

\details{
  All transcripts sharing the TSS were grouped and one
  transcription values is given for each group of transcripts.
  The example data contains a subset of about 3500 TSS located on
  chromosome 1.
}

%%\source{
%%
%%}

\examples{
data(fpkm)
head(fpkm)
}

\keyword{datasets}
