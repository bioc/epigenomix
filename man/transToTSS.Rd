\name{transToTSS}
\alias{transToTSS}
\docType{data}
\title{
  A data frame with Ensemble transcript IDs and transcriptional start sites.
}

\description{
  The data frame stores Ensemble transcript IDs and repective chromosomes,
  transcriptional start sites and strands for mus musculus (mm10).
}

\usage{data(transToTSS)}
\format{
  A data frame with 277 mouse transcripts with the following 4 variables:
  \describe{
    \item{\code{ensembl_transcript_id}}{A character giving the Ensemble
      transcript ID.}
    \item{\code{chromosome_name}}{A character with the respective
      chromomse name.}
    \item{\code{transcript_start}}{An integer storing the respective
      transcriptional start site.}
    \item{\code{strand}}{An integer storing the respective strand
      information.}
  }
}

\details{
  Given a character vector \code{transcripts} with the Ensemble
  transcript IDs, a data frame like this can be obtained via biomaRt:
  
  \code{
    library("biomaRt")
    mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    transToTSS <- getBM(attributes=c("ensembl_transcript_id",
      "chromosome_name", "transcript_start", "transcript_end", "strand"),
      filters="ensembl_transcript_id", values=transcripts, mart=mart)
  }
}

\source{
  http://www.ensembl.org
}

\examples{
data(transToTSS)
head(transToTSS)
}

\keyword{datasets}
