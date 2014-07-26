\name{getAlignmentQuality}
\alias{getAlignmentQuality}
\alias{getAlignmentQuality,character-method}

\title{Calculation of basic alignments statistics}

\description{
Calculates some basic alignment statistics for given bam files.
}

\usage{
getAlignmentQuality(bamFile, verbose = FALSE, mc.cores = 1)
}

\arguments{
  \item{bamFile}{
  A \code{character} vector with the filenames of the bam files
}
  \item{verbose}{
  If set to \code{TRUE}, some status information is written to the
  R console.
}
  \item{mc.cores}{
  Number of cores to be used.
}
}

\details{
  The given bam files should have marked duplicates and not uniquely
  mapped reads should have a quality value of 0. In detail, this
  function returns a data frame with the following columns:
  \describe{
    \item{Sample}{File name without path and suffix}
    \item{HeaderID}{ID field from bam header, if available}
    \item{HeaderSampleID}{SM field from bam header, if available}
    \item{HeaderLibraryID}{LB field from bam header, if available}
    \item{TotalReads}{Total number of reads in bam file}
    \item{MappedReads}{Number of mapped Reads}
    \item{MappedReadsRel}{MappedReads/TotalReads}
    \item{UniquelyMappedReads}{Number of mapped reads with mapping quality larger 0}
    \item{UniquelyMappedReadsRel}{UniquelyMappedReads/MappedReads}
    \item{UniquelyMappedUniqueReads}{Number of non duplicated mapped reads with mapping quality larger 0}
    \item{UniquelyMappedUniqueReadsRel}{UniquelyMappedUniqueReads/MappedReads}
    \item{NonRedundantFraction}{UniquelyMappedUniqueReads/UniquelyMappedReads}
    \item{QualMean}{Mean mapping quality of all uniquely mapped unique reads}
    \item{QualSd}{Standard deviation of the mapping quality of all uniquely mapped unique reads}
    \item{Quantile0}{0\% quantileof the mapping quality of all uniquely mapped unique reads}
    \item{Quantile25}{25\% quantile of the mapping quality of all uniquely mapped unique reads}
    \item{Quantile50}{50\% quantile of the mapping quality of all uniquely mapped unique reads}
    \item{Quantile75}{75\% quantile of the mapping quality of all uniquely mapped unique reads}
    \item{Quantile100}{100\% quantile of the mapping quality of all uniquely mapped unique reads}
    \item{Path}{Full path and file name as given in argument \code{bamFile}}
  }
}

\value{
Returns a data frame with one row for each given bam file and the
columns as listed in the details section.
}

\author{
Hans-Ulrich Klein (hklein@broadinstitute.org)
}

\examples{
\dontrun{getAlignmentQuality("myFile.bam")}
}
