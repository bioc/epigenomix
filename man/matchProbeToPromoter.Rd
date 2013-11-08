\name{matchProbeToPromoter}
\alias{matchProbeToPromoter}
\alias{matchProbeToPromoter,list,data.frame-method}

\title{
  A function assigning promoter regions to given probe IDs.
}

\description{
  This function returns a \code{GRangesList} object asigning promoter regions
  to probes. The assignment of transcripts to probes and the
  transcriptional start sites must be given as arguments.
}

\usage{
matchProbeToPromoter(probeToTranscript, transcriptToTSS, promWidth = 4000, mode = "union", fix = "center")
}

\arguments{
  \item{probeToTranscript}{
    A \code{list} with character vectors as elements. The elements' names are
    probe IDs and the character vectors store the transcript IDs assigned to
    that probe.
  }
  \item{transcriptToTSS}{
    A \code{data.frame} with four columns:
    \enumerate{
      \item Transcript ID as given in the argument \code{probeToTranscript}
      \item Chromosome
      \item Transcriptional start site in base pairs
      \item Strand
    }
  }
  \item{promWidth}{
    Width of the promoter regions in base pairs. Promoters are defined
    as \code{promWidth} base pairs upstream of the transcriptional start
    site. (default 4000bp)
  }
  \item{mode}{
    How probes with multiple transcripts should be handled. Must be either 
    "union", "keepAll" or "dropMultiple". (default "union")
  }
  \item{fix}{
    Denotes what to use as anchor when defining the promoter
    region. Must be either  "center", "start" or "end". "Center" means
    that the TSS is in the middle of the promoter, whereas "end" means
    that the promoter is placed upsream of the TSS. (default "center")
  }
}

\details{
  More than one transcript can be assigned to one probe in the given
  \code{probeToTranscript} argument. Several options how to handle
  such cases can be choosen by argument \code{mode}. "union":
  The union of all promoters is calculated and assigned to the probe.
  "keepAll": All promoters of all transcripts are assigned to the
  probe. If some transcript have identical TSSs, the same promoter
  region occurs several times. "dropMultiple": All probes that have
  more than one transcript with different TSS are removed.

  The argument \code{transcriptToTSS} must have at least 4 columns
  giving the information as described above. The column names are not
  decisive, but their position.
}

\value{
  An object of class \code{GRangesList} with one element for each probe.
  If \code{mode} is not set to "dropMultiple", \code{GRanges} may consist
  of more than one range. The names of the lists' elements are the probe
  IDs and additionally, each \code{GRanges} has a meta data column
  "probe" giving the corresponding probe ID.
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\seealso{
  \code{\link{summarizeReads}}
}

\examples{
probeToTrans <- list("101"="ENST00011", 
                     "102"=c("ENST00021", "ENST00022"),
                     "103"=NA)
transToTSS <- data.frame(
    transID=c("ENST00011", "ENST00021", "ENST00022"),
    chr=c("1", "1", "1"),
    tss=c(100000, 200000, 201000),
    strand=c("-", "+", "+"))

matchProbeToPromoter(probeToTrans, transToTSS,
    promWidth=4000, mode="union")
matchProbeToPromoter(probeToTrans, transToTSS,
    promWidth=4000, mode="keepAll")
}
