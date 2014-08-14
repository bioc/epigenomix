\name{summarizeReads}
\alias{summarizeReads}
\alias{summarizeReads,GRangesList,GRanges,character-method}
\alias{summarizeReads,GRangesList,GRanges,missing-method}
\alias{summarizeReads,GRangesList,GRangesList,character-method}
\alias{summarizeReads,GRangesList,GRangesList,missing-method}


\title{
  Count reads lying within given regions.
}

\description{
  This function takes reads from e.g. ChIP-seq experiments and regions,
  e.g. promoters of genes, and assigns the number of overlapping reads
  to that region. This method was written particularly with regard to
  histone ChIP-seq experiments. Some histone modifications mainly occur
  at transcriptional start sites and thus ChIP-seq values can be assigned
  to genes by counting the number of reads within genes' pomoter
  regions. However, some genes may have several transcript and hence
  several promoters. Different options for handling multiple promoters
  are implemented. This method is also useful when integrating microarray
  expression data and ChIP-seq data, since most array platforms are gene
  centric and have probes that measure several transcripts.
}

\usage{
  summarizeReads(object, regions, summarize)
}

\arguments{
  \item{object}{
    A \code{GRangesList} with one \code{GRanges} object for each sample
    storing the ChIP-seq reads. The names of the \code{GRangesList}
    elements are used as sample names.
  }
  \item{regions}{
    An object of type \code{GRangesList} storing the promoter regions.
    Each element can be interpreted as gene or probe that has one or
    more promoters. The names of the lists' elements are used as
    row names.
    Alternatively, regions can be a \code{GRanges} object which as then handled
    like a \code{GRangesList} object with only one region in each list element.
    names of the \code{GRanges} are used as row names in this case.
  }
  \item{summarize}{
    Defines how regions with several ranges are handled. "average"
    means that the mean count of reads across all ranges is
    assigned to the region whereas "add" means that all counts are
    simply added (default).
  }
}

\details{
  This function is usually applied after calling
  \code{\link{matchProbeToPromoter}}. When
  \code{\link{matchProbeToPromoter}} is used with \code{mode} "union",
  it is recommended to use "add". If the option "keepAll"
  had been used, one might want to use "average".

  This method uses \code{\link{countOverlaps}} and counts
  each read that overlaps with at least one base. 
%  If more flexibility is needed, \code{\link{summarizeOverlaps}}
%  might be an option. In contrast to this method,
%  \code{\link{summarizeOverlaps}} ensures that each
%  read is counted only once, even if it overlaps with different regions.
}

\value{
  Returns a \code{ChIPseqSet} with number of rows equal to the length
  of \code{regions} and number of samples equal to the length of
  \code{object}.
}

\author{
  Hans-Ulrich Klein (h.klein@uni-muenster.de)
}

\seealso{
  \code{\link{matchProbeToPromoter}}
  \code{\link{ChIPseqSet-class}}
%  \code{\link[GenomicRanges]{summarizeOverlaps}}
}

\examples{
chipSeq <- GRangesList()
chipSeq[[1]] <- GRanges(seqnames=c("1", "1", "1", "1"),
    ranges=IRanges(start=c(97900, 198200, 198600, 202500),
        end=c(98100, 198400, 198800, 202700)),
    strand=c("+", "+", "+", "+"))
chipSeq[[2]] <- GRanges(seqnames=c("1", "1", "1", "1"),
    ranges=IRanges(start=c(97900, 198200, 198600, 300000),
        end=c(98100, 198400, 198800, 300200)),
    strand=c("+", "+", "+", "+"))
names(chipSeq) = c("sample1", "sample2")

promoters <- GRanges(seqnames=c("1", "1", "1"),
    ranges=IRanges(start=c(98000, 198000, 202000),
        end=c(101999, 201999, 205999)),
    strand=c("-", "+", "+"),
    probe=c("101", "102", "102"))
promoters <- split(promoters, elementMetadata(promoters)$probe)

chipSet <- summarizeReads(chipSeq, promoters, summarize="add")
chipVals(chipSet)
}
