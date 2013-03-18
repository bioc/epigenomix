.matchProbeToPromoter <- function(probeToTranscript, transcriptToTSS, promWidth, mode) {

  if (is.null(names(probeToTranscript)) || any(duplicated(names(probeToTranscript)))) {
    stop("names of probeToTranscript must be valid unique probe names.")
  }

  if (ncol(transcriptToTSS) < 4) {
    stop("transcriptToTSS must have four columns with IDs, chromosome, strand and TSS position.")
  }
  
  # remove all probes without transcript annotation
  validP2T <- probeToTranscript[!sapply(probeToTranscript, function(transcripts) {
    length(transcripts) == 0 || is.na(transcripts[1])
  })]
  probes <- rep(names(validP2T), sapply(validP2T, length))
  validP2T <- unlist(validP2T)

  # match all probes to the row number in the transcriptToTSS data frame and
  # remove not matchable transcripts/probes
  matchedP2T <- match(validP2T, transcriptToTSS[, 1])
  indNA <- is.na(matchedP2T)
  transcriptDf <- transcriptToTSS[matchedP2T[!indNA], ]
  transcriptDf[,5] <- probes[!indNA]
  
  # keepAll: just keep all promoter regions of a probe
  if (mode == "keepAll") {
    ranges <- GRanges(seqnames=transcriptDf[, 2],
        ranges=IRanges(start=transcriptDf[, 3], width=1),
        strand=transcriptDf[, 4],
        probe=transcriptDf[, 5])
    ranges <- resize(ranges, width=promWidth, fix="center")


  # dropMultiple: remove all probes with > 1 unique promoter
  } else if (mode == "dropMultiple") {
    dupPromoters <- duplicated(transcriptDf[, -1])
    transcriptDf <- transcriptDf[!dupPromoters, ]
    ambiguousProbes <- unique(transcriptDf[duplicated(transcriptDf[, 5]), 5])
    transcriptDf <- transcriptDf[!is.element(transcriptDf[, 5], ambiguousProbes), ]
    ranges <- GRanges(seqnames=transcriptDf[, 2],
        ranges=IRanges(start=transcriptDf[, 3], width=1),
        strand=transcriptDf[, 4],
        probe=transcriptDf[, 5])
    ranges <- resize(ranges, width=promWidth, fix="center")

    
  # union: calculate the union of all promoter regions of a probe
  } else if (mode == "union") {
    dupPromoters <- duplicated(transcriptDf[, -1])
    transcriptDf <- transcriptDf[!dupPromoters, ]

    ambiguousProbes <- unique(transcriptDf[duplicated(transcriptDf[, 5]), 5])
    tDfClear <- transcriptDf[!is.element(transcriptDf[, 5], ambiguousProbes), ]
    tDfAmbiguous <- transcriptDf[is.element(transcriptDf[, 5], ambiguousProbes), ]
    tDfAmbiguous <- split.data.frame(tDfAmbiguous, tDfAmbiguous[, 5])

    rangesClear <- GRanges(seqnames=tDfClear[, 2],
        ranges=IRanges(start=tDfClear[, 3], width=1),
        strand=tDfClear[, 4],
        probe=tDfClear[, 5])
    rangesClear <- resize(rangesClear, width=promWidth, fix="center")

    rangesAmbiguous <- lapply(tDfAmbiguous, function(df) { # this takes some time :(
        ranges = GRanges(seqnames=df[, 2],
          ranges=IRanges(start=df[, 3], width=1),
          strand=df[, 4])
        ranges = resize(ranges, width=promWidth, fix="center")
        ranges = reduce(disjoin(ranges))
        elementMetadata(ranges)$probe = df[1, 5]
        return(ranges)
    })
    names(rangesAmbiguous) <- NULL
    rangesAmbiguous <- suppressWarnings(do.call(c, rangesAmbiguous))
    ranges <- c(rangesClear, rangesAmbiguous)

    
  } else {
    stop("Argument mode must be \"union\", \"keepAll\" or \"dropMultiple\".")
  }

  ranges <- sort(ranges)
  rangesList <- split(ranges, elementMetadata(ranges)$probe)

  return(rangesList)
}


setMethod("matchProbeToPromoter",
          signature=c(probeToTranscript="list", transcriptToTSS="data.frame",
            promWidth="numeric", mode="character"),
          .matchProbeToPromoter)

setMethod("matchProbeToPromoter",
          signature=c(probeToTranscript="list", transcriptToTSS="data.frame",
            promWidth="missing", mode="character"),
          function(probeToTranscript, transcriptToTSS, mode) {
            .normalizeChIP(probeToTranscript, transcriptToTSS, 4000, mode)
          })

setMethod("matchProbeToPromoter",
          signature=c(probeToTranscript="list", transcriptToTSS="data.frame",
            promWidth="numeric", mode="missing"),
          function(probeToTranscript, transcriptToTSS, promWidth) {
            .normalizeChIP(probeToTranscript, transcriptToTSS, promWidth, "union")
          })

setMethod("matchProbeToPromoter",
          signature=c(probeToTranscript="list", transcriptToTSS="data.frame",
            promWidth="missing", mode="missing"),
          function(probeToTranscript, transcriptToTSS) {
            .normalizeChIP(probeToTranscript, transcriptToTSS, 4000, "union")
          })
