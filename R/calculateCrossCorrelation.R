.calculateCrossCorrelation <- function(object, shift=c(200, 250, 300), bin=10, mode="none", minReads=10000, chrs=NA, mc.cores=1) {

  if (any(shift < 0)) {stop("shift musst be larger or equal 0")}

  # Strand information must be present in the given object.
  n <- sum(strand(object) == "*")
  if (n > 0) {
    if (n == length(object)) {
      stop("Strand information is missing for all reads.")
    } else {
      warning(paste("Strand information is missing for", n, "reads."))
    }
  }
  
  object <- resize(object, width=1)
  object <- split(object, seqnames(object))
  if (!is.na(chrs[1])) {
    object <- object[is.element(names(object), chrs)]
  }
  numReads <- sapply(object, length)
  object <- object[numReads >= minReads]
  numReads <- numReads[numReads >= minReads]

  chrCor <- mclapply(object, function(r) {
      r <- keepSeqlevels(r, value=as.character(seqnames(r)[1]))
      maxPos <- max(start(r))
      indPlus <- strand(r) == "+"
      shiftName <- as.character(shift)
      
      positive <- coverage(r[indPlus], width=maxPos)[[1]]
      negative <- coverage(r[!indPlus], width=maxPos)[[1]]

      if (bin > 1) {
        views <- Views(positive, start=seq(1, length(positive)-bin, bin), width=bin)
        positive <- Rle(viewSums(views))
        views <- Views(negative, start=seq(1, length(negative)-bin, bin), width=bin)
        negative <- Rle(viewSums(views))
        shift <- round(shift/bin)
      }

      cors <- numeric()
      for (i in 1:length(shift)) {
        positiveShift <- window(positive, start=1, end=length(positive)-shift[i])
        negativeShift <- window(negative, start=shift[i]+1, end=length(negative))
        if (mode == "none") {          # use all bins/bases
          cors[i] <- cor(positiveShift, negativeShift)
        } else if (mode == "one") {    # use only bins/bases covered on at least one strand
          ind <- !(positiveShift == 0 & negativeShift == 0)
          cors[i] <- cor(positiveShift[ind], negativeShift[ind])
        } else if (mode == "both") {   # use only bins/bases covered on both strands
          ind <- !(positiveShift == 0 | negativeShift == 0)
          cors[i] <- cor(positiveShift[ind], negativeShift[ind])
        } else {
          stop("Argument mode must be \"none\", \"one\" or \"both\".")
        }
      }

      names(cors) <- shiftName
      return(cors)
  }, mc.cores=mc.cores)

  cors <- do.call(cbind, chrCor)
  weights <- numReads / sum(numReads)
  cors <- apply(t(cors) * weights, 2, sum)
    
  return(cors)
}


setMethod("calculateCrossCorrelation",
    signature=c(object="GRanges"),
    .calculateCrossCorrelation)
