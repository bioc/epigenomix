.getAlignmentQuality <- function(bamFile, verbose=FALSE, mc.cores=1) {

  # internal function called via mclapply
  alnQualityMC <- function(bamFile, verbose) {  
    sampleName <- gsub(".bam$", "", basename(bamFile), ignore.case=TRUE)

    readGroup <- scanBamHeader(bamFile)[[1]]$text$'@RG'
    headerID <- gsub('^ID:', "", readGroup[grep('^ID:', readGroup)])
    headerSampleID <- gsub('^SM:', "", readGroup[grep('^SM:', readGroup)])
    headerLibraryID <- gsub('^LB:', "", readGroup[grep('^LB:', readGroup)])
    if (length(headerID) == 0) { headerID <- NA }
    if (length(headerSampleID) == 0) { headerSampleID <- NA }
    if (length(headerLibraryID) == 0) { headerLibraryID <- NA }
    
    if (verbose) {cat(paste("sample \"", sampleName, "\": Counting number of reads ... \n", sep=""))}
    param <- ScanBamParam()
    totalReads <- countBam(bamFile, param=param)$records

    if (verbose) {cat(paste("sample \"", sampleName, "\": Reading mapping qualities ... \n", sep=""))}
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
        what=c("rname", "pos", "strand", "mapq", "flag"))
    bam <- scanBam(bamFile, param=param)[[1]]
    mappedReads <- length(bam$pos)
    indMQ <- !(bam$mapq == 0)
    flags <- bam$flag[indMQ]
    umReads <- length(flags)
    indDup <- !bamFlagTest(flags, "isDuplicate")
    umuReads <- sum(indDup)
    mapq <- bam$mapq[indMQ][indDup]
    quantilesQual <- quantile(mapq, na.rm=TRUE)

    summary <- data.frame(
        Sample=sampleName,
        HeaderID=headerID,
        HeaderSampleID=headerSampleID,
        HeaderLibraryID=headerLibraryID,
        TotalReads=totalReads,
        MappedReads=mappedReads,
        MappedReadsRel=mappedReads/totalReads,
        UniquelyMappedReads=umReads,                    
        UniquelyMappedReadsRel=umReads/mappedReads,
        UniquelyMappedUniqueReads=umuReads,
        UniquelyMappedUniqueReadsRel=umuReads/mappedReads,
        NonRedundantFraction=umuReads/umReads,
        QualMean=mean(mapq),
        QualSd=sd(mapq),
        Quantile0=quantilesQual["0%"],
        Quantile25=quantilesQual["25%"],
        Quantile50=quantilesQual["50%"],
        Quantile75=quantilesQual["75%"],
        Quantile100=quantilesQual["100%"],
        Path=bamFile,
        stringsAsFactors=FALSE)
    
    return(summary)
  }
  
  res <- mclapply(bamFile, alnQualityMC, verbose=verbose, mc.preschedule=FALSE, mc.cores=mc.cores)
  res <- do.call(rbind, res)
  return(res)
}


setMethod("getAlignmentQuality",
    signature=c(bamFile="character"),
    .getAlignmentQuality)
