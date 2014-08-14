.summarizeReads <- function(object, regions, summarize) {

  sampleNames <- names(object)

  if (is.null(sampleNames) |
      length(unique(sampleNames)) != length(object) |
      any(is.na(sampleNames))) {
    stop("Names of argument object must be valid and unique sample names.")
  }

  if (is.null(names(regions))) {
    stop("Argument regions must have valid names.")
  }
  
  if (!is.element(summarize, c("add", "average"))) {
    stop("Argument summarize must be \"add\" or \"average\".")
  }

  counts <- NULL
  for (i in 1:length(sampleNames)) {
    sample <- object[[sampleNames[i]]]
    # countOverlaps is strand aware, so remove strand
    strand(sample) <- "*"
    co <- countOverlaps(regions, sample)
    counts <- cbind(counts, co)
  }

  rownames(counts) <- names(regions)
  colnames(counts) <- sampleNames

  if (summarize == "average" & is(regions, "GRangesList")) {
    n <- sapply(regions, length)
    counts <- apply(counts, 2, function(x) {x/n})
  }

  totalCount <- sapply(object[sampleNames], length)

  chipSet <- ChIPseqSet(chipVals=counts, rowData=regions,
      colData=DataFrame(totalCount=totalCount))

  return(chipSet)
}


setMethod("summarizeReads",
    signature=c(object="GRangesList", regions="GRangesList", summarize="character"),
    .summarizeReads)

setMethod("summarizeReads",
    signature=c(object="GRangesList", regions="GRanges", summarize="character"),
    .summarizeReads)

setMethod("summarizeReads",
    signature=c(object="GRangesList", regions="GRangesList", summarize="missing"),
    function(object, regions) {
      .summarizeReads(object=object, regions=regions, summarize="add")
    })

setMethod("summarizeReads",
    signature=c(object="GRangesList", regions="GRanges", summarize="missing"),
    function(object, regions) {
      .summarizeReads(object=object, regions=regions, summarize="add")
    })
