setMethod("ChIPseqSet", signature(chipVals="matrix", rowRanges="GRanges"),
          function(chipVals, rowRanges, colData=DataFrame(row.names=colnames(chipVals)),
                   metadata=list(), ...) {
              ssla <- new("ShallowSimpleListAssays", data=SimpleList(chipVals=chipVals))
              new("ChIPseqSet",
                  assays = ssla,
                  rowRanges = rowRanges,
                  colData = colData,
                  metadata = as.list(metadata))
          }
)

setMethod("ChIPseqSet", signature(chipVals="matrix", rowRanges="GRangesList"),
          function(chipVals, rowRanges, colData=DataFrame(row.names=colnames(chipVals)),
                   metadata=list(), ...) {
              ssla <- new("ShallowSimpleListAssays", data=SimpleList(chipVals=chipVals))
              new("ChIPseqSet",
                  assays = ssla,
                  rowRanges = rowRanges,
                  colData = colData,
                  metadata = as.list(metadata))
          }
)


setMethod("chipVals", signature(object="ChIPseqSet"),
          function(object){
            return(assays(object)$chipVals)
          }
)

setReplaceMethod("chipVals", signature(object="ChIPseqSet", value="matrix"),
                 function(object, value){
                     assays(object)$chipVals <- value
                     return(object)
                 }
)

setMethod("cpm", signature(object="ChIPseqSet"),
  function(object, libSize, log2=FALSE, priorCount=0.1) {
      if (missing(libSize)) {
        libSize <- apply(chipVals(object), 2, sum)
      }
      if (length(libSize) != ncol(object)) {
        stop("Length of libSize must equal number of samples.")
      }
      if (log2) {
	priorCountScaled <- libSize / mean(libSize) * priorCount
	libSize <- libSize + (nrow(object) * priorCountScaled)
      }
      if (log2) {
        chipVals(object) <- log2(t( (t(chipVals(object)) + priorCountScaled) * (1000000 / libSize) ))
      } else {
        chipVals(object) <- t(t(chipVals(object)) * (1000000 / libSize))
      }
      return(object)
})

setValidity("ChIPseqSet", function(object){
  if(length(assays(object)) != 1)
    return("The assays slot in ChIPseqSet objects must be of length one.")
  if(!is(assay(object), "matrix"))
    return("The chipValues slot of a ChIPseqSet object must be a matrix.")
  if(mode(assay(object)) != "numeric")
    return("The chipValues matrix of a ChIPseqSet object must contain numeric data.")
})
