setMethod("ChIPseqSet", signature(chipVals="matrix", rowData="GRanges"),
          function(chipVals, rowData, colData=DataFrame(row.names=colnames(chipVals)),
                   exptData=SimpleList(), ...) {
              ssla <- new("ShallowSimpleListAssays", data=SimpleList(chipVals=chipVals))
              new("ChIPseqSet",
                  assays = ssla,
                  rowData = rowData,
                  colData = colData,
                  exptData = exptData)
          }
)

setMethod("ChIPseqSet", signature(chipVals="matrix", rowData="GRangesList"),
          function(chipVals, rowData, colData=DataFrame(row.names=colnames(chipVals)),
                   exptData=SimpleList(), ...) {
              ssla <- new("ShallowSimpleListAssays", data=SimpleList(chipVals=chipVals))
              new("ChIPseqSet",
                  assays = ssla,
                  rowData = rowData,
                  colData = colData,
                  exptData = exptData)
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


setValidity("ChIPseqSet", function(object){
  if(length(assays(object)) != 1)
    return("The assays slot in ChIPseqSet objects must be of length one.")
  if(!is(assay(object), "matrix"))
    return("The chipValues slot of a ChIPseqSet object must be a matrix.")
  if(mode(assay(object)) != "numeric")
    return("The chipValues matrix of a ChIPseqSet object must contain numeric data.")
})
