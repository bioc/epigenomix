setValidity("ChIPseqSet", function(object) {
  msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("chipVals")))
  if (!is.element("totalCount", colnames(pData(object)))) {
    msg <- validMsg(msg,
        "totalCount must be a pheno data column in ChIPseqSet")
  } else {
    if (!is(pData(object)$totalCount, "numeric")) {
      msg <- validMsg(msg,
        "Pheno data column totalCount in ChIPseqSet must be numeric")
    }
  }    
  if (is.null(msg)) TRUE else msg
})

setMethod("chipVals", signature(object="ChIPseqSet"),
          function(object) assayDataElement(object, "chipVals"))

setReplaceMethod("chipVals", signature(object="ChIPseqSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "chipVals", value))
