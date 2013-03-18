.integrateData <- function(expr, chipseq, factor, reference) {

  if (!(is.element(factor, colnames(pData(expr))) &
      is.element(factor, colnames(pData(chipseq))))) {
    stop("Argument \"factor\" must be a pheno data column in \"expr\" and \"chipseq\".")
  }

  if (class(pData(expr)[, factor]) != "factor" | class(pData(chipseq)[, factor]) != "factor") {
    stop(paste("\"", factor, "\" must be of class factor in \"expr\" and \"chipseq\".", sep=""))
  }

  geFac <- pData(expr)[, factor]
  chipFac <- pData(chipseq)[, factor]

  if (length(levels(geFac)) != 2 | length(levels(chipFac)) != 2 |
      !all(is.element(levels(geFac), levels(chipFac)))) {
    stop(paste("\"", factor, "\" must have (the same) two levels in \"expr\" and \"chipseq\".", sep=""))
  }

  if (missing(reference)) {
    reference <- levels(geFac)[1]
  } else if (!is.element(reference, levels(geFac))) {
    stop("\"reference\" must be a level of \"factor\".")
  }

  ids <- intersect(featureNames(expr), featureNames(chipseq))
  if (length(ids) == 0) {
    stop("No commmon features.")
  }
  
  treatment <- setdiff(levels(geFac), reference)

  data <- cbind(apply(exprs(expr)[ids, geFac == treatment, drop=FALSE], 1, mean),
      apply(exprs(expr)[ids, geFac == reference, drop=FALSE], 1, mean),
      apply(chipVals(chipseq)[ids, chipFac == treatment, drop=FALSE], 1, mean),
      apply(chipVals(chipseq)[ids, chipFac == reference, drop=FALSE], 1, mean))

  deltaGe <- data[,1] - data[,2]
  deltaChip <- data[,3] - data[,4]
  
  z <- (deltaGe/sd(deltaGe)) * (deltaChip/sd(deltaChip))
  data <- cbind(data, z)

  colnames(data) <- c(paste("expr_", treatment, sep=""), paste("expr_", reference, sep=""),
       paste("chipseq_", treatment, sep=""), paste("chipseq_", reference, sep=""), "z")
  rownames(data) <- ids
  
  return(data)
}

setMethod("integrateData",
    signature=c(expr="ExpressionSet", chipseq="ChIPseqSet", factor="character", reference="character"),
    .integrateData)

setMethod("integrateData",
    signature=c(expr="ExpressionSet", chipseq="ChIPseqSet", factor="character", reference="missing"),
    .integrateData)

setMethod("integrateData",
    signature=c(expr="ExpressionSetIllumina", chipseq="ChIPseqSet", factor="character", reference="character"),
    .integrateData)

setMethod("integrateData",
    signature=c(expr="ExpressionSetIllumina", chipseq="ChIPseqSet", factor="character", reference="missing"),
    .integrateData)

