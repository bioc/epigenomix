.normalizeChIP <- function(object, method) {

  .Deprecated("normalize")

  if (!is.element(method, c("scaleTotal", "scaleRegion", "scaleMedianRegion", "quantile"))) {
    stop("Argument norm must be \"quantile\", \"scaleMedianRegion\", \"scaleRegion\" or \"scaleTotal\".")
  }

  geoMean <- function(x) {
    return(exp(mean(log(x))))
  }

  if (method == "scaleTotal") {
    if (!is.element("totalCount", colnames(colData(object)))) {
      stop("colData must contain the column totalCount.")
    }
    if (any(is.na(pData(object)$totalCount))) {
      stop("Column totalCount in ChIPseqSet must not contain NA.")
    }
    n <- colData(object)$totalCount
    s <- median(n) / n
    chipVals(object) <- t(t(chipVals(object)) * s)

  } else if (method == "scaleRegion") {
    n <- apply(chipVals(object), 2, sum)
    s <- median(n) / n
    chipVals(object) <- t(t(chipVals(object)) * s)

  } else if (method == "scaleMedianRegion") {
    meanSample <- apply(chipVals(object), 1, geoMean)
    s <- apply(chipVals(object), 2, function(x) {
        return(median(x/meanSample, na.rm=TRUE))
    })
    s <- 1/s
    chipVals(object) <- t(t(chipVals(object)) * s)

  } else if (method == "quantile") {
    ranksMat <- apply(chipVals(object), 2, order)
    readDist <- rep(0, nrow(object))
    for (i in 1:ncol(object)) {
      readDist <- readDist + chipVals(object)[ranksMat[, i], i]
    }
    readDist <- readDist / ncol(object)
    for (i in 1:ncol(object)) {
      chipVals(object)[ranksMat[, i], i] <- readDist
    }
    
  } else {
    stop("Argument norm must be \"quantile\", \"scaleMedianRegion\", \"scaleRegion\" or \"scaleTotal\".")
  }

  return(object)
}

setMethod("normalizeChIP",
    signature=c(object="ChIPseqSet", method="character"),
    .normalizeChIP)
