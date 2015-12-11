.normalize.scale <- function(mat, n) {
    s <- median(n) / n
    mat <- t(t(mat) * s)
    return(mat)
}

.normalize.scaleMedianRegion <- function(mat) {
    geoMean <- function(x) { return(exp(mean(log(x)))) }
    meanSample <- apply(mat, 1, geoMean)
    s <- apply(mat, 2, function(x) {
        return(median(x/meanSample, na.rm=TRUE))
    })
    s <- 1/s
    mat <- t(t(mat) * s)
    return(mat)
}

.normalize.quantile <- function(mat) {
    ranksMat <- apply(mat, 2, order)
    readDist <- rep(0, nrow(mat))
    for (i in 1:ncol(mat)) {
      readDist <- readDist + mat[ranksMat[, i], i]
    }
    readDist <- readDist / ncol(mat)
    for (i in 1:ncol(mat)) {
      mat[ranksMat[, i], i] <- readDist
    }
    return(mat)
}

.normalize.tmm <- function(mat, trim) {
    geoMean <- function(x) { return(exp(mean(log(x)))) }
    meanSample <- apply(mat, 1, geoMean)
    indNull <- meanSample == 0
    s <- apply(log(mat[!indNull, ]/meanSample[!indNull]), 2, mean, trim=trim)
    mat <- t(t(mat) / exp(s))
    return(mat)
}


.normalize <- function(object, method, isLogScale=FALSE, trim=0.3, totalCounts, ...) {

    if (!is.element(method, c("scale", "scaleMedianRegion", "quantile", "tmm"))) {
        stop("Argument method must be \"quantile\", \"scaleMedianRegion\", \"scale\" or \"tmm\".")
    }

    
    if (is(object, "ChIPseqSet")) {
        mat <- chipVals(object)
    }
    if (is(object, "ExpressionSet")) {
        mat <- exprs(object)
    }    
    if (isLogScale) {
        mat <- exp(mat)
    }

    
    if (method == "scale") {
        if (missing(totalCounts)) {
            n <- apply(mat, 2, sum)
        } else {
            if (totalCounts != ncol(mat)) {
                stop("Length of \"totalCounts\" must equal number of samples.")
            }
        }
        mat <- .normalize.scale(mat, n)


    } else if (method == "scaleMedianRegion") {
        mat <- .normalize.scaleMedianRegion(mat)


    } else if (method == "quantile") {
        mat <- .normalize.quantile(mat)

    } else if (method == "tmm") {
        if (trim < 0 | trim > 1) {
            stop("Argument \"trim\" must be in [0, 1].")
        }
        mat <- .normalize.tmm(mat, trim=trim)
    }
    

    if (isLogScale) {
        mat <- log(mat)
    }
    if (is(object, "ChIPseqSet")) {
        chipVals(object) <- mat
    }
    if (is(object, "ExpressionSet")) {
        exprs(object) <- mat
    }
    
    return(object)
}

setMethod("normalize",
    signature=c(object="ChIPseqSet"),
    .normalize)


setMethod("normalize",
    signature=c(object="ExpressionSet"),
    .normalize)
