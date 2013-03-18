.plotClassification <- function(object, method, ...) {

  args <- list(...)

  if (missing(method)) {
    method <- names(object@results$classification)[1]
  }

    
  # default colors
  if (!is.element("col", names(args))) {
    myCols <- sapply(components(object), function(x) {x@color})
  } else {
    myCols <- args$col
  }
  
  # set xlim and ylim if not given
  if (!is.element("xlim", names(args))) {
    args$xlim <- c(min(mmData(object)), max(mmData(object)))
  }
  if (!is.element("ylim", names(args))) {
    args$ylim <- c(1, length(components(object)))
  }

  # indices of used components
  compInds <- unique(classification(object, method))

  # first component - initial call of plots()
  if (!is.element("ylab", names(args))) {
    args$ylab <- "Mixture component"
  }
  if (!is.element("xlab", names(args))) {
    args$xlab <- "Data"
  }
  args$x <- mmData(object)[classification(object, method) == compInds[1]]
  args$y <- rep(compInds[1], length(args$x))
  args$col <- myCols[compInds[1]]
  do.call(plot, args)

  # plot other components
  if (length(compInds) > 1) {
    for (i in 2:length(compInds)) {
      x <- mmData(object)[classification(object, method) == compInds[i]]
      y <- rep(compInds[i], length(x))
      points(x, y, col=myCols[compInds[i]])
    }
  }
}

setMethod("plotClassification", signature(object="MixModel"),
          .plotClassification)
