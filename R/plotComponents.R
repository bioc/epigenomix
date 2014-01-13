.plotComponents <- function(object, density=FALSE, ...) {

  # number of points / resolution
  numPoints <- 2048

  # default colors
#  myDefCols <- c("normNull"="blue", "expNeg"="red", "expPos"="green", "gammaNeg"="yellow", "gammaPos"="brown")

  # set xlim if not given
  args <- list(...)
  if (!is.element("xlim", names(args))) {
    args$xlim <- c(floor(-max(abs(mmData(object)))), ceiling(max(abs(mmData(object)))))
  }

  # calculate data distribution, if density=TRUE
  if (density) {
    if (is.element("adjust", names(args))) {
      d <- density(mmData(object), n=numPoints, from=args$xlim[1], to=args$xlim[2], adjust=args$adjust)
      args$adjust = NULL
    } else {
      d <- density(mmData(object), n=numPoints, from=args$xlim[1], to=args$xlim[2])
    }
    xp <- d$x
    yData <- d$y
  } else {
    xp <- seq(from=args$xlim[1], to=args$xlim[2], length.out=numPoints)
  }
  
  # calculate weighted components' distributions
  yComps <- matrix(nrow=length(components(object)), ncol=numPoints)
  myCols <- rep(NA, length(components(object)))
  for (i in 1:length(components(object))) {
    yComps[i, ] <- weights(object)[i] * do.call(components(object)[[i]]@pdf, c(list(x=xp), components(object)[[i]]@parameters))
    myCols[i] <- components(object)[[i]]@color
    # for some distributions, it looks better to not plot beyond 0 
    if (components(object)[[i]]@name == "ExpNeg" | components(object)[[i]]@name == "GamNeg") {
      yComps[i, xp > 0] <- NA
    } else if (components(object)[[i]]@name == "ExpPos" | components(object)[[i]]@name == "GamPos") {
      yComps[i, xp < 0] <- NA
    }
  }
  
  # calculate mixture distribution
  yMix <- apply(yComps, 2, sum, na.rm=TRUE)

  # set ylim if not given
  if (!is.element("ylim", names(args))) {
    if (density) {
      args$ylim <- c(0, max(yMix, yData)*1.05)
    } else {
      args$ylim <- c(0, max(yMix)*1.05)
    }
  }

  # set main if not given
  if (!is.element("main", names(args))) {
    args$main <- ""
  }
  
  # set colors (mix, 1st comp, ..., nth comp, data)
  myCols <- c("black", myCols, "orange")
  if (is.element("col", names(args))) {
    myCols <- args$col
  }

  # set line types (mix, 1st comp, ..., nth comp, data)
  myLtys <- rep(1, length(myCols))
  if (is.element("lty", names(args))) {
    myLtys <- args$lty
  }

  # set xlab and ylab if not given
  if (!is.element("xlab", names(args))) {
    args$xlab <- "Data"
  }
  if (!is.element("ylab", names(args))) {
    args$ylab <- "Density"
  }

  # plot densities
  if (density) { # plot empirical density 
    args$col <- myCols[length(myCols)]
    args$type <- "l"
    args$lty <- myLtys[length(myLtys)]
    args$x <- xp
    args$y <- yData
    do.call(plot, args)
  } else {       # plot histogram of data
    args$x <- mmData(object)
    args$freq <- FALSE
    args$col <- "gray"
    if (!is.element("border", names(args))) {
       args$border <- FALSE
    }
    if (!is.element("breaks", names(args))) {
       args$breaks <- ((max(args$x) - min(args$x)) / (args$xlim[2] - args$xlim[1])) * 100
    }
    do.call(hist, args)
  }
  lines(xp, yMix, col=myCols[1], lty=myLtys[1]) # plot mixture density
  # plot null components first, then neg. and pos. components
  indFirst <- which(!is.element(sapply(components(object), function(x) {return(x@name)}),
                                c("ExpNeg", "ExpPos", "GamNeg", "GamPos")))
  indSecond <- setdiff(1:length(components(object)), indFirst)
  for (i in indFirst) {   # plot null components
    lines(xp, yComps[i, ], col=myCols[1+i], lty=myLtys[1+i])
  }
  for (i in indSecond) {     # plot neg + pos components
    lines(xp, yComps[i, ], col=myCols[1+i], lty=myLtys[1+i])
  }
}

setMethod("plotComponents", signature=c(object="MixModel"),
           .plotComponents)
