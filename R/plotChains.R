.plotChains <- function(object, chain, component, itb=1, thin=1, cols, ...) {

  compChains <- lapply(chains(object)$components, function(x) {
                     setdiff(names(x), "name")
                   })
  compNames <- sapply(chains(object)$components, function(x) {
                     x$name
                   })
  validChains <- c("pi", "dirichletParameter", "allocations", unique(unlist(compChains)))
  validComponents <- 1:length(components(object))

  if (!missing(chain) && length(chain) > 1) {
    stop("Argument chain must be of length one.")
  }
  if (!missing(chain) && !is.element(chain, validChains)) {
    stop(paste("Invalid chain argument. Must be one of" , paste(validChains, collapse=", ")))
  }
  if (!missing(component) && !all(is.element(component, validComponents))) {
    stop("Invalid component argument. Component with specified ID does not exist.")
  }
  
  # chain and component are given
  if (!missing(chain) & !missing(component)) {
    if (!(all(sapply(compChains[component], function(x) {is.element(chain, x)})) |
        (length(chain) == 1 & chain[1] == "pi"))) {
      stop("All specified chains must be present in all specified components.")
    }
  }
  # chain is missing and component is given
  if (missing(chain) & !missing(component)) {
    if (length(unique(sapply(chains(object)$components[component], function(x) {x$name}))) != 1) {
      stop("All specified components must be of the same type.")
    }
    chain <- compChains[[component[1]]]
  }
  # chain is given and component is missing
  if (!missing(chain) & missing(component)) {
    if (chain == ("pi")) {
      component <- 1:ncol(object)
    } else if (chain == "dirichletParameter") {
      component <- 0
    } else {
      component <- which(sapply(compChains, function(x) {is.element(chain, x)}))
    }
  }
  # component and chain are missing
  if (missing(chain) & missing(component)) {
    stop(paste("At least one of arguments chain and component must be specified.",
               "Valid arguments for chain are:", paste(validChains, collapse=", ")))
  }

  # extract data
  if (chain[1] == "pi") {
    data <- chains(object)[["pi"]]
    ind <- seq(itb+1, nrow(data), thin)
    data <- data[ind, component, drop=FALSE]
    ylab <- rep("Value", ncol(data))
    xlab <- rep("Iteration", ncol(data))
    main <- paste("pi - component ", component, sep="")
    if (missing(cols)) {
      cols <- 1
      rows <- ncol(data)
    } else {
      cols <- min(ncol(data), cols)
      rows <- ceiling(ncol(data)/cols)
    }
  } else if (chain[1] == "dirichletParameter") {
    data <- chains(object)[["dirichletParameter"]]
    ind <- seq(itb+1, length(data), thin)
    data <- matrix(data[ind], ncol=1)
    ylab <- "Value"
    xlab <- "Iteration"
    main <- "Dirichlet parameter"
    cols <- 1
    rows <- 1
  } else if (chain[1] == "allocations") {
    data <- chains(object)[["allocations"]]
    ind <- seq(itb+1, nrow(data), thin)
    data <- data[ind,,drop=FALSE]
    ylab <- rep("Value", ncol(data))
    xlab <- rep("Iteration", ncol(data))
    main <- paste("allocations - component ", 1:ncol(data), sep="")
    if (missing(cols)) {
      cols <- 1
      rows <- ncol(data)
    } else {
      cols <- min(ncol(data), cols)
      rows <- ceiling(ncol(data)/cols)
    }
  } else {
    data <- c()
    xlab <- character()
    ylab <- character()
    main <- character()
    for (i in component) {
      for (j in chain) {
        data <- cbind(data, chains(object)$components[[i]][[j]])
        main <- c(main, paste(j, " - component ", i, sep=""))
        ylab <- c(ylab, "Value")
        xlab <- c(xlab, "Iteration")
      }
    }
    ind <- seq(itb+1, nrow(data), thin)
    data <- data[ind,,drop=FALSE]
    if (missing(cols)) {
      cols <- length(chain)
      rows <- length(component)
    } else {
      cols <- min(ncol(data), cols)
      rows <- ceiling(ncol(data)/cols)
    }
  }
    
  args <- list(...)

  # default type is "l"
  if (!is.element("type", names(args))) {
    args$type <- "l"
  }

  # default main is set for each plot
  if (is.element("main", names(args))) {
    main <- args$main
    if (length(main) < ncol(data)) {
      main <- rep(main, ceiling(ncol(data)/length(main)))
    }
  }

  # call mfrow for mutliple plots
  if (cols > 1 | rows > 1) {
    par(mfrow=c(rows, cols))
  }

  for (i in 1:ncol(data)) {

    args$x <- 1:nrow(data)
    args$y <- data[,i]
    args$xlab <- xlab[i]
    args$ylab <- ylab[i]
    args$main <- main[i]

    do.call(plot, args)
  }
}

setMethod("plotChains", signature=c(object="MixModelBayes"), .plotChains)
