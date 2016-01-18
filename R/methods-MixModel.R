# methods for class MixModel
setMethod("classification", signature(object="MixModel", method="missing"),
          function(object) {
            object@results$classification[[1]]
          })

setMethod("classification", signature(object="MixModel", method="character"),
          function(object, method) {
            if (is.element(method, names(object@results$classification))) {
              object@results$classification[[method]]
            } else {
              stop(paste("A classification based on ", method, " does not exist for the given object.", sep=""))
            }
          })

setMethod("components", signature(object="MixModel"),
          function(object) {
            object@results$components
          })

setMethod("mmData", signature(object="MixModel"),
          function(object) {
            object@mmData
          })

setMethod("dim", signature(x="MixModel"),
          function(x){
            d <- length(mmData(x))
            c <- length(components(x))
            return(c("Number of data points" = d, "Number of components" = c))
          })

setMethod("length", signature(x="MixModel"),
          function(x){
            d <- length(mmData(x))
            names(d) <- "Number of data points"
            return(d)
          })

setMethod("listClassificationMethods", signature(object="MixModel"),
          function(object) {
            names(object@results$classification)
          })

setMethod("show", signature(object="MixModel"),
          function(object) {
            comps <- components(object)
            cat("\nMixModel object")
            cat("\n", "    ", "Number of data points:  ", length(mmData(object)), sep="")
            cat("\n", "    ", "Number of components:   ", length(comps), sep="")
            for (i in 1:length(comps)) {
              cat("\n", "        ", i, ": ", comps[[i]]@name, sep="")
              for (para in names(comps[[i]]@parameters)) {
                cat("\n", "             ", para, " = ", comps[[i]]@parameters[[para]], sep="")
              }
              cat("\n", "           ", "weight pi = ", weights(object)[i], sep="")
              cat("\n", "           ", "classified data points: ", sum(classification(object) == i), sep="")
            }
            cat("\n\n")
          })

setMethod("summary", signature(object="MixModel"),
          function(object) {
            compTypes <- sapply(components(object), function(x) {x@name})
            pi <- weights(object)
            result <- list()
            for (comp in unique(compTypes)) {
              ind <- which(compTypes == comp)
              parameterNames <- names(components(object)[ind][[1]]@parameters)
              df <- data.frame(t(sapply(components(object)[ind], function(x) {unlist(x@parameters)})))
              colnames(df) <- parameterNames
              rownames(df) <- as.character(ind)
              df$pi <- pi[ind]
              numPoints_temp <- table(classification(object))[as.character(ind)]
              df$numPoints <- numPoints_temp
              df$numPoints[is.na(numPoints_temp)] <- 0 
              result[[comp]] <- df
            }  
            return(result)
          })

setMethod("weights", signature(object="MixModel"),
          function(object) {
            object@results$pi
          })

setMethod("as.data.frame", signature(x="MixModel"),
          function(x, classificationMethod) {
              compNames <- sapply(components(x), function (c) {c@name})
              if (missing(classificationMethod)) {
                  classification <- factor(paste(classification(x), "_", compNames[classification(x)], sep=""),
                                           levels=paste(1:length(compNames), "_", compNames, sep=""))
              } else {
                  classification <- factor(paste(classification(x, method=classificationMethod), "_", compNames[classification(x, method=classificationMethod)], sep=""),
                                           levels=paste(1:length(compNames), "_", compNames, sep=""))
              }
              return(results <- data.frame(z=mmData(x), classification=classification, row.names=names(mmData(x))))
          })

# methods for class MixModelBayes
setMethod("chains", signature(object="MixModelBayes"),
          function(object) {
            object@chains
          })

setMethod("acceptanceRate", signature(object="MixModelBayes"),
          function(object) {
            mean(object@chains$dirichletParameterAcceptance)
          })

# methods for class MixModelML
setMethod("convergence", signature(object="MixModelML"),
          function(object) {
            object@convergence
          })
