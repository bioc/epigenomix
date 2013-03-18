setMethod("show", signature(object="MixtureComponent"),
          function(object) {
            cat("\nMixtureComponent object")
            cat("\n", "    ", "Name: ", object@name, sep="")
            cat("\n", "    ", "Parameters: ", sep="")
            for (i in 1:length(object@parameters)) {
              cat("\n", "      ", names(object@parameters)[i], ": ", object@parameters[[i]], sep="")
            }
            cat("\n", "    ", "Color: ", object@color, sep="")
            cat("\n\n")
          })
