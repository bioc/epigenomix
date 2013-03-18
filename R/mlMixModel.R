.mlMixModel <- function(z, normNull=c(), expNeg=c(), expPos=c(),
                        sdNormNullInit=c(), rateExpNegInit=c(), rateExpPosInit=c(),
                        piInit=c(), maxIter=500, tol=0.001) {

  # total number of components
  k <- length(c(normNull, expNeg, expPos))  

  # validate arguments
  if (k < 1) {
    stop("At least one component must be specified.")
  } 
  if (!all(is.element(1:k, c(normNull, expNeg, expPos)))) {
    stop(paste("Components' indices must be in 1, ..., ", k, ".", sep=""))
  }
  if (any(duplicated(c(normNull, expNeg, expPos)))) {
    stop(paste("Components' indices must be unique."))
  }
  if (length(piInit) != k) {
    stop("Length of argument pi must equal the number of mixture components.")
  }
  if (!(length(normNull) == length(sdNormNullInit) &
        length(expNeg) == length(rateExpNegInit) &
        length(expPos) == length(rateExpPosInit))) {
    stop("Number of given components must fit to the number of given initialization values.")
  }

  # standardize weights
  piInit <- piInit / sum(piInit)
  pi <- piInit
    
  # construct lists with components parameters
  comps = list()
  for (ind in normNull) {
    comps[[ind]] = new("MixtureComponent",
           name="NormNull",
           parameters=list(mean=0, sd=sdNormNullInit[ind==normNull]),
           pdf=dnorm,
           color="blue")
  }
  for (ind in expNeg) {
    comps[[ind]] = new("MixtureComponent",
           name="ExpNeg",
           parameters=list(rate=rateExpNegInit[ind==expNeg]),
           pdf=function(x, rate) { dexp(-x, rate) },
           color="red")
  }
  for (ind in expPos) {
    comps[[ind]] = new("MixtureComponent",
           name="ExpPos",
           parameters=list(rate=rateExpPosInit[ind==expPos]),
           pdf=dexp,
           color="green")
  }
  compsInit <- comps

  # start EM algorithm
  iter <- 1
  logLik <- -Inf
  deltaLogLik <- Inf
  while (iter < maxIter & deltaLogLik > tol) {

    # E step
    f <- matrix(nrow=length(z), ncol=k)
    for (ind in normNull) {
      f[, ind] <- dnorm(z, mean=comps[[ind]]@parameters$mean, sd=comps[[ind]]@parameters$sd)
    }
    for (ind in expNeg) {
      f[, ind] <- dexp(-z, comps[[ind]]@parameters$rate)
    }
    for (ind in expPos) {
      f[, ind] <- dexp(z, comps[[ind]]@parameters$rate)
    }
    pi.f <- t(pi * t(f))
    c <- pi.f / apply(pi.f, 1, sum)
    pi <- apply(c, 2, sum) / nrow(c)
  
    # M step
    w <- apply(c, 2, function(x) {x / sum(x)})
    for (ind in normNull) {
      comps[[ind]]@parameters$sd <- sqrt(sum(w[, ind] * z^2)) # update sigma of normal comps
    }
    for (ind in expNeg) {
      comps[[ind]]@parameters$rate <- 1 / sum(w[, ind] * -z)   # update lambda of neg. exp. comps
    }
    for (ind in expPos) {
      comps[[ind]]@parameters$rate <- 1 / sum(w[, ind] * z)    # update lambda of pos. exp. comps
    }

    iter <- iter + 1
    oldLogLik <- logLik
    logLik <- sum(log(apply(pi.f, 1, sum)))
    deltaLogLik <- logLik - oldLogLik

    print(paste("Iteration: ", iter, "  deltaLogLik: ", deltaLogLik,
        "  logLik:", logLik, sep=""))
  }

  classification <- apply(c, 1, which.max)

  mm <- new("MixModelML",
            mmData=z,
            configuration=list(inits=list(components=compsInit,pi=piInit),
              convergence=list(maxIter=maxIter, tol=tol)),
            results=list(components=comps, pi=pi, classification=list(maxDens=classification)),
            convergence=list(iterations=iter, deltaLogLik=deltaLogLik, logLik=logLik))
  
  return(mm)
}

setMethod("mlMixModel", signature=c(z="numeric"), .mlMixModel)
