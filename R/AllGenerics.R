# methods-ChIPseqSet.R
setGeneric("chipVals", function(object)
           standardGeneric("chipVals"))

setGeneric("chipVals<-", function(object, value)
           standardGeneric("chipVals<-"))

# methods-MixModel.R
setGeneric("chains", function(object)
           standardGeneric("chains"))

setGeneric("classification", function(object, method)
           standardGeneric("classification"))

setGeneric("components", function(object)
           standardGeneric("components"))

setGeneric("convergence", function(object)
           standardGeneric("convergence"))

setGeneric("mmData", function(object)
           standardGeneric("mmData"))

setGeneric("listClassificationMethods", function(object)
           standardGeneric("listClassificationMethods"))

setGeneric("weights", function(object)
           standardGeneric("weights"))

# other methods
setGeneric("bayesMixModel", function(z, normNull=c(), expNeg=c(), expPos=c(), gamNeg=c(), gamPos=c(),
                                     sdNormNullInit=c(), rateExpNegInit=c(), rateExpPosInit=c(),
                                     shapeGamNegInit=c(), scaleGamNegInit=c(), shapeGamPosInit=c(), scaleGamPosInit=c(),
                                     piInit, classificationsInit, dirichletParInit=1, shapeDir=1, scaleDir=1,
                                     shapeNorm0=c(), scaleNorm0=c(), shapeExpNeg0=c(), scaleExpNeg0=c(), shapeExpPos0=c(), scaleExpPos0=c(),
                                     shapeGamNegAlpha0=c(), shapeGamNegBeta0=c(), scaleGamNegAlpha0=c(), scaleGamNegBeta0=c(),
                                     shapeGamPosAlpha0=c(), shapeGamPosBeta0=c(), scaleGamPosAlpha0=c(), scaleGamPosBeta0=c(),
                                     itb, nmc, thin, average="mean", gammaProposalFactor, gammaShapeGrid)
           standardGeneric("bayesMixModel"),
           signature="z")

setGeneric("integrateData", function(expr, chipseq, factor, reference)
           standardGeneric("integrateData"))

setGeneric("matchProbeToPromoter", function(probeToTranscript, transcriptToTSS,
                                            promWidth, mode)
           standardGeneric("matchProbeToPromoter"))

setGeneric("mlMixModel", function(z, normNull=c(), expNeg=c(), expPos=c(),
                                  sdNormNullInit=c(), rateExpNegInit=c(), rateExpPosInit=c(),
                                  piInit=c(), maxIter=500, tol=0.001)
           standardGeneric("mlMixModel"),
           signature="z")

setGeneric("normalizeChIP", function(object, method)
           standardGeneric("normalizeChIP"))

setGeneric("plotChains", function(object, chain, component, itb=1, thin=1, cols, ...)
           standardGeneric("plotChains"))

setGeneric("plotClassification", function(object, method, ...)
           standardGeneric("plotClassification"),
           signature="object")

setGeneric("plotComponents", function(object, density=FALSE, ...)
           standardGeneric("plotComponents"),
           signature="object")

setGeneric("summarizeReads", function(object, regions, summarize)
           standardGeneric("summarizeReads"))
