setClass("ChIPseqSet",
         contains = "eSet"
         )

setClass("MixModel",
         representation = representation(
           mmData = "numeric",
           configuration = "list",
           results = "list",
           "VIRTUAL")
         )

setClass("MixModelBayes",
         representation = representation(
           chains = "list"),
         contains = "MixModel"
         )

setClass("MixModelML",
         representation = representation(
           convergence = "list"),
         contains = "MixModel"
         )

setClass("MixtureComponent",
         representation = representation(
           name = "character",
           parameters = "list",
           pdf = "function",
           color = "character")
         )
