DigitalGeneExpressionMatrix <- setClass(Class = "DigitalGeneExpressionMatrix",
                                        slots = c(dge = "data.frame")
)

# removes genes which are expressed in less than three cells
setGeneric(name = "removeLowQualityGenes",
           def = function(object) {standardGeneric("removeLowQualityGenes")})
setMethod(f = "removeLowQualityGenes",
          signature = "DigitalGeneExpressionMatrix",
          function(object) {
            return (new("DigitalGeneExpressionMatrix",
                        dge = object@dge[rownames(object@dge) %in% rownames(object@dge)[rowSums(object@dge != 0) > 3], ]))
          })

# remove cells which express less than 2000 genes
setGeneric(name = "removeLowQualityCells",
           def = function(object, ...) {standardGeneric("removeLowQualityCells")})
setMethod(f = "removeLowQualityCells",
          signature = "DigitalGeneExpressionMatrix",
          function(object) {
            return (new("DigitalGeneExpressionMatrix",
                        dge = object@dge[, names(object@dge)[colSums(object@dge != 0) > 2000]]))
          })

setGeneric(name = "computeGenesVariability",
           def = function(object, ...) {standardGeneric("computeGenesVariability")})
setMethod(f = "computeGenesVariability",
          signature = "DigitalGeneExpressionMatrix",
          function(object) {
            return (0)
          })

