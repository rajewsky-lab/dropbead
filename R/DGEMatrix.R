DigitalGeneExpressionMatrix <- setClass(Class = "DigitalGeneExpressionMatrix",
                                        slots = c(dge = "data.frame")
)

# removes genes which are expressed in less than three cells
setGeneric(name = "removeLowQualityGenes",
           def = function(object) {standardGeneric("removeLowQualityGenes")})
setMethod(f = "removeLowQualityGenes",
          signature = "DigitalGeneExpressionMatrix",
          function(object) {
            return (object@dge[!rownames(object@dge) %in% names(rowSums(object@dge)[rowSums(object@dge) < 3]), ])
          })

setGeneric(name = "removeLowQualityCells",
           def = function(object) {standardGeneric("removeLowQualityCells")})
setMethod(f = "removeLowQualityCells",
          signature = "DigitalGeneExpressionMatrix",
          function(object) {
            return (0)
          })

setGeneric(name = "computeGenesVariability",
           def = function(object) {standardGeneric("computeGenesVariability")})
setMethod(f = "computeGenesVariability",
          signature = "DigitalGeneExpressionMatrix",
          function(object) {
            return (0)
          })

