# Functions that act directly on the DGE matrix, represented by a data.frame.

#' Filter out low quality genes
#'
#' Remove genes which are expressed in less than a minimum number of cells.
#' @param object A \code{data.frame} representing DGE.
#' @param n The minimum number of cells.
#' @return A \code{data.frame} without the low quality genes.
setGeneric("removeLowQualityGenes",
           function(object, n=3) {standardGeneric("removeLowQualityGenes")})
setMethod("removeLowQualityGenes",
          "data.frame",
          function(object, n) {
            return (object[rownames(object) %in% rownames(object)[rowSums(object != 0) > n], ])
          })

#' Filter out low quality cells
#'
#' Remove cells expressing less than a minimum of genes (default value is 2000).
#' @param object A \code{data.frame} representing DGE.
#' @param n The minimum number of genes required.
#' @return A \code{data.frame} without the low quality cells.
setGeneric("removeLowQualityCells",
           function(object, n=2000, ...) {standardGeneric("removeLowQualityCells")})
setMethod("removeLowQualityCells",
          "data.frame",
          function(object, n) {
            return (object[, names(object)[colSums(object != 0) >= n]])
          })

#' Compute the average expression of each gene across all cells.
#'
#' @param object A \code{data.frame} representing DGE.
#' @return A vector of mean expression values for each gene. Computation in log space.
setGeneric("geneExpressionMean",
           function (object) {standardGeneric("geneExpressionMean")})
setMethod("geneExpressionMean",
          "data.frame",
          function(object) {
            return (log(apply(object, 1, mean)+1, 2))
          })

#' Compute the dispersion of each gene across all cells.
#'
#' @param object A \code{data.frame} representing DGE.
#' @return A vector of variance/mean expression values for each gene. Computation in log space.
setGeneric("geneExpressionDispersion",
           function (object) {standardGeneric("geneExpressionDispersion")})
setMethod("geneExpressionDispersion",
          "data.frame",
          function(object) {
            return (log(apply(object, 1, var)/apply(object, 1, mean), 2))
          })

#' Gene variability
#'
#' Compute the variability of each gene across all cells. The average expression
#' of each gene is first computed, sorted and the genes are placed into \code{b}
#' bins. The dispersion of genes into each bin is then z-normalized. A cuttoff
#' is then used to identify the top \code{h} highly variable genes. Genes with 0
#' mean expression should be removed before using this function. Computation
#' is performed in log space.
#'
#' @param object A \code{data.frame} representing DGE.
#' @param bins The number of bins.
#' @return A vector of the top highly varied genes.
setGeneric("geneExpressionVariability",
           function (object, bins=20, low=F, do.plot=F) {standardGeneric("geneExpressionVariability")})
setMethod("geneExpressionVariability",
          "data.frame",
          function(object, bins, low, do.plot) {
            x = geneExpressionMean(object)
            y = geneExpressionDispersion(object)
            x.bins = cut(x, bins)
            names(x.bins) = names(x)
            y.mean = tapply(y, x.bins, mean)
            y.sd = tapply(y, x.bins, sd)
            ret = (y - y.mean[as.numeric(x.bins)])/y.sd[as.numeric(x.bins)]
            names(ret) = names(x)
            high.var <- intersect(names(ret[ret > 2/log(2) & ret < 12]), names(x[x > 2]))
            low.var <- c()

            if (do.plot) {
              plot(x, ret, pch=19, cex=0.3)
              points(x[high.var], ret[high.var], pch=19, cex=0.5, col = 'red')
            }
            if (low) {
              return (low.var)
            }
            return (high.var)
          })

#' Identify doublets
#'
#' Subsets the DGE matrix to genes which are highly expressed and at the
#' same time lowly varied across all cells. Comparing the cells then probabilities
#' are assigned for each cell to be a doublet of the same species.
#' @param object A \code{data.frame} representing DGE.0
#' @return A \code{vector} containing the assigned probabilities for each cell.
setGeneric("identifyDoublets",
           function(object) {standardGeneric("identifyDoublets")})
setMethod("identifyDoublets",
          "data.frame",
          function(object) {
            lv.genes <- geneExpressionVariability(object, low=TRUE)
            c.genes <- intersect(names(tail(sort(geneExpressionMean(object)), length(lv.genes))), names(lv.genes))
            dge.red <- object[c.genes, ]

            return (dge.red)
          })


