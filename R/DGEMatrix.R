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
            return (apply(log(object+1, 2), 1, mean))
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
            return (apply(log(object+1, 2), 1, var)/apply(log(object+1, 2), 1, mean))
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
#' @param low Boolean. If True then the lowly varied genes are returned.
#' @return A vector of the top \code{ntop} highly varied genes.
setGeneric("geneExpressionVariability",
           function (object, bins=20, low=FALSE, all=FALSE) {standardGeneric("geneExpressionVariability")})
setMethod("geneExpressionVariability",
          "data.frame",
          function(object, bins, low, all) {
            mean.sorted <- sort(geneExpressionMean(object))
            bin.size = ceiling(dim(object)[1]/bins)
            all.genes <- c()
            hv.genes <- c()
            lv.genes <- c()

            for (bin in 0:(bins-1)) {
              if (bin == bins-1) {
                endpoint = dim(object)[1]
              } else {
                endpoint = (bin+1)*bin.size
              }
              bin.genes <- mean.sorted[(bin*bin.size+1):endpoint]
              bin.disp <- geneExpressionDispersion(object)[names(bin.genes)]
              all.genes <- c(all.genes, scale(bin.disp)[, ])
              hv.genes <- c(hv.genes, scale(bin.disp)[scale(bin.disp) > 2, ])
              lv.genes <- c(lv.genes, scale(bin.disp)[scale(bin.disp) < -1.6, ])
            }

            if (all == T) {return (all.genes)}
            if (low == F) {return (hv.genes)}
            if (low == T) {return (lv.genes)}
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


