# Functions that act directly on the DGE matrix, represented by a data.frame.

#' Remove cells by barcode
#'
setGeneric("removeCells",
           function(object, cells) {standardGeneric("removeCells")})
setMethod("removeCells",
          "data.frame",
          function(object, cells) {
            return (object[, !(names(object) %in% cells)])
          })

#' Filter out low quality cells
#'
#' Remove cells expressing less than a minimum of genes (default value is 2000).
#' @param object A \code{data.frame} representing DGE.
#' @param n The minimum number of genes required.
#' @return A \code{data.frame} without the low quality cells.
setGeneric("removeLowQualityCells",
           function(object, min.genes=1000, ...) {standardGeneric("removeLowQualityCells")})
setMethod("removeLowQualityCells",
          "data.frame",
          function(object, min.genes) {
            return (object[, names(object)[colSums(object != 0) >= min.genes]])
          })

setGeneric("keepBestCells",
           function(object, num.cells, min.num.trans=NULL) {standardGeneric("keepBestCells")})
setMethod("keepBestCells",
          "data.frame",
          function(object, num.cells, min.num.trans) {
            if (!is.null(min.num.trans)) {
              tdf <- object[, colSums(object) >= min.num.trans, drop=F]
              return (tdf[rowSums(tdf) > 0, ])
            }
            else {
              tdf <- object[, head(order(colSums(object), decreasing=T), num.cells), drop=F]
              return (tdf[rowSums(tdf) > 0, ])
            }
          })

#' Filter out low quality genes
#'
#' Remove genes which are expressed in less than a minimum number of cells.
#' @param object A \code{data.frame} representing DGE.
#' @param n The minimum number of cells.
#' @return A \code{data.frame} without the low quality genes.
setGeneric("removeLowQualityGenes",
           function(object, min.cells=3) {standardGeneric("removeLowQualityGenes")})
setMethod("removeLowQualityGenes",
          "data.frame",
          function(object, min.cells) {
            return (object[rownames(object) %in% rownames(object)[rowSums(object != 0) > min.cells], ])
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
            return (log(apply(object, 1, mean), 2))
          })

#' Compute the sum of transcript counts of each gene across all cells.
#'
#' @param object A \code{data.frame} representing DGE.
#' @return A vector of the logarithm of the sum of all UMI values for each gene.
#' A pseudocount is added to avoid infinities.
setGeneric("geneExpressionSumUMI",
           function (object) {standardGeneric("geneExpressionSumUMI")})
setMethod("geneExpressionSumUMI",
          "data.frame",
          function(object) {
            return (log(rowSums(object)+1, 2))
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

setGeneric("computeCorrelationSingleCellsVersusBulk",
           function(single.cells, bulk.data, method="pearson") {
             standardGeneric("computeCorrelationSingleCellsVersusBulk")
           })
setMethod("computeCorrelationSingleCellsVersusBulk",
          "data.frame",
          function(single.cells, bulk.data, method) {
            common.genes <- intersect(rownames(single.cells), rownames(bulk.data))
            single.cells <- single.cells[common.genes, ]
            single.cells <- max(colSums(single.cells)) * sweep(single.cells, MARGIN=2,
                                                               FUN='/', colSums(single.cells))
            single.cells <- log2(rowSums(single.cells) + 1)

            bulk <- bulk.data[common.genes, ]

            df = data.frame("sc"=single.cells, "bulk"=bulk)
            correlation = signif(cor(single.cells, bulk, method=method), 2)

            return (list(correlation, df))
          })

setGeneric("computeCorrelationCellToCellVersusBulk",
           function(single.cells, bulk.data, measure="umi", method="pearson") {
             standardGeneric("computeCorrelationCellToCellVersusBulk")
           })
setMethod("computeCorrelationCellToCellVersusBulk",
          "data.frame",
          function(single.cells, bulk.data, measure, method) {
            vectorOfCorrelations <- c()
            for (cell in names(single.cells)) {
              vectorOfCorrelations <- c(vectorOfCorrelations,
                                        computeCorrelationSingleCellsVersusBulk(single.cells[, cell, drop=F],
                                                                                bulk.data)[[1]])
            }
            return (vectorOfCorrelations)
          })

#' Compare average gene expression from single cell data against bulk data.
#'
#' @param single.cells A \code{data.frame} containing genes as rownames and columns as cells.
#' @param bulk.data A \code{data.frame} containing genes as rownames and a single column
#' with transcript counts.
#' @param log.space If True (default value) then the average of single cells data is
#' performed in log space.
#' @param ylab Label of the single cell sample.
#' @param col The color of the scatterplot points.
#' @return A scatterplot comparing the gene expression levels and with the correlation
#' among the samples computed.
setGeneric("compareSingleCellsAgainstBulk",
           function(single.cells, bulk.data, method="pearson",
                    ylab="single cell sample", col="steelblue", ...) {
             standardGeneric("compareSingleCellsAgainstBulk")
           })
setMethod("compareSingleCellsAgainstBulk",
          "data.frame",
          function(single.cells, bulk.data, method, ylab, col) {
            corr.df <- computeCorrelationSingleCellsVersusBulk(single.cells, bulk.data, method)

            corr = corr.df[[1]]
            df = corr.df[[2]]

            g <- (ggplot(df, aes(x=bulk, y=sc)) + geom_point(col=col, alpha=0.5, size = 3)
                  + theme_minimal() +  plotCommonTheme + xlab("log2(RPKM+1) [mRNA-seq]")
                  + ylab(paste0("log2(ATPM+1) [", ylab, "]"))
                  + ggtitle(paste0("R= ", corr)) + scale_x_continuous(expand=c(0.01, 0.02))
                  + theme(plot.title = element_text(size=22),
                          text=element_text(size=24, family="NumbusSan"),
                          plot.margin = unit(c(1, 1 , 0.5, 0.5), "cm"),
                          panel.border = element_rect(colour = "black", fill=NA, size=1),
                          panel.grid.major = element_blank()) )
            return (g)
          })

setGeneric("computeCellGeneFilteringFromBulk",
           function(single.cells, bulk.data, log.space=T, method="pearson",
                    min.cells=500, max.cells=2000, iteration.steps=50, do.plot=F) {
             standardGeneric("computeCellGeneFilteringFromBulk")
           })
setMethod("computeCellGeneFilteringFromBulk",
          "data.frame",
          function(single.cells, bulk.data, log.space, method,
                   min.cells, max.cells, iteration.steps, do.plot) {
            cor.list = c()
            cell.list = c()

            for (cells in seq(min.cells, max.cells, iteration.steps)) {
              temp.sample <- removeLowQualityGenes(removeLowQualityCells(single.cells, cells))
              temp.cor = computeCorrelationSingleCellsVersusBulk(temp.sample, bulk.data)[[1]]
              cor.list = c(cor.list, temp.cor)
              cell.list = c(cell.list, dim(temp.sample)[2])
              print (c(temp.cor, cells))
            }
            if (do.plot) {
              df <- data.frame("min.genes"=seq(min.cells, max.cells, iteration.steps),
                               "correlation"=cor.list,
                               "cells"=cell.list)
              g1 <- (ggplot(df, aes(x=min.genes, y=correlation)) + geom_point(size=3, col="steelblue")
                     + theme_minimal() + plotCommonGrid + plotCommonTheme)
              g2 <- (ggplot(df, aes(x=min.genes, y=cells)) + geom_point(size=3, col="steelblue")
                     + theme_minimal() + plotCommonGrid + plotCommonTheme)
              return(grid.arrange(g1, g2, ncol=2))
            }
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


