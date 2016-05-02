SingleSpeciesSample <- setClass(Class = "SingleSpeciesSample",
                                slots = c(species1 = "character",
                                          cells = "vector",
                                          genes = "vector",
                                          dge = "data.frame"))

setMethod("initialize",
          "SingleSpeciesSample",
          function (.Object, species1 = "", cells = c(), genes = c(), dge = data.frame()) {
            .Object@species1 = species1
            .Object@cells = names(dge)
            .Object@genes = rownames(dge)
            .Object@dge = dge
            .Object
          })

#' Compute genes per cell
#'
#' @param object A Single species sample.
#' @return A \code{data.frame} with cells, gene counts and species.
setGeneric("computeGenesPerCell",
           function(object, min.reads=2, ...) {standardGeneric("computeGenesPerCell")})
setMethod("computeGenesPerCell",
          "SingleSpeciesSample",
          function(object, min.reads) {
            genes <- data.frame("cells" = names(colSums(object@dge >= min.reads)),
                                "counts" = as.numeric(colSums(object@dge >= min.reads)),
                                "species" = object@species1)
            return (genes)
          })

#' Compute transcripts per cell
#'
#' @param object A Single species sample.
#' @return A \code{data.frame} with cells, transcript counts and species.
setGeneric("computeTranscriptsPerCell",
           function(object, ...) {standardGeneric("computeTranscriptsPerCell")})
setMethod("computeTranscriptsPerCell",
          "SingleSpeciesSample",
          function(object) {
            transcripts <- data.frame("cells" = object@cells,
                                      "counts" = as.numeric(colSums(object@dge)),
                                      "species" = object@species1)
            return (transcripts)
          })

setMethod("removeLowQualityCells",
          "SingleSpeciesSample",
          function(object, min.genes) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityCells(object@dge, min.genes)))
          })

setMethod("keepBestCells",
          "SingleSpeciesSample",
          function(object, num.cells, min.num.trans) {
            return (new("SingleSpeciesSample", species1=object@species1,
                        dge=keepBestCells(object@dge, num.cells, min.num.trans)))
          })

setMethod("removeLowQualityGenes",
          "SingleSpeciesSample",
          function(object, min.cells) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityGenes(object@dge, min.cells)))
          })

setMethod("geneExpressionMean",
          "SingleSpeciesSample",
          function(object) {
            return (geneExpressionMean(object@dge))
          })

setMethod("geneExpressionDispersion",
          "SingleSpeciesSample",
          function(object) {
            return (geneExpressionDispersion(object@dge))
          })

setMethod("geneExpressionVariability",
          "SingleSpeciesSample",
          function(object, bins, low, do.plot) {
            return (geneExpressionVariability(object@dge, bins, low, do.plot))
          })

#' List cells that are candidates for collapsing.
#'
#' Identify and list cells which share 11 bases in their barcodes and only the last
#' one is different. The cells are marked as candidates if and only if they're classified
#' as belonging to the same species.
#' @param A \code{SingleSpeciesSample} object.
#' @return A list of pairs od candidate cells marked for collapsing.
setGeneric("listCellsToCollapse",
           function(object, ...) {
             standardGeneric("listCellsToCollapse")
           })
setMethod("listCellsToCollapse",
          "SingleSpeciesSample",
          function (object) {
            theListOfCellPairs <- list()
            cells <- sort(object@cells)
            for (cell in 1:(length(cells)-1)) {
              if (substr(cells[cell], 1, 11) == substr(cells[cell+1], 1, 11)) {
                if (substr(cells[cell], nchar(cells[cell]), nchar(cells[cell])) == "N" |
                    substr(cells[cell+1], nchar(cells[cell+1]), nchar(cells[cell+1])) == "N") {
                  theListOfCellPairs <- c(theListOfCellPairs, list(c(cells[cell], cells[cell+1])))
                }
              }
            }
            return(theListOfCellPairs)
          })

#' Collapse cells by barcodes similarity
#'
#' Collapse cells which have barcodes differing by one mutation on the last base.
setGeneric("collapseCellsByBarcode",
           function(object, ...) {
             standardGeneric("collapseCellsByBarcode")
           })
setMethod("collapseCellsByBarcode",
          "SingleSpeciesSample",
          function(object) {
            listOfCells <- listCellsToCollapse(object)
            if (length(listOfCells) == 0) {
              return(object)
            }

            for (index in 1:length(listOfCells)) {
              object@dge <- cbind(object@dge, rowSums(object@dge[, listOfCells[[index]]]))
            }
            object@dge <- object@dge[, !names(object@dge) %in% unlist(listOfCells)]

            names(object@dge)[(length(names(object@dge)) -
                                     length(listOfCells) + 1):length(names(object@dge))] <- unlist(listOfCells)[seq(1, length(unlist(listOfCells)), 2)]
            object@cells <- names(object@dge)
            return (object)
          })

#' Compare gene expression levels between two single species samples
#'
#' @param object1 A \code{SingleSpeciesSample} object.
#' @param object2 A \code{SingleSpeciesSample} object.
#' @param name1 The name of the first sample.
#' @param name2 The name of the second sample.
#' @return A plot comparing the gene expression levels.
setGeneric("compareGeneExpressionLevels",
           function(object1, object2, name1="sample1", name2="sample2",
                    col="steelblue", method="pearson", ...) {
             standardGeneric("compareGeneExpressionLevels")
           })
setMethod("compareGeneExpressionLevels",
          "SingleSpeciesSample",
          function(object1, object2, name1, name2, col, method) {
            common.genes <- intersect(object1@genes, object2@genes)

            object1 <- object1@dge[common.genes, ]
            object2 <- object2@dge[common.genes, ]

            big.df <- data.frame("genes" = common.genes,
                                 "sample1" = log(as.numeric(rowSums(object1))+1, 2),
                                 "sample2" = log(as.numeric(rowSums(object2))+1, 2))

            cor <- signif(cor(big.df$sample1, big.df$sample2, method=method), 2)

            comp.plot <- (ggplot(data = big.df, aes(x = sample1, y = sample2))
                          + ggtitle(paste0("R= ", cor))
                          + xlab(name1) + ylab(name2)
                          + geom_point(col=col, alpha=0.5, size=2)
                          + scale_y_continuous(expand = c(0.02, 0.02)) + scale_x_continuous(expand = c(0.02, 0.02))
                          + theme_minimal() + plotCommonGrid + plotCommonTheme
                          + theme(axis.ticks.y = element_blank()))

            return (comp.plot)
          })

setMethod("computeCorrelationSingleCellsVersusBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, measure, method) {
            corr = computeCorrelationSingleCellsVersusBulk(single.cells@dge, bulk.data, measure, method)[[1]]
            df = computeCorrelationSingleCellsVersusBulk(single.cells@dge, bulk.data, measure, method)[[2]]
            return (list(corr, df))
          })

setMethod("computeCorrelationCellToCellVersusBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, measure, method) {
            return (computeCorrelationCellToCellVersusBulk(single.cells@dge, bulk.data, measure, method))
          })

setMethod("compareSingleCellsAgainstBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, measure, method, ylab, col) {
            return (compareSingleCellsAgainstBulk(single.cells@dge, bulk.data,
                                                  measure, method, ylab, col))
            })

setMethod("computeCellGeneFilteringFromBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, log.space, method,
                   min.cells, max.cells, iteration.steps) {
            computeCellGeneFilteringFromBulk(single.cells@dge, bulk.data, log.space, method,
                                             min.cells, max.cells, iteration.steps)
            })

setGeneric("assignCellCyclePhases",
           function(object, gene.file="~/Desktop/things/git/dropseq/data/cell_cycle_genes.xlsx", do.plot=T, ...) {
             standardGeneric("assignCellCyclePhases")
           })
setMethod("assignCellCyclePhases",
          "SingleSpeciesSample",
          function(object, gene.file, do.plot) {
            require(xlsx)
            file.ext = gsub(".*\\.(.*)$","\\1", gene.file)
            
			if (file.ext == "xlsx") {
				if (object@species1 == "human") {
					cc_genes <- read.xlsx(gene.file, sheetIndex = 2, stringsAsFactors = F)
				}

				if (object@species1 == "mouse") {
					cc_genes <- read.xlsx(gene.file, sheetIndex = 3, stringsAsFactors = F)
				}

			} else if (file.ext == "csv") {
			    cc_genes = read.table(gene.file, header=T, stringsAsFactors=F, sep=",")
		    }  

            g.g1s <- gsub(" ", "", cc_genes$G1.S)
            g.g1s <- g.g1s[!is.na(g.g1s)]
            g.s <- gsub(" ", "", cc_genes$S)
            g.s <- g.s[!is.na(g.s)]
            g.g2m <- gsub(" ", "", cc_genes$G2.M)
            g.g2m <- g.g2m[!is.na(g.g2m)]
            g.m <- gsub(" ", "", cc_genes$M)
            g.m <- g.m[!is.na(g.m)]
            g.mg1 <- gsub(" ", "", cc_genes$M.G1)
            g.mg1 <- g.mg1[!is.na(g.mg1)]

            phase_score = data.frame("G1.S"=rep(0, length(object@cells)),
                                     "S"=rep(0, length(object@cells)),
                                     "G2.M"=rep(0, length(object@cells)),
                                     "M"=rep(0, length(object@cells)),
                                     "M.G1"=rep(0, length(object@cells)))

            for (i in 1:length(object@cells)) {
              phase_score$G1.S[i] = mean(log(object@dge[g.g1s, i]+1, 2), na.rm=T)
              phase_score$S[i] = mean(log(object@dge[g.s, i]+1, 2), na.rm=T)
              phase_score$G2.M[i] = mean(log(object@dge[g.g2m, i]+1, 2), na.rm=T)
              phase_score$M[i] = mean(log(object@dge[g.m, i]+1, 2), na.rm=T)
              phase_score$M.G1[i] = mean(log(object@dge[g.mg1, i]+1, 2), na.rm=T)
            }

            phase_score_norm <- data.frame(apply(phase_score, 2, scale))
            phase_score_norm2 <- data.frame(t(apply(phase_score_norm, 1, scale)))
            names(phase_score_norm2) <- names(phase_score_norm)
            phases <- apply(phase_score_norm2, 1, which.max)

            df1 <- phase_score_norm2[phases == 1, ]
            df1 <- round(df1[order(df1$G1.S, decreasing = T), ], 2)
            df1 <- df1[order(df1$G1.S, df1$S, df1$G2.M, df1$M, df1$M.G1, decreasing = T), ]

            df2 <- phase_score_norm2[phases == 2, ]
            df2 <- round(df2[order(df2$S, decreasing = T), ], 2)
            df2 <- df2[order(df2$S, df2$G1.S, df2$G2.M, df2$M, df2$M.G1, decreasing = T), ]

            df3 <- phase_score_norm2[phases == 3, ]
            df3 <- round(df3[order(df3$G2.M, decreasing = T), ], 2)
            df3 <- df3[order(df3$G2.M,  df3$M, df3$S, df3$M.G1, df3$G1.S, decreasing = T), ]

            df4 <- phase_score_norm2[phases == 4, ]
            df4 <- round(df4[order(df4$M, decreasing = T), ], 2)
            df4 <- df4[order(df4$M, df4$M.G1, df4$G2.M, df4$S, df4$G1.S, decreasing = T), ]

            df5 <- phase_score_norm2[phases == 5, ]
            df5 <- round(df5[order(df5$M.G1, decreasing = T), ], 2)
            df5 <- df5[order(df5$M.G1, df5$G1.S, df5$G2.M, df5$S, df5$M,  decreasing = T), ]

            heatmap_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#cc4140"))
            if (do.plot) {
              heatmap.2(as.matrix(t(rbind(df1, df2, df3, df4, df5))), trace='none', Rowv=F, Colv=F, dendrogram='none', labCol=F, col=heatmap_palette(10))
            }

            return(rbind(df1, df2, df3, df4, df5))
          })
