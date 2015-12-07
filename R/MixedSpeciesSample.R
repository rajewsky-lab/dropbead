MixedSpeciesSample <- setClass(Class = "MixedSpeciesSample",
                         slots = c(species2 = "character"),
                         contains = "SingleSpeciesSample"
                         )

setGeneric(name = "splitDgeByGenesOfSpecies",
           def = function(object) {standardGeneric("splitDgeByGenesOfSpecies")})
setMethod(f = "splitDgeByGenesOfSpecies",
          signature = "DigitalGeneExpressionMatrix",
          function(object) {
            object.species2 <- object@dge[grep("[a-z]", rownames(object@dge)), ]
            object.species1 <- object@dge[setdiff(1:(dim(object@dge)[1]), grep("[a-z]", rownames(object@dge))), ]
            return (list(object.species1, object.species2))
          })

#' Classify cells to species
#'
#'
#' Classify the cells to species according to the transcripts they express.
#' The differentiation of species is performed internally and according to
#' the \code{species1} and \code{species2} labels of the sample. If the
#' classificiation is not confident enough, the cell is characterized as doublet.
#' @param object A \code{MixedSpeciesSample} object.
#' @param threshold The threshold which the ratio of transcripts of one species
#' over the other has to surpass in order to succesfully assign a cell to a species.
#' Below that threshold the cell is characterized as a doublet.
#' @return A \code{data.frame} containing the cell barcodes, the number of transcripts
#' per species and the characterization of the cell.
setGeneric(name = "classifyCellsAndDoublets",
           def = function(object, threshold=0.9) {standardGeneric("classifyCellsAndDoublets")})
setMethod(f = "classifyCellsAndDoublets",
          signature = "MixedSpeciesSample",
          function(object, threshold) {
            object.species1 <- splitDgeByGenesOfSpecies(object@dge)[[1]]
            object.species2 <- splitDgeByGenesOfSpecies(object@dge)[[2]]

            df <- data.frame("cell" = names(object@dge@dge),
                             "s1" = as.vector(colSums(object.species1)),
                             "s2" = as.vector(colSums(object.species2)),
                             "species" = "", stringsAsFactors = F)
            for (i in 1:dim(df)[1]) {
               if (df$s1[i] < quantile(df$s1, 0.25) & df$s2[i] < quantile(df$s2, 0.25)) {
                 df$species[i] = "undefined"
                 next
               }
              if (df$s1[i] > df$s2[i] & df$s1[i]/(df$s1[i] + df$s2[i]) > threshold) {
                df$species[i] = object@species1
                next
              }
              if (df$s2[i] > df$s1[i] & df$s2[i]/(df$s1[i] + df$s2[i]) > threshold) {
                df$species[i] = object@species2
                next
                }
              df$species[i] = "mixed"
            }
            names(df) <- c("cell", object@species1, object@species2, "species")

            return (df)
          })

setGeneric(name = "splitDgeByGenesAndCellsOfSpecies",
           def = function(object, threshold) {standardGeneric("splitDgeByGenesAndCellsOfSpecies")})
setMethod(f = "splitDgeByGenesAndCellsOfSpecies",
          signature = "MixedSpeciesSample",
          function(object, threshold) {
            object.species1 <- splitDgeByGenesOfSpecies(object@dge)[[1]]
            object.species2 <- splitDgeByGenesOfSpecies(object@dge)[[2]]

            df <- classifyCellsAndDoublets(object, threshold)

            object.species1 <- object.species1[, df[df$species == object@species1, ]$cell]
            object.species2 <- object.species2[, df[df$species == object@species2, ]$cell]

            return (list(object.species1, object.species2))
          })

setGeneric(name = "splitMixedSpeciesSampleToSingleSpecies",
           def = function(object, threshold) {standardGeneric("splitMixedSpeciesSampleToSingleSpecies")})
setMethod(f = "splitMixedSpeciesSampleToSingleSpecies",
          signature = "MixedSpeciesSample",
          function(object, threshold) {
            s1 <- new("SingleSpeciesSample",
                      species1 = object@species1,
                      dge = new("DigitalGeneExpressionMatrix",
                                dge = splitDgeByGenesAndCellsOfSpecies(object, threshold)[[1]]))
            s2 <- new("SingleSpeciesSample",
                      species1 = object@species2,
                      dge = new("DigitalGeneExpressionMatrix",
                                dge = splitDgeByGenesAndCellsOfSpecies(object, threshold)[[2]]))
            return (list(s1, s2))
          })

setMethod(f = "computeGenesPerCell",
          signature = "MixedSpeciesSample",
          function(object, threshold) {
            return(rbind.fill(lapply(splitMixedSpeciesSampleToSingleSpecies(object, threshold), computeGenesPerCell)))
          })

setMethod(f = "computeTranscriptsPerCell",
          signature = "MixedSpeciesSample",
          function(object, threshold) {
            return (rbind.fill(lapply(splitMixedSpeciesSampleToSingleSpecies(object, threshold), computeTranscriptsPerCell)))
          })

setGeneric(name = "listCellsToCollapse",
           def = function(object, threshold) {
             standardGeneric("listCellsToCollapse")
           })
setMethod(f = "listCellsToCollapse",
          signature = "MixedSpeciesSample",
          function (object, threshold) {
            theListOfCellPairs <- list()
            first.df <- classifyCellsAndDoublets(object, threshold)
            first.df <- first.df[order(first.df$cell), ]

            cells <- first.df$cell
            for (cell in 1:(length(cells)-1)) {
              if (substr(cells[cell], 1, 11) == substr(cells[cell+1], 1, 11)) {
                if (substr(cells[cell], nchar(cells[cell]), nchar(cells[cell])) == "N" |
                    substr(cells[cell+1], nchar(cells[cell+1]), nchar(cells[cell+1])) == "N") {
                  if (first.df$species[cell] == first.df$species[cell+1]) {
                    theListOfCellPairs <- c(theListOfCellPairs, list(c(cells[cell], cells[cell+1])))
                  }
                }
              }
            }
            return(theListOfCellPairs)
          })

setGeneric(name = "collapseCellsByBarcode",
           def = function(object, threshold) {
             standardGeneric("collapseCellsByBarcode")
           })
setMethod(f = "collapseCellsByBarcode",
          signature = "MixedSpeciesSample",
          function(object, threshold) {
            listOfCells <- listCellsToCollapse(object, threshold)

            for (index in 1:length(listOfCells)) {
              object@dge@dge <- cbind(object@dge@dge, rowSums(object@dge@dge[, listOfCells[[index]]]))
            }
            object@dge@dge <- object@dge@dge[, !names(object@dge@dge) %in% unlist(listOfCells)]

            names(object@dge@dge)[(length(names(object@dge@dge)) -
                                     length(listOfCells) + 1):length(names(object@dge@dge))] <- unlist(listOfCells)[seq(1, length(unlist(listOfCells)), 2)]
            return (object)
          })

setGeneric(name = "compareGeneExpressionLevels",
           def = function(object1, threshold1, object2, threshold2) {
             standardGeneric("compareGeneExpressionLevels")
             })
setMethod(f = "compareGeneExpressionLevels",
          signature = "MixedSpeciesSample",
          function(object1, threshold1, object2, threshold2) {
            object1.species1 <- splitDgeByGenesAndCellsOfSpecies(object1, threshold1)[[1]]
            object1.species2 <- splitDgeByGenesAndCellsOfSpecies(object1, threshold1)[[2]]
            object2.species1 <- splitDgeByGenesAndCellsOfSpecies(object2, threshold2)[[1]]
            object2.species2 <- splitDgeByGenesAndCellsOfSpecies(object2, threshold2)[[2]]

            common.genes.species1 <- intersect(row.names(object1.species1), row.names(object2.species1))
            common.genes.species2 <- intersect(row.names(object1.species2), row.names(object2.species2))

            object1.species1 <- object1.species1[row.names(object1.species1) %in% common.genes.species1, ]
            object1.species2 <- object1.species2[row.names(object1.species2) %in% common.genes.species2, ]
            object2.species1 <- object2.species1[row.names(object2.species1) %in% common.genes.species1, ]
            object2.species2 <- object2.species2[row.names(object2.species2) %in% common.genes.species2, ]

            big.df <- data.frame("genes" = c(common.genes.species1, common.genes.species2),
                                 "sample1" = c(log(as.numeric(rowSums(object1.species1)), 2),
                                               log(as.numeric(rowSums(object1.species2)), 2)),
                                 "sample2" = c(log(as.numeric(rowSums(object2.species1)), 2),
                                               log(as.numeric(rowSums(object2.species2)), 2)),
                                 "species" = c(rep(object1@species1, length(common.genes.species1)),
                                               rep(object1@species2, length(common.genes.species2))))
            big.df <- big.df[big.df$sample1 >= 0, ]
            big.df <- big.df[big.df$sample2 >= 0, ]

            cor.species1 <- signif(cor(big.df$sample1[big.df$species == object1@species1],
                                       big.df$sample2[big.df$species == object1@species1],
                                       method="spearman"), 2)
            cor.species2 <- signif(cor(big.df$sample1[big.df$species == object1@species2],
                                       big.df$sample2[big.df$species == object1@species2],
                                       method="spearman"), 2)

            levels(big.df$species) <- c(paste0(object1@species1, " (R=", cor.species1, ")"),
                                        paste0(object1@species2, " (R=", cor.species2, ")"))

            comp.plot <- (ggplot(data = big.df, aes(x = sample1, y = sample2))
                          + xlab(expression(paste(log[2]~'transcripts (sample1)')))
                          + ylab(expression(paste(log[2]~'transcripts (sample2)')))
                          + geom_point(aes(col = species), alpha = 0.5, size = 2)
                          + facet_grid(~ species, labeller = ) + guides(col = F)
                          + scale_y_continuous(expand = c(0, 0))
                          + theme_minimal() + plotCommonGrid + plotCommonTheme
                          + scale_color_manual(values = c("steelblue", "firebrick"))
                          + theme(axis.ticks.y = element_blank()))

            return (comp.plot)
          })


#
# Plotting functions
#
setGeneric(name = "plotHeatmapCorrelationMatrixDGE",
           def = function(object) {standardGeneric("plotHeatmapCorrelationMatrixDGE")})
setMethod(f = "plotHeatmapCorrelationMatrixDGE",
          signature = "MixedSpeciesSample",
          function(object) {
            heatmap_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))
            heatmap.2(cor(as.matrix(object@dge@dge)), trace = "none", labRow = F,
                      labCol = F, dendrogram = "row", col=heatmap_palette(20))
          })

setGeneric(name = "plotCellTypes",
           def = function(object) {standardGeneric("plotCellTypes")})
setMethod(f = "plotCellTypes",
          signature = "data.frame",
          function(object) {
            cell.types.plot <- (ggplot(data = object[object$species != "undefined", ],
                                       aes_string(names(object)[2], names(object)[3], col = names(object)[4]))
                                + geom_point(size = 4, alpha = 0.4) + theme_classic()
                                + xlab(paste(names(object)[2], "transcripts"))
                                + ylab(paste(names(object)[3], "transcripts"))
                                + scale_color_manual(values = c("steelblue", "purple", "firebrick"),
                                                     labels = c(paste0(names(object)[2], " (", table(object$species)[1], ")"),
                                                                paste0("mixed (", table(object$species)[2], ")"),
                                                                paste0(names(object)[3], " (", table(object$species)[3], ")")))
                                + geom_point(data = object[object$species == "undefined", ],
                                             aes_string(names(object)[2], names(object)[3], col = names(object)[4]),
                                             col = "grey", size = 4, alpha = 0.4)
                                + plotCommonTheme
                                + theme(legend.title = element_blank(),
                                        legend.position = c(0.8, 0.8))
            )
            return(cell.types.plot)
          })

setGeneric(name = "plotGenesHistogram",
           def = function(object) {standardGeneric("plotGenesHistogram")})
setMethod(f = "plotGenesHistogram",
          signature = "data.frame",
          function(object) {
            species1 = names(table(object$species))[1]
            species2 = names(table(object$species))[2]

            plotLocalCommon <- (theme_minimal() + plotCommonGrid + plotCommonTheme
                                + theme(axis.text.y = element_blank(), axis.ticks = element_blank()))

            genes.species1 <- (ggplot(object[object$species == species1, ], aes(counts))
                               + xlab(paste(species1, "genes per", species1, "cell")) + ylab("")
                               + geom_histogram(aes(y = ..density..), fill = "steelblue", alpha = 0.8)
                                + geom_density(col = "black", size = 1) + plotLocalCommon)

            genes.species2 <- (ggplot(object[object$species == species2, ], aes(counts))
                               + xlab(paste(species2, "genes per", species2 , "cell")) + ylab("")
                               + geom_histogram(aes(y = ..density..), fill = "firebrick", alpha = 0.8)
                               + geom_density(col = "black", size = 1) + plotLocalCommon)

            return (grid.arrange(genes.species1, genes.species2, ncol = 2))
          })

setGeneric(name = "plotTranscriptsHistogram",
           def = function(object) {standardGeneric("plotTranscriptsHistogram")})
setMethod(f = "plotTranscriptsHistogram",
          signature = "data.frame",
          function(object) {
            species1 = names(table(object$species))[1]
            species2 = names(table(object$species))[2]

            plotLocalCommon <- (theme_minimal() + plotCommonGrid + plotCommonTheme
                                + theme(axis.text.y = element_blank(), axis.ticks = element_blank()))

            transcripts.species1 <- (ggplot(object[object$species == species1, ], aes(counts))
                               + xlab(paste(species1, "transcripts per", species1, "cell")) + ylab("") + plotLocalCommon
                               + geom_histogram(aes(y = ..density..), fill = "steelblue", alpha = 0.8)
                               + geom_density(col = "black", size = 1))

            transcripts.species2 <- (ggplot(object[object$species == species2, ], aes(counts))
                               + xlab(paste(species2, "transcripts per", species2, "cell")) + ylab("") + plotLocalCommon
                               + geom_histogram(aes(y = ..density..), fill = "firebrick", alpha = 0.8)
                               + geom_density(col = "black", size = 1))

            return (grid.arrange(transcripts.species1, transcripts.species2, ncol = 2))
          })





