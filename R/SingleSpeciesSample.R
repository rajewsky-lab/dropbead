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
           function(object, ...) {standardGeneric("computeGenesPerCell")})
setMethod("computeGenesPerCell",
          "SingleSpeciesSample",
          function(object) {
            genes <- data.frame("cells" = names(colSums(object@dge != 0)),
                                "counts" = as.numeric(colSums(object@dge != 0)),
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

setMethod("removeLowQualityGenes",
          "SingleSpeciesSample",
          function(object, n) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityGenes(object@dge, n)))
          })

setMethod("removeLowQualityCells",
          "SingleSpeciesSample",
          function(object, n) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityCells(object@dge, n)))
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

            for (index in 1:length(listOfCells)) {
              object@dge <- cbind(object@dge, rowSums(object@dge[, listOfCells[[index]]]))
            }
            object@dge <- object@dge[, !names(object@dge) %in% unlist(listOfCells)]

            names(object@dge)[(length(names(object@dge)) -
                                     length(listOfCells) + 1):length(names(object@dge))] <- unlist(listOfCells)[seq(1, length(unlist(listOfCells)), 2)]
            object@cells <- names(object@dge)
            return (object)
          })

