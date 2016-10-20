#' Merges multiple samples into a big DGE matrix
#'
#' @param object A \code{list} of DGEs.
#' @return A \code{data.frame} containing all DGEs merged.
setGeneric("mergeMultipleDGEs",
           function(object) {standardGeneric("mergeMultipleDGEs")})
setMethod("mergeMultipleDGEs",
          "list",
          function(object) {
            list.of.dge <- list()

            for (sample in 1:length(object)) {
              list.of.dge[[sample]] <- cbind(object[[sample]]@dge, "genes"=as.character(object[[sample]]@genes))
            }

            merged.df <- Reduce(function(...) merge(..., by="genes", all=T), list.of.dge)
            rownames(merged.df) <- merged.df[, 1]
            merged.df <- merged.df[, 2:dim(merged.df)[2]]
            merged.df[is.na(merged.df)] <- 0
            merged.df <- merged.df[order(rownames(merged.df)), ]
            return (merged.df)
          })

#' Computes the inflection point / cell number.
#'
#' @param object A vector containing the reads per cell as given by the Drop-seq pipeline (BAMTagHistogram).
#' @param max.cells The number of cells to consider as an upper bound and to compute the inflection point.
#' @param n The nubmer of inflection points to return.
#' @return The cell number.
setGeneric("findInflectionPoint",
           function(object, max.cells=10000, n=1) {standardGeneric("findInflectionPoint")})
setMethod("findInflectionPoint",
          "numeric",
          function(object, max.cells, n) {
            dvec = c()

            # Normalize the vector
            object = cumsum(object[1:max.cells])
            object = object/max(object)
            L = length(object)

            # The distance of the endpoints
            dmax = dist(rbind(c(1/L, object[1]), c(1, object[L])), method="euclidean")

            # Find the point with maximum distance (minimum angle)
            for (i in 1:L) {
              d1 = dist(rbind(c(1/L, object[1]), c(i/L, object[i])), method="euclidean")
              d2 = dist(rbind(c(i/L, object[i]), c(1, object[L])), method="euclidean")

              if (d1 == 0 | d2 == 0) {next}   # Avoid singularities

              dvec = c(dvec, abs((d1^2 + d2^2 - dmax^2)/(2*d1*d2)))
            }
            return (order(dvec)[1:n])
          })
