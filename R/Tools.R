#' Merges multiple samples into a big DGE matrix
#
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

# Returns n best points of the elbow of a curve
findElbowOfCurve <- function(x, n) {
  dvec = c()

  # Normalize the vector
  x = x/max(x)
  L = length(x)

  # The distance of the endpoints
  dmax = dist(rbind(c(1/L, x[1]), c(1, x[L])), method="euclidean")

  # Find the point with maximum distance (minimum angle)
  for (i in 1:L) {
    d1 = dist(rbind(c(1/L, x[1]), c(i/L, x[i])), method="euclidean")
    d2 = dist(rbind(c(i/L, x[i]), c(1, x[L])), method="euclidean")

    if (d1 == 0 | d2 == 0) {next}       # Avoid singularities

    dvec = c(dvec, abs((d1^2 + d2^2 - dmax^2)/(2*d1*d2)))
  }
  return (order(dvec)[1:n])
}
