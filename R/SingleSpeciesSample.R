SingleSpeciesSample <- setClass(Class = "SingleSpeciesSample",
                                slots = c(species1 = "character",
                                          dge = "DigitalGeneExpressionMatrix")
                                )

setGeneric(name = "computeGenesPerCell",
           def = function(object, ...) {standardGeneric("computeGenesPerCell")})
setMethod(f = "computeGenesPerCell",
          signature = "SingleSpeciesSample",
          function(object) {
            genes <- data.frame("cells" = names(colSums(object@dge@dge != 0)),
                                "counts" = as.numeric(colSums(object@dge@dge != 0)),
                                "species" = object@species1)
            return (genes)
          })

setGeneric(name = "computeTranscriptsPerCell",
           def = function(object, ...) {standardGeneric("computeTranscriptsPerCell")})
setMethod(f = "computeTranscriptsPerCell",
          signature = "SingleSpeciesSample",
          function(object) {
            transcripts <- data.frame("cells" = names(colSums(object@dge@dge)),
                                      "counts" = as.numeric(colSums(object@dge@dge)),
                                      "species" = object@species1)
            return (transcripts)
          })

setMethod(f = "removeLowQualityGenes",
          signature = "SingleSpeciesSample",
          function(object) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityGenes(object@dge)))
          })

setMethod(f = "removeLowQualityCells",
          signature = "SingleSpeciesSample",
          function(object) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityCells(object@dge)))
          })

#
# Plotting methods
#
# Set common theme options
#
plotCommonTheme <- theme(text = element_text(family = "IPAPMincho", size = 24),
                         axis.title.y = element_text(vjust = 1.2),
                         axis.title.x = element_text(vjust = 0))

plotCommonGrid <- theme(panel.grid.major.x = element_blank(),
                        panel.grid.major.y = element_line(colour = "grey50", size = 0.1),
                        panel.grid.minor = element_blank())

setGeneric(name = "plotViolinGenes",
           def = function(object) {standardGeneric("plotViolinGenes")})
setMethod(f = "plotViolinGenes",
          signature = "data.frame",
          function(object) {
            v <- (ggplot(object, aes_string(x = names(object)[3], y = names(object)[2], fill = names(object)[3]))
                  + geom_violin() + scale_fill_manual(values = c("steelblue", "firebrick"))
                  + scale_y_continuous(limits = c(0, max(object[, 2]))) + ylab("Number of genes per cell")
                  + xlab("") + geom_boxplot(width = 0.1, outlier.size = 0.5) + guides(fill = F)
                  + theme_minimal() + plotCommonTheme + plotCommonGrid
                  + theme(axis.ticks = element_blank())
            )
            return (v)
          })

setGeneric(name = "plotViolinTranscripts",
           def = function(object) {standardGeneric("plotViolinTranscripts")})
setMethod(f = "plotViolinTranscripts",
          signature = "data.frame",
          function(object) {
            v <- (ggplot(object, aes_string(x = names(object)[3], y = names(object)[2], fill = names(object)[3]))
                  + geom_violin() + scale_fill_manual(values = c("steelblue", "firebrick"))
                  + scale_y_continuous(limits = c(0, max(object[, 2]))) + ylab("Number of Transcripts per cell")
                  + xlab("") + geom_boxplot(width = 0.1, outlier.size = 0.5) + guides(fill = F)
                  + theme_minimal() + plotCommonTheme + plotCommonGrid
                  + theme(axis.ticks = element_blank())
            )
            return (v)
          })
