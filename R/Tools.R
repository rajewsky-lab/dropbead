#' Merge multiple samples into a big DGE matrix
#'
#' Merge the DGEs of multiple samples into a big one that can be used for
#' further analysis.
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
              list.of.dge[[sample]] <- cbind(object[[sample]]@dge,
                                             "genes"=as.character(object[[sample]]@genes))
            }

            merged.df <- Reduce(function(...) merge(..., by="genes", all=T),
                                list.of.dge)
            rownames(merged.df) <- merged.df[, 1]
            merged.df <- merged.df[, 2:dim(merged.df)[2]]
            merged.df[is.na(merged.df)] <- 0
            merged.df <- merged.df[order(rownames(merged.df)), ]
            return (merged.df)
          })

#' Estimate the cell number
#'
#' Estimate the number of STAMPS that are present in the data (see Drop-seq
#' computational cookbook for further details).
#'
#' @param object A vector containing the reads per cell as given by the Drop-seq
#' pipeline (BAMTagHistogram).
#' @param max.cells The number of cells to consider as an upper bound and to
#' compute the inflection point (default is 15000).
setGeneric("estimateCellNumber",
           function(object, max.cells=15000) {
           standardGeneric("estimateCellNumber")
           })
setMethod("estimateCellNumber",
          "numeric",
          function(object, max.cells) {
             xdata <- (1:max.cells)/max.cells
             cs <- cumsum(df$V1[1:max.cells]/sum(df$V1[1:max.cells]))
             m <- (cs[max.cells] - cs[1]) / (max.cells - xdata[1])
             dists <- sapply(1:max.cells, function(cell) sin(atan(m)) * abs(cs[cell]-xdata[cell]))
             return (which.max(dists))
           })

#' Extract statistics from mapped reads.
#'
#' Extract statistics from mapped reads, such as total number of reads,
#' uniquely mapped reaqds, percentages of intronic/intergenic etc. The function
#' makes use of the R package Rsamtools (part of Bioconductor).
#'
#' @param object A string describing the exact location of the BAM file for the
#' sample under consideration.
#' @param yieldSize Large BAM files are processed in batches because of memory
#' issues. The parameter \code{yieldSize} controls the batch size. Default value: 10^6.
#' @param sample.name A string describing the sample.
#' @return A \code{data.frame} containing the computed statistics.
setGeneric("extractReadStatistics",
           function(object, yieldSize=10^6, sample.name=NULL) {
             standardGeneric("extractReadStatistics")
           })
setMethod("extractReadStatistics",
          "character",
          function(object, yieldSize, sample.name) {
            bam <- scanBam(BamFile(object), yieldSize=yieldSize,
                           param=ScanBamParam(
                             what='mapq',
                             tag=c('XF', 'GE', 'XC')
                           ))[[1]]

            non.ge.types <- table(bam$tag$XF[which(bam$mapq == 255 & is.na(bam$tag$GE))])

            df <- data.frame('total_reads'=round(length(bam$mapq)/10^6, 1),
                             'uniq_mapped'=round(sum(bam$mapq == 255, na.rm = T)/10^6, 1),
                             'intronic'=round(non.ge.types['INTRONIC']/10^6, 1),
                             'intergenic'=round(non.ge.types['INTERGENIC']/10^6, 1),
                             'coding'=round(non.ge.types['CODING']/10^6, 1),
                             'UTR'=round(non.ge.types['UTR']/10^6, 1),
                             'total_barcodes'=round(length(unique(bam$tag$XC[which(bam$mapq == 255)]))/10^6, 1))
            rownames(df) <- sample.name
            return (df)
          })

#' Extract downstream statistics
#'
#' Extract downstream basic statistics such as cell number, median numger of reads,
#' genes, UMIs per cell. The \code{data.table} package is used for fast reading
#' of the DGE matrices.
#'
#' @param object The exact location of the folder containing the analyzed data.
#' @param dge.name The file containing the DGE. If none is provided, the function assumes
#' that DGE is stored under the name \code{dge.txt.gz}.
#' @param reads.by.cell The \code{data.frame} containing the reads per cell. If none is
#' provided, the function tries to read the file \code{out_readcounts.txt.gz}.
#' @param min.umi The minimum number of UMIs to keep a cell and not discard it.
#' @param read.stats A \code{data.frame} containing read statistics, as it is produced
#' by the function \code{extractReadStatistics}.
setGeneric("extractDownstreamStatistics",
           function(object, dge.name=NULL, reads.by.cell=NULL, min.umi=250, read.stats=NULL) {
             standardGeneric("extractDownstreamStatistics")
           })
setMethod("extractDownstreamStatistics",
          "character",
          function(object, dge.name, reads.by.cell, min.umi, read.stats) {
            if (is.null(dge.name)) {
              dge <- data.frame(fread(paste0('zcat < ', object, 'dge.txt.gz')), row.names = 1)
            }
            else {
              dge <- data.frame(fread(paste0('zcat < ', object, dge.name)), row.names = 1)
            }
            if (is.null(reads.by.cell)) {
              reads.cell <- data.frame(fread(paste0('zcat < ', object, 'out_readcounts.txt.gz')))
            }
            else {
              reads.cell <- data.frame(fread(paste0('zcat < ', object, reads.by.cell)))
            }

            dge <- dge[, names(dge)[which(colSums(dge) >= min.umi)]]
            reads.cell <- reads.cell[reads.cell[, 2] %in% names(dge)[which(colSums(dge) >= min.umi)], ]

            df = data.frame('cells'=dim(dge)[2],
                            'reads'=round(median(reads.cell[, 1])),
                            'genes'=round(median(colSums(dge > 0))),
                            'umis'=round(median(colSums(dge))),
                            'sum.umi'=round(sum(dge)/10^6, 1),
                            'PCR'=round(median(reads.cell[, 1]/colSums(dge)), 1))
            if (!is.null(read.stats)){
              df$PCR <- round(median((reads.cell[, 1]*(1-sum(read.stats[, 3:6])/read.stats[, 2]))/colSums(dge)), 1)
            }
            return (df)
          })
