# Load the R package
library(dropbead)

# Cumulative fraction of reads as function of cell barcodes plots ('knee' plots)
# Assumes that the out.readcounts.txt.gz has been generated for each lbrary and read
# as live.reads, fixed.reads, fixed.1wk.reads, fixed.3wk.reads
cutoff=5000
plotCumulativeFractionOfReads(live.reads, cutoff=cutoff) + geom_vline(xintercept=730, linetype="longdash")
plotCumulativeFractionOfReads(fixed.reads, cutoff=cutoff) + geom_vline(xintercept=1261, linetype="longdash")
plotCumulativeFractionOfReads(fixed.1wk.reads, cutoff=cutoff) + geom_vline(xintercept=307, linetype="longdash")
plotCumulativeFractionOfReads(fixed.3wk.reads, cutoff=cutoff) + geom_vline(xintercept=1656, linetype="longdash")

# Load the 4 DGEs into dge.live, dge.fixed, dge.fixed.1wk, dge.fixed.3wk
# Create the mixed species samples
live <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=dge.live)
fixed <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=dge.fixed)
fixed.1wk <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=dge.fixed.1wk)
fixed.3wk <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=dge.fixed.3wk)

# Collapse cell barcodes that could have eluded due to barcode corrections
live <- collapseCellsByBarcode(live)
fixed <- collapseCellsByBarcode(fixed)
fixed.1wk <- collapseCellsByBarcode(fixed.1wk)
fixed.3wk <- collapseCellsByBarcode(fixed.3wk)

# Separation species plots
plotCellTypes(classifyCellsAndDoublets(live, min.trans=3500, threshold=0.9))
plotCellTypes(classifyCellsAndDoublets(fixed, min.trans=3500, threshold=0.9))
plotCellTypes(classifyCellsAndDoublets(fixed.1wk, min.trans=3500, threshold=0.9))
plotCellTypes(classifyCellsAndDoublets(fixed.3wk, min.trans=3500, threshold=0.9))

# Keep only cells with at least 3,500 UMIs
min.umi.threshold=3500
live <- keepBestCells(live, min.num.trans=min.umi.threshold)
fixed <- keepBestCells(fixed, min.num.trans=min.umi.threshold)
fixed.1wk <- keepBestCells(fixed.1wk, min.num.trans=min.umi.threshold)
fixed.3wk <- keepBestCells(fixed.3wk, min.num.trans=min.umi.threshold)

# Violin plots - genes
plotViolin(computeGenesPerCell(keepBestCells(live, num.cells=100), min.reads=1), attribute="genes")
plotViolin(computeGenesPerCell(keepBestCells(fixed, num.cells=100), min.reads=1), attribute="genes")
plotViolin(computeGenesPerCell(keepBestCells(fixed.1wk, num.cells=100), min.reads=1), attribute="genes")
plotViolin(computeGenesPerCell(keepBestCells(fixed.3wk, num.cells=100), min.reads=1), attribute="genes")

# Violin plots - UMIs
plotViolin(computeTranscriptsPerCell(keepBestCells(live, num.cells=100)), attribute="transcripts (UMIs)")
plotViolin(computeTranscriptsPerCell(keepBestCells(fixed, num.cells=100)), attribute="transcripts (UMIs)")
plotViolin(computeTranscriptsPerCell(keepBestCells(fixed.1wk, num.cells=100)), attribute="transcripts (UMIs)")
plotViolin(computeTranscriptsPerCell(keepBestCells(fixed.3wk, num.cells=100)), attribute="transcripts (UMIs)")

# Correlations of gene expression levels between single cell libraries
compareGeneExpressionLevels(live, fixed, threshold1=0.9, threshold2=0.9, name1="Live", name2="Fixed")
compareGeneExpressionLevels(fixed.1wk, fixed, threshold1=0.9, threshold2=0.9, name1="Fixed2: 1 week", name2="Fixed")
compareGeneExpressionLevels(fixed.3wk, fixed, threshold1=0.9, threshold2=0.9, name1="Fixed2: 3 weeks", name2="Fixed")
compareGeneExpressionLevels(fixed.1wk, fixed.3wk, threshold1=0.9, threshold2=0.9, name1="Fixed2: 1 week", name2="Fixed: 3 weeks")

# Split the samples to obtain the HEK cells
live.h <- splitMixedSpeciesSampleToSingleSpecies(live)[[1]]
fixed.h <- splitMixedSpeciesSampleToSingleSpecies(fixed)[[1]]
fixed.1wk.h <- splitMixedSpeciesSampleToSingleSpecies(fixed.1wk)[[1]]
fixed.3wk.h <- splitMixedSpeciesSampleToSingleSpecies(fixed.3wk)[[1]]

# Correlations with bulk
# Assume that a bulk sample has been loaded (data.frame with genes as rownames
# and RPKMs in the first column) as bulk
compareSingleCellsAgainstBulk(live.h, bulk.data=bulk, ylab="Live")
compareSingleCellsAgainstBulk(ds13.50fix.h, bulk.data=bulk, ylab="Fixed")
compareSingleCellsAgainstBulk(ds19.50.h, bulk.data=bulk, ylab="Fixed2: 1wk")
compareSingleCellsAgainstBulk(ds21.h, bulk.data=bulk, ylab="Fixed2: 3wk")

# Compute mitochondrial percentages
live.mtch <- computeMitochondrialPercentage(live.h)
fixed.mtch <- computeMitochondrialPercentage(fixed.h)
fixed.1wk.mtch <- computeMitochondrialPercentage(fixed.1wk.h)
fixed.3wk.mtch <- computeMitochondrialPercentage(fixed.3wk.h)

# Plot non-mitochondrial reads
non.mtch <- plotMitochondrialContent(list(100-live.mtch, 100-fixed.mtch, 100-fixed.1wk.mtch, 100-fixed.3wk.mtch),
                                     c("Live", "Fixed", "Fixed.12k", "Fixed.3wk"), log_scale=F)
non.mtch + ylab("% of non-mitochondrial reads") + scale_y_continuous(limits=c(70, 100))

# Plot correlation between cell number and number of genes and UMIs
# Plot will look different from that of the paper because the libraries are not downsampled
# to the same depth.
number.of.cells <- c(length(fixed.1wk.h@cells), length(live.h@cells), length(fixed.3wk.h@cells), length(fixed.h@cells))
genes <- data.frame("cells"=number.of.cells, "name"=c("1wk", "L", "3wk", "F"),
                    "median"=c(median(computeGenesPerCell(fixed.3wk.h, min.reads=1)$counts),
                               median(computeGenesPerCell(live.h, min.reads=1)$counts),
                               median(computeGenesPerCell(fixed.3wk.h, min.reads=1)$counts),
                               median(computeGenesPerCell(fixed.h, min.reads=1)$counts)))

umis <- data.frame("cells"=number.of.cells, "name"=c("1wk", "L", "3wk", "F"),
                   "median"=c(median(computeTranscriptsPerCell(fixed.1wk.h)$counts),
                              median(computeTranscriptsPerCell(live.h)$counts),
                              median(computeTranscriptsPerCell(fixed.3wk.h)$counts),
                              median(computeTranscriptsPerCell(fixed.h)$counts)))

(ggplot(genes, aes(cells, median)) + geom_point(size=3) + geom_text(aes(label=name), hjust=0.2, vjust=-0.7, size=9)
  + geom_smooth(method="glm") + theme_minimal() +  plotCommonTheme + ylab("Number of genes"))
(ggplot(umis, aes(cells, median)) + geom_point(size=3) + geom_text(aes(label=name), hjust=0.2, vjust=-0.7, size=9)
  +  geom_smooth(method="glm") + theme_minimal() + plotCommonTheme + ylab("Number of transcripts (UMIs)"))
