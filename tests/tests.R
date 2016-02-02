library(data.table)
library(dropseq)

# load a DGE to test
d <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/macosko/SpeciesMix_ThousandSTAMPs_50cellspermicroliter/sampled136/dge.txt.gz"), row.names = 1)
m <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=d)
sh <- splitMixedSpeciesSampleToSingleSpecies(m, 0.9)[[1]]
sm <- splitMixedSpeciesSampleToSingleSpecies(m, 0.9)[[2]]

d50 <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/macosko/SpeciesMix_ThousandSTAMPs_50cellspermicroliter/SRR1748411/dge.txt.gz"), row.names = 1)
m50 <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=d50)

do <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/hek3t3Pilot002/hek3t3_pilot/NR_JS0001/hek3t3pilot2/dge.txt.gz"), row.names = 1)
mo <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=do)
sho <- splitMixedSpeciesSampleToSingleSpecies(mo, 0.9)[[1]]
smo <- splitMixedSpeciesSampleToSingleSpecies(mo, 0.9)[[2]]

d.hek.50 <- data.frame(fread("zcat < /mydaten/hek3t3fly/Rajewsky/NR_CK_020/ds_009_50/dge.txt.gz"), row.names = 1)
o.50 <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=d.hek.50)

createSamplesForTesting <- function() {
  genes <- c(sample(letters, 25), sample(paste0(letters, letters), 25),
             sample(LETTERS, 25), sample(paste0(LETTERS, LETTERS), 25))
  cells <- c()
  for (i in 1:510) {
    cells <- c(cells, paste(sample(c("A", "C", "G", "T"), 12, replace=T), collapse=""))
  }
  for (m in c(57, 172, 311, 489)) {
    cells[m] <- paste0(substr(cells[m-1], 1, 11), "N")
  }
  dge.human.cells <- data.frame(matrix(0, ncol=240, nrow=50))
  dge.mouse.cells <- data.frame(matrix(0, ncol=260, nrow=50))
  dge.full <- data.frame(matrix(0, ncol=510, nrow=100))
  for (i in 1:10) {
    dge.full[, i] <- round(runif(100, 0, 2*i))
  }
  for (i in 11:240) {
    dge.full[, i] <- c(sample(4, 50, replace = T, prob = c(0.6, .15, .15, .1))-1, sample(i/2, 50, replace = T)-1)
    dge.human.cells[, i-10] <- dge.full[51:100, i]
  }
  for (i in 241:250) {
    dge.full[, i] <- c(sample(4, 50, replace = T, prob = c(0.6, .15, .15, .1))-1, sample(i/2, 50, replace = T)-1)
    dge.human.cells[, i-10] <- dge.full[51:100, i]
  }
  for (i in 251:280) {
    dge.full[, i] <- c(sample((i-240)/2, 50, replace = T)-1, sample(4, 50, replace = T, prob = c(0.6, .15, .15, .1))-1)
    dge.mouse.cells[, i-250] <- dge.full[1:50, i]
  }
  for (i in 281:510) {
    dge.full[, i] <- c(sample((i-240)/2, 50, replace = T)-1, sample(4, 50, replace = T, prob = c(0.6, .15, .15, .1))-1)
    dge.mouse.cells[, i-250] <- dge.full[1:50, i]
  }
  rownames(dge.full) <- genes
  rownames(dge.human.cells) <- genes[51:100]
  rownames(dge.mouse.cells) <- genes[1:50]
  names(dge.full) <- cells
  names(dge.human.cells) <- names(dge.full)[11:250]
  names(dge.mouse.cells) <- names(dge.full)[251:510]

  single.human <- new("SingleSpeciesSample", species1="human", dge=dge.human.cells)
  single.mouse <- new("SingleSpeciesSample", species1="mouse", dge=dge.mouse.cells)
  mixed <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=dge.full)

  return (list(mixed, single.human, single.mouse))
}

testSplittingMixedToSingleSpecies <- function() {
  in.silico.samples <- createSamplesForTesting()
  mixed <- in.silico.samples[[1]]
  single.human <- in.silico.samples[[2]]
  single.mouse <- in.silico.samples[[3]]

  # Test number of cells for each species
  if (!table(classifyCellsAndDoublets(mixed, 0.6, min.trans=0)$species)[[1]] == dim(single.human@dge)[2]) {
    print ("Error! Number of human cells doesn't match.")
  } else {
      print ("Testing number of human cells... Passed")
    }
  if (!table(classifyCellsAndDoublets(mixed, 0.6, min.trans=0)$species)[[3]] == dim(single.mouse@dge)[2]) {
    print ("Error! Number of mouse cells doesn't match.")
  } else {
    print ("Testing number of mouse cells... Passed")
  }
  # Test number of genes per cell per species
  if (sum(computeGenesPerCell(mixed, 0.6, min.reads = 0)[1:240, 'counts'] == computeGenesPerCell(single.human, min.reads = 0)$counts) != 240) {
    print ("Error! The number of genes in human cells don't match")
  } else {
    print ("Testing number of genes in human cells... Passed")
  }
  if (sum(computeGenesPerCell(mixed, 0.6, min.reads=0)[241:500, 'counts'] == computeGenesPerCell(single.mouse, min.reads = 0)$counts) != 260) {
    print ("Error! The number of genes in mouse cells don't match")
  } else {
    print ("Testing number of genes in mouse cells... Passed")
  }
  # Test number of transcripts per cell per species
  if (sum(computeTranscriptsPerCell(mixed, 0.6)[1:240, 'counts'] == computeTranscriptsPerCell(single.human)$counts) != 240) {
    print ("Error! The number of transcripts in human cells don't match")
  } else {
    print ("Testing number of transcripts in human cells... Passed")
  }
  if (sum(computeTranscriptsPerCell(mixed, 0.6)[241:500, 'counts'] == computeTranscriptsPerCell(single.mouse)$counts) != 260) {
    print ("Error! The number of transcripts in mouse cells don't match")
  } else {
    print ("Testing number of transcripts in mouse cells... Passed")
  }
  # Testing collapsing of cells by barcodes
  collapseDf <- classifyCellsAndDoublets(collapseCellsByBarcode(mixed, 0.6), 0.6)
  collapseDfS <- collapseCellsByBarcode(single.human)@cells
  if (sum(collapseDf[collapseDf$species == "human", ]$cell ==  collapseCellsByBarcode(single.human)@cells) == 238) {
    print ("Testing collapsing human cells by barcode... Passed")
  } else {
      print ("Error! Collapsing human cells by barcodes is not the same for single and mixed species.")
  }
  if (sum(collapseDf[collapseDf$species == "mouse", ]$cell ==  collapseCellsByBarcode(single.mouse)@cells) == 258) {
    print ("Testing collapsing mouse cells by barcode... Passed")
  } else {
    print ("Error! Collapsing mouse cells by barcodes is not the same for single and mixed species.")
  }
  }

testSplittingMixedToSingleSpecies()

plotViolin(computeGenesPerCell(m))
length(m@dge[, 1][m@dge[, 1] >= 10])
head(computeGenesPerCell(splitMixedSpeciesSampleToSingleSpecies(m)[[1]], min.reads = 1))

computeGenesPerCell(new("SingleSpeciesSample", species1="human", dge=m@dge[, 1, drop=F]))

hist(computeGenesPerCell(keepBestCells(o.50, 1000), min.reads = 2)$counts, breaks = 100, col ='red',
     xlab = 'genes per cell (10 reads required to call a gene expressed)', main='')


hist(computeGenesPerCell(m, min.reads = 10)$counts, breaks = 100, col ='red',
     xlab = 'genes per cell (10 reads required to call a gene expressed)', main='')


plotViolin(computeGenesPerCell(keepBestCells(m, min.num.trans = 500)))
plotCellTypes(classifyCellsAndDoublets(m, min.trans = 500))
plotCellTypes(classifyCellsAndDoublets(keepBestCells(m, min.num.trans = 500)))


mg1 <- data.frame(fread("/data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/macosko/SpeciesMix_ThousandSTAMPs_50cellspermicroliter/sampled136/dge_summary.txt"))
mg2 <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/macosko/SpeciesMix_ThousandSTAMPs_50cellspermicroliter/sampled136/out_readcounts.txt.gz"))
names(mg2) <- c("reads", "CELL_BARCODE")
mq <- merge(mg2, mg1, by="CELL_BARCODE")

g1 <- data.frame(fread("/mydaten/hek3t3fly/ZinzenWahle/NR_CK_019/ds_008_100/dge_summary.txt"))
g2 <- data.frame(fread("zcat < /mydaten/hek3t3fly/ZinzenWahle/NR_CK_019/ds_008_100/out_readcounts.txt.gz"))
names(g2) <- c("reads", "CELL_BARCODE")
q <- merge(g2, g1, by="CELL_BARCODE")

g1 <- data.frame(fread("/mydaten/hek3t3fly/Rajewsky/NR_CK_020/ds_009_50/dge_summary.txt"))
g2 <- data.frame(fread("zcat < /mydaten/hek3t3fly/Rajewsky/NR_CK_020/ds_009_50/out_readcounts.txt.gz"))
names(g2) <- c("reads", "CELL_BARCODE")
q <- merge(g2, g1, by="CELL_BARCODE")

compareGeneExpressionLevels(m, mo, threshold1=0.9, threshold2=0.9)

























