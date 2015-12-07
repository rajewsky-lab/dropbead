library(data.table)
library(dropseq)

# load a DGE to test
d <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/macosko/hek_3t3_12_5/SRR1873277/dge.txt.gz"), row.names = 1)
m <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=d)
sh <- splitMixedSpeciesSampleToSingleSpecies(m, 0.9)[[1]]
sm <- splitMixedSpeciesSampleToSingleSpecies(m, 0.9)[[2]]

do <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/hek3t3Pilot002/hek3t3_pilot/NR_JS0001/hek3t3pilot2/dge.txt.gz"), row.names = 1)
mo <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=do, cells=names())
sho <- splitMixedSpeciesSampleToSingleSpecies(mo, 0.9)[[1]]
sho.f <- removeLowQualityCells(removeLowQualityGenes(sho))
smo <- splitMixedSpeciesSampleToSingleSpecies(mo, 0.9)[[2]]

plotCellTypes(classifyCellsAndDoublets(mo))
plotViolinGenes(computeGenesPerCell(m))

qm <- identifyDoublets(sh@dge)
q <- identifyDoublets(sho@dge)
q.f <- identifyDoublets(sho.f@dge)

rowSums(q)

hist(geneExpressionMean(sho), breaks = 200, col = 'red')
hist(geneExpressionMean(sh), breaks = 200, col = 'red')

sort(geneExpressionMean(sho))

rowSums(q.f)

dim(sho@dge)

q.f

q <- new("SingleSpeciesSample", species1 = "human", dge = do)
q@genes

