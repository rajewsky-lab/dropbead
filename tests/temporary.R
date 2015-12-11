library(data.table)
library(dropseq)

dm <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/macosko/hek_3t3_12_5/SRR1873277/dge.txt.gz"), row.names = 1)
mspm <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=dm)

sh <- splitMixedSpeciesSampleToSingleSpecies(mspm)[[1]]
sh.f <- removeLowQualityCells(removeLowQualityGenes(sh, n=31), n=4000)



geneExpressionVariability(sh.f, do.plot=T)




