library(data.table)
library(dropseq)

dm <- data.frame(fread("zcat < /data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/macosko/hek_3t3_12_5/SRR1873277/dge.txt.gz"), row.names = 1)
mspm <- new("MixedSpeciesSample", species1="human", species2="mouse", dge=dm)

sh <- splitMixedSpeciesSampleToSingleSpecies(mspm)[[1]]
sh.f <- removeLowQualityCells(removeLowQualityGenes(sh, n=31), n=4000)

head(nbt@data[, 1:2])
log(head(sh.f@dge[, 1:2])+1, 2)
dim(nbt@data)
dim(sh.f@dge)

qn <- geneExpressionVariability(sh, all=T)



pq <- prcomp(t(sh.f@dge[names(q), ]), scale. = T,  center = T)
plot(pq$x)



summary(p)


length(q)


plot(qn)

plot(geneExpressionMean(sh@dge)[sort(names(qn))], qn[sort(names(qn))])


head(qn[sort(names(qn))], 10)
head(geneExpressionMean(sh@dge)[sort(names(qn))], 10)

