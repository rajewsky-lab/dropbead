library(data.table)

om <- data.frame(fread("zcat < macosko/hek_3t3_12_5/SRR1873277/out_readcounts.txt.gz"))
omT <- data.frame(fread("zcat < macosko/SpeciesMix_ThousandSTAMPs_50cellspermicroliter/SRR1748411/out_readcounts.txt.gz"))
o1 <- data.frame(fread("zcat < hek3t3Pilot001/hek3t3_pilot/NR_JS0001/sample1/out_readcounts.txt.gz"))
o2 <- data.frame(fread("zcat < hek3t3Pilot002/hek3t3_pilot/NR_JS0001/hek3t3pilot2/out_readcounts.txt.gz"))
o3a <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003a/sample/out_readcounts.txt.gz"))
o3b <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003b/sample/out_readcounts.txt.gz"))
o3c <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003c/sample/out_readcounts.txt.gz"))

sum(om$V1[1:590])/sum(om$V1)
sum(omT$V1[1:1027])/sum(omT$V1)
sum(o1$V1[1:439])/sum(o1$V1)
sum(o2$V1[1:449])/sum(o2$V1)
sum(o3a$V1[1:492])/sum(o3a$V1)
sum(o3b$V1[1:583])/sum(o3b$V1)
sum(o3c$V1[1:597])/sum(o3c$V1)

setwd("/data/BIO3/home/nkarais/Work/@@/dropseq_cell/data/")

dm <- data.frame(fread("zcat < macosko/hek_3t3_12_5/SRR1873277/dge.txt.gz"), row.names = 1)
d1 <- data.frame(fread("zcat < hek3t3Pilot001/hek3t3_pilot/NR_JS0001/sample1/dge.txt.gz"), row.names = 1)
d2 <- data.frame(fread("zcat < hek3t3Pilot002/hek3t3_pilot/NR_JS0001/hek3t3pilot2/dge.txt.gz"), row.names = 1)
da <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003a/sample/dge.txt.gz"), row.names = 1)
db <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003b/sample/dge.txt.gz"), row.names = 1)
dc <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003c/sample/dge.txt.gz"), row.names = 1)

library(dropseq)

dgematrixm <- new("DigitalGeneExpressionMatrix", dge = dm)
dgematrix1 <- new("DigitalGeneExpressionMatrix", dge = d1)
dgematrix2 <- new("DigitalGeneExpressionMatrix", dge = d2)
dgematrixa <- new("DigitalGeneExpressionMatrix", dge = da)
dgematrixb <- new("DigitalGeneExpressionMatrix", dge = db)
dgematrixc <- new("DigitalGeneExpressionMatrix", dge = dc)

mspm <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixm)
msp1 <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrix1)
msp2 <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrix2)
mspa <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixa)
mspb <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixb)
mspc <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixc)


plotCellTypes(classifyCellsAndDoublets(msp1, 0.9))
plotViolinGenes(computeGenesPerCell(mspa, 0.9))
compareGeneExpressionLevels(mspa, 0.9, msp1, 0.9)

plotCellTypes(classifyCellsAndDoublets(mspb, 0.9))
plotViolinGenes(computeGenesPerCell(mspb, 0.9))
compareGeneExpressionLevels(mspa, 0.9, mspb, 0.9)

plotCellTypes(classifyCellsAndDoublets(mspc, 0.9))
plotViolinGenes(computeGenesPerCell(mspc, 0.9))
compareGeneExpressionLevels(mspb, 0.9, mspc, 0.9)

#### Testing

plotCellTypes(classifyCellsAndDoublets(msp1, 0.9))
t1 <- computeGenesPerCell(msp1, 0.9)
t2 <- computeTranscriptsPerCell(msp1, 0.9)

plotViolinGenes(computeGenesPerCell(msp1, 0.9))

s.t <- splitMixedSpeciesSampleToSingleSpecies(msp1, 0.9)[[2]]


quantile(rowSums(s.t@dge@dge), 0.9)


subset(s.t@dge@dge, rowSums(s.t@dge@dge) >= 247)








