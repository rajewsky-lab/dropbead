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

melvir <- data.frame(fread("zcat < ~/Desktop/tempFlies/sampled25/dge.txt.gz"), row.names = 1)
mv <- new("MixedSpeciesSample", species1="melanogaster", species2="virilis", dge=melvir)
plotCellTypes(classifyCellsAndDoublets(mv, 0.51))



dm <- data.frame(fread("zcat < macosko/hek_3t3_12_5/SRR1873277/dge.txt.gz"), row.names = 1)
dmT <- data.frame(fread("zcat < macosko/SpeciesMix_ThousandSTAMPs_50cellspermicroliter/SRR1748411/dge.txt.gz"), row.names = 1)
dmH <- data.frame(fread("zcat < macosko/SpeciesMix_HundredSTAMPs_50cellpermicroliter/SRR1748412/dge.txt.gz"), row.names = 1)
d1 <- data.frame(fread("zcat < hek3t3Pilot001/hek3t3_pilot/NR_JS0001/sample1/dge.txt.gz"), row.names = 1)
d2 <- data.frame(fread("zcat < hek3t3Pilot002/hek3t3_pilot/NR_JS0001/hek3t3pilot2/dge.txt.gz"), row.names = 1)
da <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003a/sample/dge.txt.gz"), row.names = 1)
db <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003b/sample/dge.txt.gz"), row.names = 1)
dc <- data.frame(fread("zcat < hek3t3Pilot003/hek3t3Pilot003/NR_JS0003c/sample/dge.txt.gz"), row.names = 1)

library(dropseq)

dgematrixmT <- new("DigitalGeneExpressionMatrix", dge = dmT)
dgematrix1 <- new("DigitalGeneExpressionMatrix", dge = d1)
dgematrix2 <- new("DigitalGeneExpressionMatrix", dge = d2)
dgematrixa <- new("DigitalGeneExpressionMatrix", dge = da)
dgematrixb <- new("DigitalGeneExpressionMatrix", dge = db)
dgematrixc <- new("DigitalGeneExpressionMatrix", dge = dc)

mspm <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = new("DigitalGeneExpressionMatrix", dge = dm))
mspmT <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixmT)
mspmH <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse",
             dge = new("DigitalGeneExpressionMatrix", dge = dmH))
msp1 <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrix1)
msp2 <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrix2)
mspa <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixa)
mspb <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixb)
mspc <- new("MixedSpeciesSample", species1 = "human", species2 = "mouse", dge = dgematrixc)2/83*100



plotCellTypes(classifyCellsAndDoublets(msp1, 0.9))
plotViolinGenes(computeGenesPerCell(mspa, 0.9))
compareGeneExpressionLevels(mspa, 0.9, msp1, 0.9)

plotCellTypes(classifyCellsAndDoublets(mspb, 0.9))
plotViolinGenes(computeGenesPerCell(mspb, 0.9))
compareGeneExpressionLevels(mspa, 0.9, mspb, 0.9)

plotCellTypes(classifyCellsAndDoublets(mspc, 0.9))
plotViolinGenes(computeGenesPerCell(mspc, 0.9))
compareGeneExpressionLevels(mspb, 0.9, mspc, 0.9)
?computeGenesPerCell
#### Testing
sh <- splitMixedSpeciesSampleToSingleSpecies(mspm, 0.9)[[1]]
sh.f <- removeLowQualityCells(removeLowQualityGenes(sh))
dim(sh@dge@dge)
dim(sh.f@dge@dge)

plotViolinGenes(computeGenesPerCell(sh))
plotViolinGenes(computeGenesPerCell(sh.f))

hist(apply(log(sh.f@dge@dge+1, 2), 1, var), breaks = 150, col = 'red')
hist(apply(log(sh.f@dge@dge+1, 2), 1, mean), breaks = 150, col = 'red')

hist(apply(log(sh.f@dge@dge+1, 2), 1, var)/apply(log(sh.f@dge@dge+1, 2), 1, mean),
     breaks = 150, col = 'red')




df <- identifyDoublets(sh.f@dge)



hist(colSums(df), breaks = 100, col = 'red')
colSums(df)

hist(apply(log(df, 2) - as.numeric(geneExpressionMean(new("DigitalGeneExpressionMatrix", dge=df))), 2, var),
     breaks = 100, col = 'red')
hist(apply(log(df, 2) - as.numeric(geneExpressionMean(new("DigitalGeneExpressionMatrix", dge=df))), 2, mean),
     breaks = 100, col = 'red')

df
q <- (log(df+1, 2) - apply(log(df+1, 2), 1, mean))
q



apply(q, 2, mean)[apply(q, 2, mean) > 1.5]
q[, names(apply(q, 2, mean)[apply(q, 2, mean) > 1.5])]

hist(apply(q, 2, mean), breaks = 100, col ='red')
boxplot(apply(q, 2, mean))$out

df


heatmap.2(cor(as.matrix(log(df+1,2))), trace = "none", labRow = F,
          labCol = F, dendrogram = "row", col=heatmap_palette(20))

hclust(cor(as.matrix(log(df+1,2))))
as.matrix(log(df+1,2))

cell.list <- list()
for (i in 2:200) {
  for (j in (i+1):201) {
    if (sum(abs(df[, 1] - (df[, i] + df[, j]))) < 750) {
      cell.list <- append(cell.list, list(c(i, j)))
    }
  }
}
cell.list


df[, 1]
df[, 5] + df[, 50]
df[, 19] + df[, 28]
df[, 9] + df[, 22]

sum(abs(df[, 1] - (df[, 5] + df[, 50])))
sum(abs(df[, 1] - (df[, 19] + df[, 22])))

sh.f@dge@dge[, 1]
sh.f@dge@dge[, 8] + sh.f@dge@dge[, 50] - sh.f@dge@dge[, 1]

t1 <- sh.f@dge@dge[!rownames(sh.f@dge@dge) %in% names(rowSums(sh.f@dge@dge[, c(1, 5, 50)])[rowSums(sh.f@dge@dge[, c(1, 5, 50)]) < 10]), c(1, 5, 50)]
t2 <- sh.f@dge@dge[!rownames(sh.f@dge@dge) %in% names(rowSums(sh.f@dge@dge[, c(1, 19, 28)])[rowSums(sh.f@dge@dge[, c(1, 19, 28)]) < 10]), c(1, 19, 28)]
t3 <- sh.f@dge@dge[!rownames(sh.f@dge@dge) %in% names(rowSums(sh.f@dge@dge[, c(1, 19, 22)])[rowSums(sh.f@dge@dge[, c(1, 19, 22)]) < 10]), c(1, 19, 22)]
t4 <- sh.f@dge@dge[!rownames(sh.f@dge@dge) %in% names(rowSums(sh.f@dge@dge[, c(1, 9, 22)])[rowSums(sh.f@dge@dge[, c(1, 9, 22)]) < 10]), c(1, 9, 22)]

hist(abs(t1[, 1] - t1[, 2] - t1[, 3]), breaks = 150, col = 'red')
hist(abs(t2[, 1] - t2[, 2] - t2[, 3]), breaks = 150, col = 'red')
hist(abs(t3[, 1] - t3[, 2] - t3[, 3]), breaks = 150, col = 'red')
hist(abs(t4[, 1] - t4[, 2] - t4[, 3]), breaks = 150, col = 'red')

summary(abs(t1[, 1] - t1[, 2] - t1[, 3]))
summary(abs(t2[, 1] - t2[, 2] - t2[, 3]))
summary(abs(t4[, 1] - t4[, 2] - t4[, 3]))

# testing normality
library("MVN")
mardiaTest(t(log(df[c(1:12, 15:18), ]+1, 2)), qqplot=T)
uniPlot(t(log(df[c(1:12, 15:18), ]+1, 2)), type = "histogram")
uniNorm(t(log(df+1, 2)))
par(mfrow =  c(1,1))
mvOutlier(t(log(df+1, 2)), method="quan", label = F)

