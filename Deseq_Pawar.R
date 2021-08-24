setwd("/Users/yalegenomecenter/Desktop")
countData <- read.csv("count_matrix.txt", header=T, row.names=1, sep="\t") 
head(countData)
library("DESeq2")
proband-sibling
countData <- countData[,1:38]
colData <- DataFrame(condition=factor(c("proband","proband", "proband", "proband", "proband", "proband","proband","proband", "proband", "proband", "proband", "proband","proband","proband", "proband", "proband", "proband", "proband", "proband", "sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, file = "proband-sibling.csv", row.names = TRUE)


proband-control 
countData <- countData[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)]
colData <- DataFrame(condition=factor(c("proband","proband", "proband", "proband", "proband", "proband","proband","proband", "proband", "proband", "proband", "proband","proband","proband", "proband", "proband", "proband", "proband", "proband", "control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, file = "proband-control.csv", row.names = TRUE)


sibling-control
countData <- countData[,c(20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)]
colData <- DataFrame(condition=factor(c("sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling", "sibling","sibling","sibling", "sibling", "control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, file = "sibling-control.csv", row.names = TRUE)