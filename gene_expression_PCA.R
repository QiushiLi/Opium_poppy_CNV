library("DESeq2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library/")
sessionInfo()
packageVersion("DESeq2")

setwd("./")

# read in raw count matrix generaed by stringtie
countMatrix <- read.table ("gene_counts_matrix.txt", header=TRUE, row.names=1)

head(countMatrix)

#prepare info matrix

#form coldata

coldata = data.frame(row.names = c('PS_MC', 'PS_S1', 'PS_S2', 'PS_S3', 'PS_S4', 'PS_S5', 'PS_S6', 'PS_S7', 'PS_YC', 'PS_leaf', 'PS_root', 'PS_stem', 'PslatexRNA1'), group=c('PS_MC', 'PS_S1', 'PS_S2', 'PS_S3', 'PS_S4', 'PS_S5', 'PS_S6', 'PS_S7', 'PS_YC', 'PS_leaf', 'PS_root', 'PS_stem', 'PslatexRNA1'), lane=c('5','3','3','3','3','4','4','4','5','5','5','4','6'))


#model batch effect

dds<-DESeqDataSetFromMatrix(countData=countMatrix,colData=coldata,design = ~lane)
dds <- dds[ rowSums(counts(dds)) > 0, ]

#matrix of normalized counts


vsd <- varianceStabilizingTransformation(dds)
matrix2 <- (assay(vsd))
write.table(matrix2, file = "DESeq_NormalizedCounts_vsd.txt", sep = "\t", quote=F)


#principal component plot of samples

pdf("PCA_samples_vsd.pdf", width = 8, height = 8)
plotPCA(vsd, intgroup=c("group"))
dev.off()

