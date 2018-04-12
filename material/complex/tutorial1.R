## ------------------------------------------------------------------------
sampleNames <- c("trapnell_counts_C1_R1", "trapnell_counts_C1_R2", "trapnell_counts_C1_R3", "trapnell_counts_C2_R1", "trapnell_counts_C2_R2", "trapnell_counts_C2_R3")

sampleFiles <- c("trapnell_counts_C1_R1.tab", "trapnell_counts_C1_R2.tab", "trapnell_counts_C1_R3.tab", "trapnell_counts_C2_R1.tab", "trapnell_counts_C2_R2.tab", "trapnell_counts_C2_R3.tab")

sampleConditions <- c("C1", "C1", "C1", "C2", "C2", "C2")

## ------------------------------------------------------------------------
sampleNames
sampleFiles
sampleConditions

## ------------------------------------------------------------------------
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)

sampleTable

## ---- message=FALSE------------------------------------------------------
library("DESeq2")

## ------------------------------------------------------------------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design= ~ condition)

## ------------------------------------------------------------------------
ddsHTSeq <- DESeq(ddsHTSeq)

## ------------------------------------------------------------------------
resHTSeq <- results(ddsHTSeq)

head(resHTSeq)

## ------------------------------------------------------------------------
orderedRes <- resHTSeq[ order(resHTSeq$padj), ]

write.csv(as.data.frame(orderedRes), file="trapnell_C1_VS_C2.DESeq2.csv")

## ------------------------------------------------------------------------
normCounts <- counts(ddsHTSeq, normalized=TRUE)

head(normCounts)

write.csv(as.data.frame(orderedRes), file="trapnell_normCounts.DESeq2.csv")

## ------------------------------------------------------------------------
plotDispEsts(ddsHTSeq)

## ------------------------------------------------------------------------
hist(resHTSeq$pvalue, breaks=0:50/50, xlab="p value", main="Histogram of nominal p values")

## ------------------------------------------------------------------------
plotMA(resHTSeq)

shrunk_res <- lfcShrink(dds = ddsHTSeq, res = resHTSeq, coef = 2)
plotMA(shrunk_res)

## ------------------------------------------------------------------------
plot(resHTSeq$log2FoldChange, -log10(resHTSeq$pvalue), xlab="log2 Fold-change", ylab="-log P-value", pch=20, cex=0.5)
points(resHTSeq$log2FoldChange[ resHTSeq$padj<0.05 ], -log10(resHTSeq$pvalue[ resHTSeq$padj<0.05 ]), col="red", pch=20, cex=0.5)
abline(v=0, h=-log10(0.05), lty="dashed", col="grey")

## ------------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)

plotPCA(vsd)

## ------------------------------------------------------------------------
rawCounts <- counts(ddsHTSeq, normalized = TRUE)

# headmap of distances
heatmap(as.matrix(dist(t(assay(vsd)))), 
        main="Clustering of euclidean distances", scale="none")

# heatmap of correlations
heatmap(as.matrix(cor(log10(rawCounts + 1), method="pearson")), 
        main="Clustering of Pearson correlations", scale="none")

## ---- fig.height=8, fig.width=5------------------------------------------
library(gplots)

diffgenes <- rownames(resHTSeq)[ which(resHTSeq$padj < 0.05) ]
diffcounts <- normCounts[ diffgenes, ]

heatmap.2(diffcounts, 
          labRow = "", 
          trace = "none", density.info = "none",
          scale = "row",
          lhei = c(1,5),
          distfun = function(x) as.dist(1 - cor(t(x))))

## ------------------------------------------------------------------------
library(pheatmap)

# select the 20 most differentially expressed genes
select <- row.names(orderedRes[1:20, ])

# transform the counts to log10
log10_normCounts <- log10(normCounts + 1)

# get the values for the selected genes
values <- log10_normCounts[ select, ]

pheatmap(values,
         scale = "none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_row = 8,
         annotation_names_col = FALSE,
         gaps_col = c(3,6),
         display_numbers = TRUE,
         number_format = "%.2f",         
         height=12,
         width=6)

## ------------------------------------------------------------------------
sampleFiles <- c("trapnell_counts_C1_R1.tab", "trapnell_counts_C1_R2.tab", "trapnell_counts_C1_R3.tab", "trapnell_counts_C2_R1.tab", "trapnell_counts_C2_R2.tab", "trapnell_counts_C2_R3.tab")

tabs <- lapply(sampleFiles, function(x) read.table(x, col.names = c("Gene", x)))
countdata <- Reduce(f = function(x, y) merge(x, y, by="Gene"), x = tabs)

head(countdata)

## ------------------------------------------------------------------------
library(edgeR)

#countdata <- read.table("trapnell_counts.tab",sep="\t",header=T)
#row.names(countdata) <- countdata$Gene
#countdata$Gene<-NULL
row.names(countdata) <- countdata$Gene
countdata$Gene<-NULL
mygroups <- c("C1","C1","C1","C2","C2","C2")

w <- which(rowSums(countdata) > 0)
countdata <- countdata[ w, ]

y <- DGEList(counts=countdata[,1:6], genes=rownames(countdata), group = mygroups)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
topgenes<-topTags(et,n=dim(countdata)[[1]])
plot(topgenes$table$logFC, -log10(topgenes$table$FDR), col=ifelse(topgenes$table$FDR<0.05,"red","black"),main="FDR volcano plot",xlab="log2FC",ylab="-log10(FDR)")
hist(topgenes$table$PValue,breaks=20,xlab="P Value",col="royalblue",ylab="Frequency",main="P-value distribution")

