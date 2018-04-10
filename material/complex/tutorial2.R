## ------------------------------------------------------------------------
rawdata <- read.delim("edgeR_example1_Tuch.tab", sep = "\t")

## ------------------------------------------------------------------------
head(rawdata)
dim(rawdata)
summary(rawdata)

## ------------------------------------------------------------------------
rawcounts <- rawdata[, 2:7]
genes <- rawdata[, 1]

## ------------------------------------------------------------------------
all(is.na(rawcounts) == FALSE)

## ------------------------------------------------------------------------
colSums(rawcounts)

barplot(colSums(rawcounts), ylab="Total number of reads", xlab="Sample ID")

## ------------------------------------------------------------------------
library(ggplot2)
library(reshape2)

ggplot(melt(rawcounts, measure.vars = 1:6), aes(y=value, x=variable, col=variable)) + 
  geom_boxplot() +                                    # we want a boxplot
  scale_y_log10()                                     # y scale is log10

ggplot(melt(rawcounts, measure.vars = 1:6), aes(x=value, col=variable)) + 
  geom_density() +                                    # we want a density plot
  scale_x_log10()                                     # x scale is log10

## ------------------------------------------------------------------------
log10_rawcounts <- log10(rawcounts + 1)

cor(log10_rawcounts, method="pearson")

## ------------------------------------------------------------------------
heatmap(as.matrix(cor(log10_rawcounts, method="pearson")), 
        main="Clustering of Pearson correlations", scale="none")
heatmap(as.matrix(cor(log10_rawcounts, method="spearman")), 
        main="Clustering of Spearman correlations", scale="none")

## ------------------------------------------------------------------------
library(edgeR)

## ------------------------------------------------------------------------
y <- DGEList(counts=rawcounts, genes=genes)
y <- calcNormFactors(y)

y$samples

## ------------------------------------------------------------------------
cpms <- as.data.frame(cpm(y, normalized.lib.sizes = TRUE))

ggplot(melt(cpms, measure.vars = 1:6), aes(x=value, col=variable)) + 
  geom_density() +                                    # we want a density plot
  scale_x_log10()                                     # x scale is log10

ggplot(melt(cpms, measure.vars = 1:6), aes(y=value, x=variable, col=variable)) + 
  geom_boxplot() +                                    # we want a boxplot
  scale_y_log10()                                     # y scale is log10

## ------------------------------------------------------------------------
plotMDS(y)

## ------------------------------------------------------------------------
Tissue <- factor(c("N","T","N","T","N","T"))

design <- model.matrix(~Tissue)
rownames(design) <- colnames(y)

design

## ------------------------------------------------------------------------
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)

## ------------------------------------------------------------------------
summary(decideTestsDGE(lrt))

## ------------------------------------------------------------------------
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

## ------------------------------------------------------------------------
result <- as.data.frame(topTags(lrt, n = nrow(rawcounts)))

head(result)

write.table(result, file = "edgeR_Tuch_Tumor_vs_NonTumor.csv", sep="\t", row.names = FALSE)

## ------------------------------------------------------------------------
Patient <- factor(c(8, 8, 33, 33, 51, 51))
Tissue <- factor(c("N","T","N","T","N","T"))

design <- model.matrix(~Patient + Tissue)
rownames(design) <- colnames(y)

design

y <- DGEList(counts=rawcounts, genes=genes)
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)

de <- decideTestsDGE(lrt)
summary(de)

plotMD(lrt)
abline(h=c(-1, 1), col="blue")

result_paired <- as.data.frame(topTags(lrt, n = nrow(rawcounts)))
write.table(result_paired, file = "edgeR_Tuch_Tumor_vs_NonTumor_paired.csv", sep="\t", row.names = FALSE)

