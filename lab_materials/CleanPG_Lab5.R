library(ggplot2)
library(BiocManager)
library (DESeq2)
library(tximport)

# Set Working Directory
wdir <- setwd("/Users/jordanmanchego/Desktop/")

# Read in count data
cd <- read.table('countdata.tsv', header = TRUE, sep = "\t")
cd$gene_id <- sub(".v6.1", "", cd$gene_id)
cd$gene_id <- sub("Soltu.DM.", "", cd$gene_id)
head(cd)

# Read in metadata
md <- read.table('metadata.tsv', header = TRUE, sep = "\t")
colnames(md)[1]="genotype"

md$genotype <- factor(md$genotype)
levels(md$genotype)

md$condition <- factor(md$condition)
levels(md$condition)

md$day_of_drought <- factor(md$day_of_drought)
levels(md$day_of_drought)

md


###############################################################################
# CREATE DESEQ2 OBJECT
dds <- DESeqDataSetFromMatrix(countData = cd,
                              colData = md,
                              design=~condition,
                              tidy = TRUE)


###############################################################################
# DATA NORMALIZATION

# Prefilter
dds #before
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds #after

dds <- estimateSizeFactors(dds)
dds$sizeFactor


###############################################################################
# QUALITY CONTROL AND VISUALIZATION

# Plot Estimated Dispersions
dds <- estimateDispersions(dds)
plotDispEsts(dds)

# Variance Stabilizing Transformation (vsd)
vsd <- vst(dds, blind=FALSE)

# Regularize Logarithm (rlog)
rld <- rlog(dds, blind=FALSE)

# Calculate Distance Between Samples
sampleDists <- dist(t(assay(vsd)))

# Plot Sample Heatmap
library(pheatmap)
library(RColorBrewer)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd), vsd$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Plot PCA
plotPCA(vsd, intgroup=c("genotype"))


###############################################################################
# CALLING DIFFERENTIAL EXPRESSION

# Wald Test
dds <- nbinomWaldTest(dds)
resultsTable <- as.data.frame(results(dds, contrast=c("condition","control","drought")))
nrow(na.omit(resultsTable[resultsTable$padj < 0.05,]))
# OUT: [1] 4943

head(resultsTable)

# Plotting Results: Global Changes
# Volcano/MA Plot
plotMA(results(dds, contrast=c("condition", "control", "drought")), ylim=c(-2,2))








