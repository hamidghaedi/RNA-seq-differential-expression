#________Packgaes______________#
library("TCGAbiolinks") # bioconductor package
library("SummarizedExperiment") # bioconductor package
library("DESeq2") # bioconductor package
library("IHW") # bioconductor package
library("biomaRt") # bioconductor package
library("apeglm") # bioconductor package
library("pheatmap") # CRAN package
library("RColorBrewer") # CRAN package
library("PCAtools") # bioconductor package
library(reshape2) # CRAN package


#_______ Donloading_Data_______#
query_TCGA = GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

GDCdownload(query = query_TCGA)

dat <- GDCprepare(query = query_TCGA, save = TRUE, save.filename = "exp.rda")


# exp matrix
rna <- as.data.frame(SummarizedExperiment::assay(dat))
# clinical data
clinical <- data.frame(dat@colData)

# Check how many tumor and control sample are there:
## data are stored under "definition" column in the clinical dataset.
table(clinical$definition)

# Also from sample id it is possible to count normal and tumor samples:
table(substr(colnames(rna),14,14))

#The count matrix (rna dataset) and the rows of the column data (clinical dataset) MUST be in the same order.

# to see whether all rows of clinical are present in rna datset
all(rownames(clinical) %in% colnames(rna))

# whether they are in the same order:
all(rownames(clinical) == colnames(rna))

# if not reorder them by:
#rna <- rna[, rownames(clinical)]
#all(rownames(clinical) == colnames(rna))

#_______Making_Expression_Object__________#

#We will use the column “definition”, as the grouping variable for gene expression analysis. 

# replace spaces with "_" in levels of definition column
clinical$definition <-  gsub(" ", "_", clinical$definition)

# making the definition column as factor
clinical$definition <- as.factor(clinical$definition)
# relevling factor to ensure tumors would be compared to normal tissue.
levels(clinical$definition)
#
clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")

# Making DESeqDataSet object which stores all experiment data
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = clinical,
                              design = ~ definition)
#dds

# prefilteration: it is not necessary but recommended to filter out low expressed genes

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# data tranfromation
vsd <- vst(dds, blind=FALSE)

# making PC object
p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)

# create PCA plot for PCA1 and PCA2
biplot(p, colby = "definition", lab = NULL, legendPosition = 'right')

# Fol all top 10 possible combination 
pairsplot(p,
          components = getComponents(p, c(1:10)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'definition',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
normal_idx <- substr(colnames(assay(vsd)),14,14) == "1"
n_sample <- assay(vsd)[, c(normal_idx) ]
colnames(n_sample) <- paste("NT_", substr(colnames(n_sample),1,12))

# Dissimilarity matrix calculation
sampleDists <- dist(t(n_sample))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- c(colnames(n_sample), colnames(t_sample))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap visualization
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#_______DE_analysis________#
 
dds <- DESeq(dds) #This would take some time

#default method
res <- results(dds, alpha = 0.05,  altHypothesis = "greaterAbs", lfcThreshold = 1.5) # alpha controls FDR rate

#result Log fold change shrinkage method (suitable for logfc based visualization)
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC.Ordered<-resLFC[with(resLFC, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

# converting Ensebl id to Gene symbole using biomart
ens2symbol<-function(ids){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- getBM(filters= "ensembl_gene_id", 
                 attributes= c("ensembl_gene_id","hgnc_symbol"),
                 values=ids, mart= mart)
  return(genes)
}

df <- ens2symbol(row.names(res))

res_df <- as.data.frame(res)                 
res_df$ensembl_gene_id <- row.names(res_df)
res_df <- merge(df,res_df, by = "ensembl_gene_id")
resOrdered<-res_df[with(res_df, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

#saving the results
write.csv(res_df, 
          file= paste0(resultsNames(dds)[2], ".csv")

          
#result with Independent hypothesis weighting
resIHW <- results(dds, filterFun=ihw, alpha = 0.05, altHypothesis = "greaterAbs", lfcThreshold = 1.5)
resIHW_df <- as.data.frame(resIHW)                 
resIHW_df$ensembl_gene_id <- row.names(resIHW_df)
resIHW_df <- merge(df,resIHW_df, by = "ensembl_gene_id")

resIHWOrdered <- resIHW_df[with(resIHW_df, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

write.csv(resIHW_df, 
          file= paste0("IHW",resultsNames(dds)[2], ".csv"))
#ploting genes differentially expressed 
ylim <- c(-6.5,6.5)
drawLines <- function() abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)
plotMA(resLFC, ylim=ylim); drawLines()

# obtaining normalized count read
dys_reg <- assay(vsd)[which(row.names(assay(vsd)) %in% resOrdered$ensembl_gene_id [1:25]), ]

#melting dataset for visualization
melted_norm_counts <-data.frame(melt(dys_reg))
colnames(melted_norm_counts) <- c("gene", "samplename", "normalized_counts")
melted_norm_counts$group <- ifelse(melted_norm_counts$samplename %in% colnames(assay(vsd))[normal_idx], "Normal", "Tumor")

# write.table to import in python
write.table(melted_norm_counts, file = "data.csv", sep = ",", row.names = F)

