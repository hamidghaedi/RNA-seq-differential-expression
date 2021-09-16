# RNA-seq-differential-expression
#### Differential gene expression analysis on RNA-seq data

## RNA-seq

RNA-seq is emerging as a powerful genome-wide gene expression analysis tools in genomic studies. Technically an RNA-seq experiment start from library preparation, sequencing, and data analysis. Data analysis is usually start from quality control, read mapping to a reference genome. From RNA reads aligned to a reference genome, a count matrix will generate that is the main file for differential expression analysis. For differential gene expression analysis from this count RNA matrix, there are two main R package which are widely used in this field: [DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)  and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html).

In this article I will share my workflow for performing differential gene expression analysis on  RNA-seq  count matrix. The data used here are coming from bladder cancer data from TCGA (TCGA-BLCA ). Here is a brief overview on analysis steps:

1.	Obtaining  data
2.	Data transformation, quality control and data normalization
3.	Differential gene expression analysis 
4.	Visualization

### 1. Obtaining data

```R
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


#_______ Downloading_Data_______#
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
```
### 2. Data transformation, quality control and  normalization
For differential expression analysis, using the raw counts is just fine. However before performing any analysis it is good idea to get a general picture of  data shape and distribution. To this aim it is necessary to transform data. For data transformation we will use variance stabilizing transformations (VST) method. This methods is described in the [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). Transformed data is useful to do some quality control checking on expression data. Quality control metrics would be used to check *sample level QC* and *gene level QC*. The DESeq2 package workflow performs gene level QC automatically. This includes removing genes that are very unlikely to be differentially expressed in the study design. Such genes are possibly those with zero counts in all samples, those with a hugely skewed read count and low mean normalized counts. 

By sample level QC, we can get how well the biological replicates(samples in different groups) cluster together. Actually sample level QC isapplied to answer the questions “How are samples similar to each other in termsof expression values?” , “ Which samples are different from each other?”,  ”Do thesedifferences align with our expectation from study design?” There are two statistical tools useful  to gain insight into sample level QC: *PCA* and *Hierarchical clustering*. 

```R
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
```
![alt text](https://github.com/hamidghaedi/RNA-seq-differential-expression/blob/master/pc1.PNG)

As depicted by the above plot, samples are somehow clustered together based on the “definition” variable levels: “Solid_Tissue_Normal”, and “Primary_solid_Tumor”.  Its helpful to see how other PCs -other than PC1 and PC2 would clusters samples together. There are no PC combinations that can perfectly cluster samples based on the “definition” variable and this is may be because of how RNA-seq data are. 


```R
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
```

![alt text](https://github.com/hamidghaedi/RNA-seq-differential-expression/blob/master/pc2.png)

*Hierarchical clustering* is suitable to visualize similarities between samples with regard to gene expression profile. At the same time it will provide information about how gene expression profiles are different between samples. Here we have a total of 433 samples , so visualizing all these sample in a plot would not be very informative. Here I will compare 19 normal samples for making a sample-sample heatmap for visualizing hierarchical clustering 

```R
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
```
![alt text](https://github.com/hamidghaedi/RNA-seq-differential-expression/blob/master/hc.png)


### 3. Differential gene expression analysis
```R
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
```
### 4. Visualization
There are several way for gene expression analysis results  visualization. 

*MA-plot* visualize gene significantly (blue dots with p < 0.05) up-regulated (log2FC > 1.5) or down-regulated (log2FC < -1.5). 

```R
#ploting genes differentially expressed 
ylim <- c(-6.5,6.5)
drawLines <- function() abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)
plotMA(resLFC, ylim=ylim); drawLines()
```
![alt text](https://github.com/hamidghaedi/RNA-seq-differential-expression/blob/master/ma-plot.png)

*Top 25 dysregulated gene*
```R

# obtaining normalized count read
dys_reg <- assay(vsd)[which(row.names(assay(vsd)) %in% resOrdered$ensembl_gene_id [1:25]), ]

#melting dataset for visualization
melted_norm_counts <-data.frame(melt(dys_reg))
colnames(melted_norm_counts) <- c("gene", "samplename", "normalized_counts")
melted_norm_counts$group <- ifelse(melted_norm_counts$samplename %in% colnames(assay(vsd))[normal_idx], "Normal", "Tumor")

# write.table to import in python
write.table(melted_norm_counts, file = "data.csv", sep = ",", row.names = F)
```
It is very interesting to visualize gene count data by violin plot. It is very informative! At the time of writing , there is no easy way to plot values in R "dodge" fashion. So I exported data into a file to plot this file in ```python```. The jupyter notebook code could be find in the repostory.

```python
import pandas as pd
import seaborn as sns
from matplotlib.pyplot import figure

#reading data
data = pd.read_csv('data.csv') 
data.head()

## plotting data
figure(num=None, figsize=(14, 6), dpi=300, facecolor='w', edgecolor='k')

graph = sns.violinplot(x="gene", y="normalized_counts", hue="group", data=data, palette="Set3", split=True)
graph.set_xticklabels(graph.get_xticklabels(), rotation=45, horizontalalignment='right')
```
![alt text](https://github.com/hamidghaedi/RNA-seq-differential-expression/blob/master/violion.png)

