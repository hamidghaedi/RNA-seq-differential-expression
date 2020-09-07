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
```
### 2. Data transformation, quality control and  normalization
For differential expression analysis, using the raw counts is fine. However before performing any analysis it is good idea to have a general picture about the data shape and pattern. To this aim it is necessary to transform data. For data transformation we will use variance stabilizing transformations (VST). This methods is finely described in the [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). Transformed data is helpful to do some quality control checking onexpression data. Quality control metrics would be used to check *sample level QC* and *gene level QC*. The DESeq2 package workflow performs gene level QCautomatically. This includes removing genes that are very unlikely to be differentiallyexpressed in the study design. Such genes are possibly those with zero countsin all samples, those with a hugely skewed read count and low mean normalizedcounts. 

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

```R
##An eigencor plot could help to find correlation between different PC and clinical variables
eigencorplot(p, metavars = c("shortLetterCode","paper_mRNA.cluster","paper_Age.at.diagnosis","paper_AJCC.pathologic.tumor.stage", "race","gender"))
# sample to sample heatmap
sampleDists <- dist(t(assay(vsd)))

#Heat map only for normal sammples
normal_idx <- substr(colnames(assay(vsd)),14,14) == "1"
norm_sample <- assay(vsd)[, normal_idx]
sampleDists <- dist(t(norm_sample))
```
