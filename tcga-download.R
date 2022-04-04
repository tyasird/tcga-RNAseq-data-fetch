## Install packages ##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("biomaRt")
install.packages("xlsx", dependencies = TRUE)
install.packages("DT", dependencies = TRUE)
install.packages("SummarizedExperiment", dependencies = TRUE)

## Load packages ##
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(xlsx)
library(biomaRt)


## Obtain TCGA dataset ##
query <- GDCquery(project = "TCGA-THCA",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          file.type  = "htseq.counts",
                          experimental.strategy = "RNA-Seq")

GDCdownload(query, method = "api")
data <- GDCprepare(query)

THCAMatrix <- assay(data)

thcaexp <- as.data.frame(THCAMatrix)
write.csv2(thcaexp, file =  "THCAexpressiondata.csv",col.names = TRUE, row.names = TRUE, sep = "")


# id conversion ##
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensemblids = row.names(THCAMatrix)[1:nrow(thcaexp)]
genesymbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = ensemblids,
                     mart = ensembl)
