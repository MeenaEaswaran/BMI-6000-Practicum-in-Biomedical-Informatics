#Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

#Run for the first time
#BiocManager::install("glmGamPoi") #needed for sctransform faster run

#To fix errors
#install.packages("devtools")
#devtools::install_github("satijalab/seurat", ref = "fix/v.5.3.1")

#install UMAP if not already installed
#reticulate::py_install(packages='umap-learn')  

#install for the first time
#install.packages("harmony")

#install the first time 
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')

#requires presto installation to speed up the find markers setup-do it the first time
#install.packages("devtools")
#devtools::install_github("immunogenomics/presto")
#install.packages('SeuratObject')
#install.packages('Matrix')
#remotes::install_github("satijalab/seurat")

#install for the first time
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("scRNAseq")

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.21")
#remove.packages("GenomeInfoDb")

#BiocManager::install("GenomeInfoDb", force = TRUE)
#install.packages("pheatmap")
  

# load the libraries if required
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(harmony)
library(metap)
library(multtest)
library(presto)
library(SingleR)
library(celldex)
library(scRNAseq)
library(SingleCellExperiment)
library(pheatmap)
options(future.globals.maxSize = 1e9)
#library(DoubletFinder)

#Read Harmony integrated Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
objs<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithtSNE_seurat_30dims.rds")

#Check default assay
DefaultAssay(object = objs) 
#DefaultAssay(objs) <- "SCT" #set if not

#SingleR can use raw or normalized counts from our dataset
#to convert the Seurat object to a SingleCellExperiment object.
objs.sce = as.SingleCellExperiment(objs, assay = "SCT", layers="data" ) #this uses the normalized count (data) after scTransform
# to use raw counts
#sce = as.SingleCellExperiment(objs, assay = "SCT", layers="counts" )

#survey avaialable reference datasets in celldex
library(celldex)
surveyReferences() # https://bioconductor.org/packages/3.21/data/experiment/vignettes/celldex/inst/doc/userguide.html
#BiocManager::install("scrapper")
library(scrapper)

#blue encode first
blue_encode <- celldex::BlueprintEncodeData() #blueprint and ENCODE. #use this as it as stroma and immune cells
unique(blue_encode$label.main)
table(blue_encode$label.main)

unique(blue_encode$label.fine)
table(blue_encode$label.fine)

# 2) SingleR at the CLUSTER level (recommended first pass)
pred.clust <- SingleR(
  test     = objs.sce,
  ref      = blue_encode,
  labels   = blue_encode$label.main,       # or blue_encode$label.fine
  clusters = objs$seurat_clusters
)

# Inspect: cluster -> label mapping
print(pred.clust[, c("labels","pruned.labels")])

cl2lab <- setNames(pred.clust$pruned.labels, rownames(pred.clust))
labs_by_cell <- cl2lab[as.character(objs$seurat_clusters)]
names(labs_by_cell) <- Cells(objs)  # set names to real cell barcodes

objs <- AddMetaData(objs, metadata = labs_by_cell, col.name = "SingleR_ENCODE_cluster")

table(objs$seurat_clusters, objs$SingleR_ENCODE_cluster)

DimPlot(objs, reduction = "umap.harmony", group.by = "SingleR_ENCODE_cluster",
        label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle("UMAP — SingleR (BlueprintEncode, cluster-level)")

####################################################################################################
#Human primary cell atlas (HPCA) -works
HPCA <- celldex::HumanPrimaryCellAtlasData()
unique(HPCA$label.main)
table(HPCA$label.main)
unique(HPCA$label.fine)
table(HPCA$label.fine)

# 2) SingleR at the CLUSTER level (recommended first pass)
pred.clust_HPCA <- SingleR(
  test     = objs.sce,
  ref      = HPCA,
  labels   =  HPCA$label.main,       
  clusters = objs$seurat_clusters
)

# Inspect: cluster -> label mapping
print(pred.clust_HPCA[, c("labels","pruned.labels")])

cl2lab_1 <- setNames(pred.clust_HPCA$pruned.labels, rownames(pred.clust_HPCA))
labs_by_cell <- cl2lab_1[as.character(objs$seurat_clusters)]
names(labs_by_cell) <- Cells(objs)  # set names to real cell barcodes

objs <- AddMetaData(objs, metadata = labs_by_cell, col.name = "SingleR_HPCA_cluster")

table(objs$seurat_clusters, objs$SingleR_HPCA_cluster)

DimPlot(objs, reduction = "umap.harmony", group.by = "SingleR_HPCA_cluster",
        label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle("UMAP — SingleR (HPCA, cluster-level)")