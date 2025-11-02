# Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

#install.packages("future.callr")
# load the libraries if required
library(Seurat)
library(SeuratObject)
library(SeuratDisk) 
library(tidyverse)
library(ggplot2)
#library(sctransform)
library(glmGamPoi)
library(harmony)
library(metap)
library(multtest)
library(presto)
library(SingleR)
library(celldex)
#library(scRNAseq)
#library(SingleCellExperiment)
library(pheatmap)
options(future.globals.maxSize = 1e9)
#library(CellChat)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(NMF)
library(compiler)
library(parallel)
#library(future.callr)
#library(DoubletFinder)
options(stringsAsFactors = FALSE)
library(reticulate)
library(SCENIC)
library(SCopeLoomR)

#Load integrated object 
objs <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithfinalcelltypeannotations_seurat_30dims.rds")


DefaultAssay(objs) # if needed: DefaultAssay(objs) <- "SCT"
DefaultAssay(objs) <- "RNA" #need raw counts for SCENIC
DefaultAssay(objs) #recheck; should be RNA

objs <- JoinLayers(objs, overwrite = TRUE) #RNA assay should have all the counts joined now
DefaultAssay(objs) #recheck; should be RNA

# Save as h5Seurat
SaveH5Seurat(objs, filename = "6_harmony_integratedwithtSNE_seurat_30dims.h5Seurat")

# Convert to h5ad
Convert("6_harmony_integratedwithtSNE_seurat_30dims.h5Seurat",
        dest = "h5ad",
        assay = "RNA",          # or "RNA" if that’s your main assay
        overwrite = TRUE)


