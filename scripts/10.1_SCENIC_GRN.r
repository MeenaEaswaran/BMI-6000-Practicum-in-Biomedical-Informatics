# Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

# SCENIC (Single cell regulatory network interference and clustering)

#Introduction to SCENIC
# SCENIC is a tool to simultaneously reconstruct gene regulatory networks and identify stable cell states from single-cell RNA-seq data. The gene regulatory network is inferred based on co-expression and DNA motif analysis, and then the network activity is analyzed in each cell to identify the recurrent cellular states.

#Required, SCENIC R is based on three R packages (for first time only)
BiocManager::install(c("AUCell", "RcisTarget"), force=TRUE)
BiocManager::install(c("GENIE3"), force = TRUE) # Can be replaced by GRNBoost (Linux & MacOS)


#Optional (but highly recommended):
#To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools"), force = TRUE)
#For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"), force = TRUE)
#To support parallel execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"), force = TRUE)

remotes::install_github("bokeh/rbokeh", force = TRUE)

#SCENIC
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC", force = TRUE) 
packageVersion("SCENIC")

devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE, force = TRUE)


#Species-specific databases: human, mouse and Drosophila (use human)
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs
#These should be downloaded from the old folders in cisTarget databases as the SCENIC R version precedes the latest pySCENIC version

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

#install cellchatDB and dependencies for the first time
#devtools::install_github("jinworks/CellChat")
#install.packages('NMF')
#devtools::install_github("jokergoo/circlize")
#devtools::install_github("jokergoo/ComplexHeatmap")
#UMAP-learn needed but already installed
#install.packages("future.callr")
# load the libraries if required
library(Seurat)
library(SeuratObject)
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

DimPlot(objs, reduction = "umap.harmony", label = TRUE)
##Join layers
objs <- JoinLayers(objs, overwrite = TRUE) #RNA assay should have all the counts joined now
DefaultAssay(objs) #recheck; should be RNA

#Get sparse count matrix
exprMat <- GetAssayData(objs, assay='RNA', slot = 'counts')
exprMat <- as(exprMat, "dgCMatrix")
head(exprMat)

#If exprMAt gives vector memory limit issue, do the below, IF REQUIRED.

# SINCE THE VECTOR MEMORY EXHAUST LIMIT REACHES WHEN RUNNING PCA FOR TUMOR SAMPLES, RSTUDIO LIMIT HAS TO BE RESET USING TERMINAL IN RSTUDIO AND THEN APPLICATION HAS TO BE RESTARTED. FOLLOW THE CODE/INSTRUCTIONS ON STACKOVERFLOW AS SHOWN BELOW

#https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos

library(usethis) 
usethis::edit_r_environ()
#when the tab opens up in R studio, add this to the 1st line: R_MAX_VSIZE=100Gb (or whatever memory you wish to allocate).
#Re-start R and/or restart computer and run the R command again that gave you the memory error.
#RERUN ABOVE COMMANDS after loading relevant libraries.

# cell information
Idents(objs)
cellInfo <- data.frame(seuratCluster=Idents(objs))

# Add common fields if they exist in meta.data
wanted_cols <- c("Condition", "Patient", "Sample", "SingleR_ENCODE", "seurat_clusters", "orig.ident")
have_cols   <- intersect(wanted_cols, colnames(objs@meta.data))
have_cols
if (length(have_cols) > 0) {
  cellInfo <- cbind(cellInfo, objs@meta.data[, have_cols, drop = FALSE])
}

# Ensure exprMat and cellInfo are properly aligned
# (cellInfo rows must match exprMat column names)
common <- intersect(colnames(exprMat), rownames(cellInfo))
exprMat  <- exprMat[, common, drop = FALSE]
cellInfo <- cellInfo[common, , drop = FALSE]

#Select genes expressed in more than 1% of all cells in the expression matrix.
loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
exprMat_filter <- exprMat[loci1, ]

#get a 2D embedding for SCope (Harmony UMAP)
emb_mat <- Embeddings(objs, "umap.harmony")
# align to loom column order
emb_mat <- emb_mat[colnames(exprMat_filter), , drop = FALSE]
default.embedding <- data.frame(
  X = emb_mat[,1],
  Y = emb_mat[,2],
  row.names = rownames(emb_mat))
default.embedding.name <- "UMAP (Harmony)"

# ---- Build the loom (with Harmony UMAP so SCope can visualize)
loom_path <- "~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/10_scenic_input.loom"
loom <- build_loom(
  file.name               = loom_path,
  dgem                    = exprMat_filter,          # sparse dgCMatrix
  title                   = "LSCC scRNA-seq",
  genome                  = "hg19",                  # switch to "mm10" if mouse
  default.embedding       = default.embedding,
  default.embedding.name  = default.embedding.name,
  chunk.size              = 1000,
  display.progress        = TRUE,
  loom.spec.version       = 3
)

#Add cell annotations
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)