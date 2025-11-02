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

# load the libraries if required
library(Seurat)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(harmony)
options(future.globals.maxSize = 1e9)
#library(DoubletFinder)

# Read CCA integrated Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
objs<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_CCA_integrated_seurat_30dims.rds")

#Check default assay
DefaultAssay(object = objs) 

#Segregation of clusters by cell cycle phase
# Explore whether clusters segregate by cell cycle phase
DimPlot(objs,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Save image
ggsave(
  filename = "5_CCAintegratedUMAP_coloredbycellcyclescores.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

#Segregation of clusters by various sources of metrics (inital QC metrics) in CCA integrated data
metrics <-  c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "S.Score", "G2M.Score")
FeaturePlot(
  objs,
  reduction  = "umap.cca", 
  features   = metrics,
  pt.size    = 0.4, 
  order      = TRUE, # <- replace sort.cell with order
  label     = TRUE, 
  min.cutoff = "q10"
)

# Save image
ggsave(
  filename = "5_CCAintegratedUMAP_coloredbyinitialQCmetrics.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

rm(objs) #before running the same with RPCA integrated data

#######################################################################

# Read RPCA integrated Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
objs<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_RPCA_integrated_seurat_30dims.rds")

#Check default assay
DefaultAssay(object = objs) 

#Segregation of clusters by cell cycle phase
# Explore whether clusters segregate by cell cycle phase
DimPlot(objs,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Save image
ggsave(
  filename = "5_RPCAintegratedUMAP_coloredbycellcyclescores.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

#Segregation of clusters by various sources of metrics (inital QC metrics) in RPCA integrated data
metrics <-  c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "S.Score", "G2M.Score")
FeaturePlot(
  objs,
  reduction  = "umap.RPCA", 
  features   = metrics,
  pt.size    = 0.4, 
  order      = TRUE, # <- replace sort.cell with order
  label     = TRUE, 
  min.cutoff = "q10"
)

# Save image
ggsave(
  filename = "5_RPCAintegratedUMAP_coloredbyinitialQCmetrics.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

rm(objs) #before running the same with harmony integrated data

#######################################################################

# Read harmony integrated Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
objs<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_harmony_integrated_seurat_30dims.rds")

#Check default assay
DefaultAssay(object = objs) 

#Segregation of clusters by cell cycle phase
# Explore whether clusters segregate by cell cycle phase
DimPlot(objs,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Save image
ggsave(
  filename = "5_harmonyintegratedUMAP_coloredbycellcyclescores.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

#Segregation of clusters by various sources of metrics (inital QC metrics) in harmony integrated data
metrics <-  c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "S.Score", "G2M.Score")
FeaturePlot(
  objs,
  reduction  = "umap.harmony", 
  features   = metrics,
  pt.size    = 0.4, 
  order      = TRUE, # <- replace sort.cell with order
  label     = TRUE, 
  min.cutoff = "q10"
)

# Save image
ggsave(
  filename = "5_harmonyintegratedUMAP_coloredbyinitialQCmetrics.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

