#Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

#Run for the first time
#BiocManager::install("glmGamPoi") #needed for sctransform faster run

#To fix errors
#install.packages("devtools")
#devtools::install_github("satijalab/seurat", ref = "fix/v.5.3.1")

#install UMAP if not already installed
#reticulate::py_install(packages='umap-learn')  

# load the libraries if required
library(Seurat)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
options(future.globals.maxSize = 1e9)
#library(DoubletFinder)

# Read Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
objs<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_all_singlet_objects_combined.rds")

#Check default assay
DefaultAssay(object = objs) 

options(future.globals.maxSize = Inf)
objs <- SCTransform(objs)
objs <- RunPCA(objs, npcs = 30, verbose = F)
objs <- IntegrateLayers(
  object = objs,
  method = CCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
## IMPORTANT ### PRIOR TO RUNNING FOR CASES PLEASE READ BELOW
# SINCE THE VECTOR MEMORY EXHAUST LIMIT REACHES WHEN RUNNING PCA FOR TUMOR SAMPLES, RSTUDIO LIMIT HAS TO BE RESET USING TERMINAL IN RSTUDIO AND THEN APPLICATION HAS TO BE RESTARTED. FOLLOW THE CODE/INSTRUCTIONS ON STACKOVERFLOW AS SHOWN BELOW

#https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos

library(usethis) 
usethis::edit_r_environ()
#when the tab opens up in R studio, add this to the 1st line: R_MAX_VSIZE=100Gb (or whatever memory you wish to allocate).
#Re-start R and/or restart computer and run the R command again that gave you the memory error.
#RERUN ABOVE COMMANDS after loading relevant libraries.

#Re-run LDA/PCA/UMAPs on the CCA integrated space

#Should lready have 'integrated.dr' from CCAIntegration
Reductions(objs)    # should include "integrated.dr"

# Rebuild graph + clusters + UMAP from the CCA-integrated reduction
objs <- FindNeighbors(objs, reduction = "integrated.dr", dims = 1:30)
objs <- FindClusters(objs, resolution = 0.6)

# Give the UMAP a distinct name
objs <- RunUMAP(
  objs,
  reduction      = "integrated.dr",
  dims           = 1:30,
  reduction.name = "umap.cca"   
)

# UMAP colored by clusters / by sample (uses the CCA UMAP)
p_clusters <- DimPlot(
  objs, reduction = "umap.cca",
  group.by = "seurat_clusters", label = TRUE, repel = TRUE
) + ggtitle("Clusters (CCA integrated)")

p_samples <- DimPlot(
  objs, reduction = "umap.cca",
  group.by = "orig.ident"
) + ggtitle("By patient samples")

p_clusters + p_samples

# Save image
ggsave(
  filename = "5_CCAintegratedUMAP_coloredbyClusters_originalsamples.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

# Add tidy metadata (Condition / Patient / Sample) if not already present
objs$Condition <- ifelse(grepl("_N$", objs$orig.ident), "Normal", "Tumor")
objs$Patient   <- sub("^patient_(\\d+)_.*$", "P\\1", objs$orig.ident)
objs$Sample    <- objs$orig.ident

# Split UMAP by condition 
p_split <- DimPlot(
  objs,
  reduction = "umap.cca",
  group.by  = "seurat_clusters", label = FALSE, repel = TRUE, #to have no cluster labels on the plots
  split.by  = "Condition", ncol = 2
) + ggtitle("CCA integrated: clusters split by condition")
p_split

# Save image
ggsave(
  filename = "5_CCAintegratedUMAP_coloredbycondition_acrossallsamples.jpg",
  units = "in", width = 12, height = 8, dpi = 300
)

# color by Sample on the same CCA UMAP
DimPlot(objs, reduction = "umap.cca", group.by = "Sample") +
  ggtitle("UMAP colored by Sample (CCA integrated)")

# Save integrated seurat object
saveRDS(objs, "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_CCA_integrated_seurat_30dims.rds")
