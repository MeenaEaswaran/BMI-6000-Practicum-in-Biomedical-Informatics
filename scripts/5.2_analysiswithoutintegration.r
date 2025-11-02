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

## IMPORTANT ### PRIOR TO RUNNING FOR CASES PLEASE READ BELOW
# SINCE THE VECTOR MEMORY EXHAUST LIMIT REACHES WHEN RUNNING PCA FOR TUMOR SAMPLES, RSTUDIO LIMIT HAS TO BE RESET USING TERMINAL IN RSTUDIO AND THEN APPLICATION HAS TO BE RESTARTED. FOLLOW THE CODE/INSTRUCTIONS ON STACKOVERFLOW AS SHOWN BELOW
#https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos
library(usethis) 
usethis::edit_r_environ()
#when the tab opens up in R studio, add this to the 1st line: R_MAX_VSIZE=100Gb (or whatever memory you wish to allocate).
#Re-start R and/or restart computer and run the R command again that gave you the memory error.
#RERUN ABOVE COMMANDS after loading relevant libraries.

objs <- FindNeighbors(objs, dims = 1:30, reduction = "pca")
objs <- FindClusters(objs, resolution = 2, cluster.name = "unintegrated_clusters")
objs <- RunUMAP(objs, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


# UMAP colored by clusters
p_clusters <- DimPlot(objs, reduction = "umap.unintegrated", group.by = "unintegrated_clusters",
                      label = TRUE, repel = TRUE) + ggtitle("Clusters unintegrated")

# UMAP colored by original sample (patient_X_{N|T})
p_samples  <- DimPlot(objs, reduction = "umap.unintegrated", group.by = "orig.ident") +
  ggtitle("By patient samples unintegrated")
p_clusters + p_samples

#TO save image 
ggsave("5_unitegrated_coloredbyClusters_originalsamples.jpg", units="in", width=12, height=8, dpi=300)
jpeg("5_unitegrated_coloredbyClusters_originalsamples.jpg", units="in", width=12, height=8, res=300)
dev.off()

#Add tidy metadata: Condition & Patient
# From your orig.ident factor: patient_1_N, patient_1_T, ...
objs$Condition <- ifelse(grepl("_N$", objs$orig.ident), "Normal", "Tumor")
objs$Patient   <- sub("^patient_(\\d+)_.*$", "P\\1", objs$orig.ident)
objs$Sample    <- objs$orig.ident  # keep a clean alias


#Split UMAP by condition (Normal vs Tumor)
# Clusters in Normal vs Tumor panels
DimPlot(objs, reduction = "umap.unintegrated",
        group.by = "unintegrated_clusters", label = FALSE, repel = TRUE, #change label to False if no cluster numbers
        split.by = "Condition", ncol = 2) +
  ggtitle("Unintegrated Clusters split by Condition")
#TO save image 
ggsave("5_unintegratedUMAP_coloredbycondition_acrossallsamples.jpg", units="in", width=12, height=8, dpi=300)
jpeg("5_unintegratedUMAP_coloredbycondition_acrossallsamples.jpg", units="in", width=12, height=8, res=300)
dev.off()

# Save integrated seurat object
saveRDS(objs, "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_unintegrated_seurat_30dims.rds")
