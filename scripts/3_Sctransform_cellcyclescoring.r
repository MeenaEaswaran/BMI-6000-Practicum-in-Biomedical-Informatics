#Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

##### Standard pre-processing workflow (includes 4 steps)
# 1. Quality control and selecting cells for further analysis
# Check QC on individual samples and then do QC on joint samples
# 2. Data normalization 
# 3. Identification of highly variable features (feature selection)
# 4. Data scaling

#2,3,4 steps to be replaced by scTransform

#Run for the first time
#BiocManager::install("glmGamPoi") #needed for sctransform faster run

#To fix errors
#install.packages("devtools")
#devtools::install_github("satijalab/seurat", ref = "fix/v.5.3.1")

# load the libraries if required
library(Seurat)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(glmGamPoi)

# Read Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes-SKIP IF NOT REQUIRED)
#All controls
patient_1_N_filtered <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_1_N_filtered.RDS")
patient_2_N_filtered <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_2_N_filtered.RDS")
patient_3_N_filtered <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_3_N_filtered.RDS")
#All cases
patient_1_T_filtered <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_1_T_filtered.RDS")
patient_2_T_filtered <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_2_T_filtered.RDS")
patient_3_T_filtered <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_3_T_filtered.RDS")

#Perform SCTransform for each dataset
#Controls-patient_1_N
patient_1_N <- SCTransform(
  patient_1_N_filtered,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb"),
  verbose = FALSE,
)
#assays should have sc transformed data and RNA. Same in meta.data

#Controls-patient_2_N
patient_2_N <- SCTransform(
  patient_2_N_filtered,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb"),
  verbose = FALSE,
)

#Controls-patient_3_N
patient_3_N <- SCTransform(
  patient_3_N_filtered,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb"),
  verbose = FALSE,
)

#Cases-patient_1_T
patient_1_T <- SCTransform(
  patient_1_T_filtered,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb"),
  verbose = FALSE,
)

#Cases-patient_2_T
patient_2_T <- SCTransform(
  patient_2_T_filtered,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb"),
  verbose = FALSE,
)

#Cases-patient_3_T
patient_3_T <- SCTransform(
  patient_3_T_filtered,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb"),
  verbose = FALSE,
)

# Cell cycle scoring on SCT assay (requires normalized data for each sample)

#Control-Patient_1_N
DefaultAssay(patient_1_N) <- "SCT"
cc <- Seurat::cc.genes.updated.2019
s.genes.p1N   <- intersect(cc$s.genes,   rownames(patient_1_N))
g2m.genes.p1N <- intersect(cc$g2m.genes, rownames(patient_1_N))
patient_1_N <- CellCycleScoring(
  patient_1_N,
  s.features   = s.genes.p1N,
  g2m.features = g2m.genes.p1N,
  set.ident    = FALSE
)

# Perform PCA and color by cell cycle phase
patient_1_N <- RunPCA(patient_1_N)
DimPlot(
  patient_1_N,
  reduction = "pca",
  group.by  = "Phase", 
  split.by = "Phase") + ggtitle("Patient_1_N_PCA after SCTransform with mt/ribo/hb regression")

#TO save image 
ggsave("3_patient_1_N_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_1_N_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

#VlnPlot(patient_1_N, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", ncol = 3, pt.size = .1)
# FeatureScatter(patient_1_N, "S.Score", "G2M.Score", group.by = "Phase")

#Re-run SCTransofrm to see any differences wuth these cell cycle score
#Switch back to the "RNA" assay, which still holds your original raw counts.
#Run SCTransform() again, this time with the extra regression variables.
#that second run overwrites/updates the "SCT" assay with the new variance-stabilized values.

DefaultAssay(patient_1_N) <- "RNA"  # SCT reads counts from RNA as it is a new run
patient_1_N <- SCTransform(
  patient_1_N,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb","S.Score","G2M.Score"),
  verbose = FALSE
)

patient_1_N <- RunPCA(patient_1_N)
DimPlot(patient_1_N,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase" )+ ggtitle("Patient_1_N_PCA after SCTransform with mt/ribo/hb/cc regression")

#TO save image 
ggsave("3_patient_1_N_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_1_N_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()
------------------------------------------------------------------------------------------------------------------------------------------------
#Control-Patient_2_N
DefaultAssay(patient_2_N) <- "SCT"
s.genes.p2N   <- intersect(cc$s.genes,   rownames(patient_2_N))
g2m.genes.p2N <- intersect(cc$g2m.genes, rownames(patient_2_N))
patient_2_N <- CellCycleScoring(
  patient_2_N,
  s.features   = s.genes.p2N,
  g2m.features = g2m.genes.p2N,
  set.ident    = FALSE
)

# Perform PCA and color by cell cycle phase
patient_2_N <- RunPCA(patient_2_N)
DimPlot(
  patient_2_N,
  reduction = "pca",
  group.by  = "Phase", 
  split.by = "Phase") + ggtitle("Patient_2_N_PCA after SCTransform with mt/ribo/hb regression")

#TO save image 
ggsave("3_patient_2_N_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_2_N_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Re-run SCTransofrm to see any differences wuth these cell cycle score
DefaultAssay(patient_2_N) <- "RNA"  # SCT reads counts from RNA
patient_2_N <- SCTransform(
  patient_2_N,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb","S.Score","G2M.Score"),
  verbose = FALSE
)

patient_2_N <- RunPCA(patient_2_N)
DimPlot(patient_2_N,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase" )+ ggtitle("Patient_2_N_PCA after SCTransform with mt/ribo/hb/cc regression")

#TO save image 
ggsave("3_patient_2_N_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_2_N_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

------------------------------------------------------------------------------------------------------------------------------------------------
#Control-Patient_3_N
DefaultAssay(patient_3_N) <- "SCT"
s.genes.p3N   <- intersect(cc$s.genes,   rownames(patient_3_N))
g2m.genes.p3N <- intersect(cc$g2m.genes, rownames(patient_3_N))
patient_3_N <- CellCycleScoring(
  patient_3_N,
  s.features   = s.genes.p3N,
  g2m.features = g2m.genes.p3N,
  set.ident    = FALSE
)

# Perform PCA and color by cell cycle phase
patient_3_N <- RunPCA(patient_3_N)
DimPlot(
  patient_3_N,
  reduction = "pca",
  group.by  = "Phase", 
  split.by = "Phase") + ggtitle("Patient_3_N_PCA after SCTransform with mt/ribo/hb regression")

#TO save image 
ggsave("3_patient_3_N_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_3_N_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Re-run SCTransofrm to see any differences wuth these cell cycle score
DefaultAssay(patient_3_N) <- "RNA"  # SCT reads counts from RNA
patient_3_N <- SCTransform(
  patient_3_N,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb","S.Score","G2M.Score"),
  verbose = FALSE
)

patient_3_N <- RunPCA(patient_3_N)
DimPlot(patient_3_N,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase" )+ ggtitle("Patient_3_N_PCA after SCTransform with mt/ribo/hb/cc regression")

#TO save image 
ggsave("3_patient_3_N_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_3_N_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()
------------------------------------------------------------------------------------------------------------------------------------------------
#Cases-Patient_1_T
DefaultAssay(patient_1_T) <- "SCT"
s.genes.p1T   <- intersect(cc$s.genes,   rownames(patient_1_T))
g2m.genes.p1T <- intersect(cc$g2m.genes, rownames(patient_1_T))
patient_1_T <- CellCycleScoring(
  patient_1_T,
  s.features   = s.genes.p1T,
  g2m.features = g2m.genes.p1T,
  set.ident    = FALSE
)

# Perform PCA and color by cell cycle phase
patient_1_T <- RunPCA(patient_1_T)
DimPlot(
  patient_1_T,
  reduction = "pca",
  group.by  = "Phase", 
  split.by = "Phase") + ggtitle("Patient_1_T_PCA after SCTransform with mt/ribo/hb regression")

#TO save image 
ggsave("3_patient_1_T_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_1_T_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Re-run SCTransofrm to see any differences wuth these cell cycle score
DefaultAssay(patient_1_T) <- "RNA"  # SCT reads counts from RNA
patient_1_T <- SCTransform(
  patient_1_T,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb","S.Score","G2M.Score"),
  verbose = FALSE
)

patient_1_T <- RunPCA(patient_1_T)
DimPlot(patient_1_T,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase" )+ ggtitle("Patient_1_T_PCA after SCTransform with mt/ribo/hb/cc regression")

#TO save image 
ggsave("3_patient_1_T_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_1_T_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

------------------------------------------------------------------------------------------------------------------------------------------------
#Cases-Patient_2_T
DefaultAssay(patient_2_T) <- "SCT"
s.genes.p2T   <- intersect(cc$s.genes,   rownames(patient_2_T))
g2m.genes.p2T <- intersect(cc$g2m.genes, rownames(patient_2_T))
patient_2_T <- CellCycleScoring(
  patient_2_T,
  s.features   = s.genes.p2T,
  g2m.features = g2m.genes.p2T,
  set.ident    = FALSE
)

# Perform PCA and color by cell cycle phase
patient_2_T <- RunPCA(patient_2_T)
DimPlot(
  patient_2_T,
  reduction = "pca",
  group.by  = "Phase", 
  split.by = "Phase") + ggtitle("Patient_2_T_PCA after SCTransform with mt/ribo/hb regression")

#TO save image 
ggsave("3_patient_2_T_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_2_T_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Re-run SCTransofrm to see any differences with these cell cycle score
DefaultAssay(patient_2_T) <- "RNA"  # SCT reads counts from RNA
patient_2_T <- SCTransform(
  patient_2_T,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb","S.Score","G2M.Score"),
  verbose = FALSE
)

patient_2_T <- RunPCA(patient_2_T)
DimPlot(patient_2_T,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase" )+ ggtitle("Patient_2_T_PCA after SCTransform with mt/ribo/hb/cc regression")

#TO save image 
ggsave("3_patient_2_T_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_2_T_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

------------------------------------------------------------------------------------------------------------------------------------------------
#Cases-Patient_3_T
DefaultAssay(patient_3_T) <- "SCT"
s.genes.p3T   <- intersect(cc$s.genes,   rownames(patient_3_T))
g2m.genes.p3T <- intersect(cc$g2m.genes, rownames(patient_3_T))
patient_3_T <- CellCycleScoring(
  patient_3_T,
  s.features   = s.genes.p3T,
  g2m.features = g2m.genes.p3T,
  set.ident    = FALSE
)

# Perform PCA and color by cell cycle phase
patient_3_T <- RunPCA(patient_3_T)
DimPlot(
  patient_3_T,
  reduction = "pca",
  group.by  = "Phase", 
  split.by = "Phase") + ggtitle("Patient_3_T_PCA after SCTransform with mt/ribo/hb regression")

#TO save image 
ggsave("3_patient_3_T_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_3_T_cellcycle_PCA_preSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Re-run SCTransofrm to see any differences with these cell cycle score
DefaultAssay(patient_3_T) <- "RNA"  # SCT reads counts from RNA
patient_3_T <- SCTransform(
  patient_3_T,
  vars.to.regress = c("percent.mt","percent.ribo","percent.hb","S.Score","G2M.Score"),
  verbose = FALSE
)

patient_3_T <- RunPCA(patient_3_T)
DimPlot(patient_3_T,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase" )+ ggtitle("Patient_3_T_PCA after SCTransform with mt/ribo/hb/cc regression")

#TO save image 
ggsave("3_patient_3_T_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, dpi=300)
jpeg("3_patient_3_T_cellcycle_PCA_postSCtransformregression.jpg", units="in", width=10, height=8, res=300)
dev.off()

------------------------------------------------------------------------------------------------------------------------------------------
# Save Seurat objects 
saveRDS(patient_1_N, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_1_N_scTransformed.RDS")
saveRDS(patient_2_N, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_2_N_scTransformed.RDS")
saveRDS(patient_3_N, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_3_N_scTransformed.RDS")
saveRDS(patient_1_T, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_1_T_scTransformed.RDS")
saveRDS(patient_2_T, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_2_T_scTransformed.RDS")
saveRDS(patient_3_T, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_3_T_scTransformed.RDS")
