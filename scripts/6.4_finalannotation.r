# Set the working directory
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

#Load integrated object 
objs <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithtSNE_seurat_30dims.rds")

DefaultAssay(objs) # if needed: DefaultAssay(objs) <- "SCT"

#Convert to SingleCellExperiment for SingleR 
objs.sce <- as.SingleCellExperiment(objs, assay = "SCT", layers = "data")

#Reference: celldex Blueprint+ENCODE
surveyReferences()  # just to list what's available
blue_encode <- celldex::BlueprintEncodeData()

#SingleR annotation (MAIN labels)
annot_ENCODE <- SingleR(test = objs.sce, ref = blue_encode, labels = blue_encode$label.main)

#Quality heatmap (pruned labels)
#jpeg("6_harmonyintegrated_singleR_ENCODEpredictions_quality.jpg",
#     units = "in", width = 15, height = 12, res = 300)
#plotScoreHeatmap(annot_ENCODE)
#dev.off()

#add pruned labels to Seurat metadata
objs <- AddMetaData(objs, annot_ENCODE$pruned.labels, col.name = "SingleR_ENCODE")

#set1: Relabel / recode, cluster overrides, remove Neutrophils
cat("\n== BEFORE relabel ==\n")
print(sort(table(objs$SingleR_ENCODE), decreasing = TRUE))

newlab <- as.character(objs$SingleR_ENCODE)

#helper: recode a set of exact labels (case-insensitive) to a single target
recode_to <- function(targets, to) {
  idx <- which(!is.na(newlab) & tolower(newlab) %in% tolower(targets))
  if (length(idx)) newlab[idx] <<- to
}

#Epithelial cells
recode_to(
  c("Skeletal muscle", "Smooth muscle", "Adipocytes", "Keratinocytes", "HSC",
   "Erythrocytes"),    
  "Epithelial cells"
)

#Cluster 10 → Epithelial cells (force override)
if ("seurat_clusters" %in% colnames(objs@meta.data)) {
  cl10 <- which(as.character(objs$seurat_clusters) == "10")
  if (length(cl10)) newlab[cl10] <- "Epithelial cells"
}

#Fibroblasts
recode_to(
  c("Chrondrocytes", "Chondrocytes", "Mesangial cells", "Pericytes"),
  "Fibroblasts"
)

#T-cells
recode_to(
  c("CD8+ T-cells", "CD4+ T-cells"),
  "T-cells"
)

#Endothelial cells
recode_to(c("Myocytes"), "Endothelial cells")

#Cluster 20 → Endothelial cells
if ("seurat_clusters" %in% colnames(objs@meta.data)) {
  cl20 <- which(as.character(objs$seurat_clusters) == "20")
  if (length(cl20)) newlab[cl20] <- "Endothelial cells"
}

#Myeloid cells
recode_to(
  c("Monocytes", "Macrophages", "DC"),
  "Myeloid cells"
)

#Write back
objs$SingleR_ENCODE <- factor(newlab)

#Remove Neutrophils entirely (drop those cells)
if ("Neutrophils" %in% levels(objs$SingleR_ENCODE)) {
  objs <- subset(objs, subset = SingleR_ENCODE != "Neutrophils")
  objs$SingleR_ENCODE <- droplevels(objs$SingleR_ENCODE)
}

cat("\n== AFTER relabel & removal ==\n")
print(sort(table(objs$SingleR_ENCODE), decreasing = TRUE))

#REMOVE NA FROM PLOTS 
plot_obj <- subset(objs, subset = !is.na(SingleR_ENCODE))
plot_obj$SingleR_ENCODE <- droplevels(plot_obj$SingleR_ENCODE)

cat("\n== PLOT OBJECT LABELS (no NA) ==\n")
print(sort(table(plot_obj$SingleR_ENCODE), decreasing = TRUE))

#Plots with updated labels
p_tsne <- DimPlot(
  plot_obj, reduction = "tSNE.harmony", group.by = "SingleR_ENCODE",
  label = FALSE, repel = FALSE, label.size = 3
) + ggtitle("Manual + SingleR automated combined Cell Type Annotation t-SNE")
p_tsne

ggsave("6_harmonyintegratedtSNE_manual+singleR_ENCODE_final.jpg", p_tsne,
       units = "in", width = 12, height = 8, dpi = 300)

library(patchwork)

# Safety: create plot_obj if it doesn't exist
if (!exists("plot_obj")) {
  plot_obj <- subset(objs, subset = !is.na(SingleR_ENCODE))
  plot_obj$SingleR_ENCODE <- droplevels(plot_obj$SingleR_ENCODE)
}

# Use SingleR_ENCODE as identities
Idents(plot_obj) <- "SingleR_ENCODE"

cell_types <- c(
  "Epithelial cells","T-cells","B-cells",
  "Myeloid cells","NK cells","Fibroblasts","Endothelial cells"
)

# Keep only types present; skip the rest
present_types <- intersect(cell_types, levels(Idents(plot_obj)))
missing_types  <- setdiff(cell_types, present_types)
if (length(missing_types)) {
  message("Skipping missing types: ", paste(missing_types, collapse = ", "))
}

# Build per-type highlight plots on tSNE
plots_by_type <- lapply(present_types, function(ct) {
  sel <- WhichCells(plot_obj, idents = ct)
  DimPlot(
    plot_obj,
    reduction = "tSNE.harmony",
    cells.highlight = sel,
    cols = "grey85",          # non-highlighted cells
    cols.highlight = "red",   # highlighted cells
    raster = TRUE,
    pt.size = 0.2
  ) + ggtitle(ct) + NoLegend()
})

# Arrange (3 columns)
p_tsne <- wrap_plots(plots_by_type, ncol = 3)

# Show & save
p_tsne
ggsave("6_harmonyintegratedtSNE_manual+singleR_ENCODE_final_celltypehighlight.jpg", p_tsne,
       width = 12, height = 8, units = "in", dpi = 300)

#UMAPS
p_umap <- DimPlot(
  plot_obj, reduction = "umap.harmony", group.by = "SingleR_ENCODE",
  label = TRUE, repel = TRUE, label.size = 3
) + ggtitle("Manual + SingleR automated combined Cell Type Annotation UMAP")
p_umap
ggsave("6_harmonyintegratedUMAP_manual+singleR_ENCODE_final.jpg", p_umap,
       units = "in", width = 12, height = 8, dpi = 300)

# Save harmony integrated Seurat object with final cell type annotations
saveRDS(plot_obj, "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithfinalcelltypeannotations_seurat_30dims.rds")


