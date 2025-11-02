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
objs <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithfinalcelltypeannotations_seurat_30dims.rds")

DefaultAssay(objs) # if needed: DefaultAssay(objs) <- "SCT"

#Perform SCT normalized data based DE based on Seurat recommendations 

#To run differential expression, we make use of ‘corrected counts’ that are stored in the data slot of the the SCT assay. #Corrected counts are obtained by setting the sequencing depth for all the cells to a fixed value and reversing the learned regularized negative-binomial regression model. Prior to performing differential expression, we first run PrepSCTFindMarkers, which ensures that the fixed value is set properly. 
#Then we use FindMarkers(assay="SCT") to find differentially expressed genes.
objs<-PrepSCTFindMarkers(objs,verbose=T)

#Once this is done, we can perform differential expression testing for each cluster compared to all other clusters using FindAllMarkers(). only.pos = TRUE will only return marker genes with an avg_log2FC.

#DefaultAssay(objs) <-"SCT" #reset in case changed

#All DEGs for specific cell type
# 1) Subset to the target annotation: Epithelial cells 
objs_epi <- subset(objs, subset = SingleR_ENCODE == "Epithelial cells")
Idents(objs_epi) <- objs_epi$Condition
table(Idents(objs_epi))

markers_epithelialcells <- FindMarkers(
  object          = objs_epi,
  ident.1         = "Tumor",
  ident.2         = "Normal", 
  assay = "SCT",
  slot ="data", #default
  fc.slot = "data", #default
  only.pos        = FALSE, #default # returns pos and neg markers
  min.pct         = 0.25, 
  logfc.threshold = 0, #to get all genes
  verbose = TRUE, #default
  test.use        = "wilcox",   #default  
  recorrect_umi = FALSE, #default is TRUE remove to not get error message
  base = 2, #default  
  max.cells.per.ident = Inf, #default
    )
markers_epithelialcells
# Save unfiltered DE for epithelial cells
markers_epithelialcells$gene <- rownames(markers_epithelialcells)
utils::write.csv(markers_epithelialcells, "RAW_unfiltered_DE_Tumor_vs_Normal_allepithelial_SCT.csv", row.names = FALSE)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.05 
#logFC already applied
markers_epithelialcells_de_sig <- markers_epithelialcells[markers_epithelialcells$p_val_adj <= 0.05 & abs(markers_epithelialcells$avg_log2FC) > 0.58,] 
markers_epithelialcells_de_sig

# Save final filtered and significant DE for epithelial cells
markers_epithelialcells_de_sig$gene <- rownames(markers_epithelialcells_de_sig)
utils::write.csv(markers_epithelialcells_de_sig, "SIG_filtered_DE_Tumor_vs_Normal_allepithelial_SCT.csv", row.names = FALSE)

#Visualization by dot plot or violin plot of top 20 upp or downregulated genes
#Split/subset by direction
up_epi <- subset(markers_epithelialcells_de_sig,
                           p_val_adj <= 0.05 & avg_log2FC > 0.58)

down_epi <- subset(markers_epithelialcells_de_sig,
                 p_val_adj <= 0.05 & avg_log2FC < -0.58)

#Pick top 20 genes in each set
top_up_genes_epi   <- head(rownames(up_epi),   min(20, nrow(up_epi)))
top_down_genes_epi <- head(rownames(down_epi), min(20, nrow(down_epi)))

#DotPlots
p_up_epi <- DotPlot(
  objs_epi,
  features = rev(top_up_genes_epi),   # reversed so top hits appear on top of plot
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Upregulated DEGs in Epithelial cells-Tumor vs. Normal")
p_up_epi
ggsave("7_harmonyintegrated_epithelialcells_top20upregulatedDEG_dotplot.jpg", p_up_epi,
       units = "in", width = 12, height = 8, dpi = 300)

p_down_epi <- DotPlot(
  objs_epi,
  features = rev(top_down_genes_epi),
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Downregulated DEGs in Epithelial cells-Tumor vs. Normal")
p_down_epi
ggsave("7_harmonyintegrated_epithelialcells_top20downregulatedDEG_dotplot.jpg", p_down_epi,
       units = "in", width = 12, height = 8, dpi = 300)

#Violin plot
#Upregulated in Tumor (top 20)
p_up_vln_epi <- VlnPlot(
  objs_epi,
  features = top_up_genes_epi,
  group.by = "Condition",   # x-axis: Normal vs Tumor
  assay    = "SCT",
  pt.size  = 0.1
)
p_up_vln_epi
ggsave("7_harmonyintegrated_epithelialcells_top20upregulatedDEG_violinplot.jpg", p_up_vln_epi,
       units = "in", width = 12, height = 12, dpi = 300)

#Downregulated in Tumor (i.e., Normal-up; top 20)
p_down_vln_epi <- VlnPlot(
  objs_epi,
  features = top_down_genes_epi,
  group.by = "Condition",
  assay    = "SCT",
  pt.size  = 0.1
)
p_down_vln_epi
ggsave("7_harmonyintegrated_epithelialcells_top20downpregulatedDEG_violinplot.jpg", p_down_vln_epi,
       units = "in", width = 12, height = 12, dpi = 300)
######################################################################################################################

# 2) Subset to the target annotation: Fibroblasts
objs_fibro <- subset(objs, subset = SingleR_ENCODE == "Fibroblasts")
Idents(objs_fibro) <- objs_fibro$Condition
table(Idents(objs_fibro))

markers_fibro <- FindMarkers(
  object          = objs_fibro,
  ident.1         = "Tumor",
  ident.2         = "Normal", 
  assay = "SCT",
  slot ="data", #default
  fc.slot = "data", #default
  only.pos        = FALSE, #default # returns pos and neg markers
  min.pct         = 0.25, 
  logfc.threshold = 0, #to get all genes
  verbose = TRUE, #default
  test.use        = "wilcox",   #default  
  recorrect_umi = FALSE, #default is TRUE remove to not get error message
  base = 2, #default  
  max.cells.per.ident = Inf, #default
)
markers_fibro

# Save unfiltered DE for fibroblasts
markers_fibro$gene <- rownames(markers_fibro)
utils::write.csv(markers_fibro, "RAW_unfiltered_DE_Tumor_vs_Normal_allfibroblasts_SCT.csv", row.names = FALSE)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.05 
#logFC already applied
# Keep genes with adjusted p < 0.05 and absolute log2FC > 0.58
markers_fibro_de_sig <- markers_fibro[
  markers_fibro$p_val_adj <= 0.05 & abs(markers_fibro$avg_log2FC) > 0.58,
]
markers_fibro_de_sig

# Save final filtered and significant DE for fibroblasts
markers_fibro_de_sig$gene <- rownames(markers_fibro_de_sig)
utils::write.csv(markers_fibro_de_sig, "SIG_filtered_DE_Tumor_vs_Normal_allfibroblasts_SCT.csv", row.names = FALSE)

#Visualization by dot plot or violin plot of top 20 up or downregulated genes
#Split/subset by direction
up_fibro <- subset(markers_fibro_de_sig,
                 p_val_adj <= 0.05 & avg_log2FC > 0.58)

down_fibro <- subset(markers_fibro_de_sig,
                   p_val_adj <= 0.05 & avg_log2FC < -0.58)

#Pick top 20 genes in each set
top_up_genes_fibro   <- head(rownames(up_fibro),   min(20, nrow(up_fibro)))
top_down_genes_fibro <- head(rownames(down_fibro), min(20, nrow(down_fibro)))

#DotPlots
p_up_fibro <- DotPlot(
  objs_fibro,
  features = rev(top_up_genes_fibro),   # reversed so top hits appear on top of plot
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Upregulated DEGs in Fibroblasts-Tumor vs. Normal")
p_up_fibro
ggsave("7_harmonyintegrated_fibroblasts_top20upregulatedDEG_dotplot.jpg", p_up_fibro,
       units = "in", width = 12, height = 8, dpi = 300)

p_down_fibro <- DotPlot(
  objs_fibro,
  features = rev(top_down_genes_fibro),
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Downregulated DEGs in Fibroblasts-Tumor vs. Normal")
p_down_fibro
ggsave("7_harmonyintegrated_fibroblasts_top20downregulatedDEG_dotplot.jpg", p_down_fibro,
       units = "in", width = 12, height = 8, dpi = 300)

#Violin plot
#Upregulated in Tumor (top 20)
p_up_vln_fibro <- VlnPlot(
  objs_fibro,
  features = top_up_genes_fibro,
  group.by = "Condition",   # x-axis: Normal vs Tumor
  assay    = "SCT",
  pt.size  = 0.1
)
p_up_vln_fibro
ggsave("7_harmonyintegrated_fibroblasts_top20upregulatedDEG_violinplot.jpg", p_up_vln_fibro,
       units = "in", width = 12, height = 12, dpi = 300)

#Downregulated in Tumor (i.e., Normal-up; top 20)
p_down_vln_fibro <- VlnPlot(
  objs_fibro,
  features = top_down_genes_fibro,
  group.by = "Condition",
  assay    = "SCT",
  pt.size  = 0.1
)
p_down_vln_fibro
ggsave("7_harmonyintegrated_fibroblasts_top20downpregulatedDEG_violinplot.jpg", p_down_vln_fibro,
       units = "in", width = 12, height = 12, dpi = 300)

######################################################################################################################

# 3) Subset to the target annotation: Myeloid cells
objs_myl <- subset(objs, subset = SingleR_ENCODE == "Myeloid cells")
Idents(objs_myl) <- objs_myl$Condition
table(Idents(objs_myl))

markers_myl <- FindMarkers(
  object          = objs_myl,
  ident.1         = "Tumor",
  ident.2         = "Normal", 
  assay = "SCT",
  slot ="data", #default
  fc.slot = "data", #default
  only.pos        = FALSE, #default # returns pos and neg markers
  min.pct         = 0.25, 
  logfc.threshold = 0, #to get all genes
  verbose = TRUE, #default
  test.use        = "wilcox",   #default  
  recorrect_umi = FALSE, #default is TRUE remove to not get error message
  base = 2, #default  
  max.cells.per.ident = Inf, #default
)
markers_myl
# Save unfiltered DE for Myeloid cells
markers_myl$gene <- rownames(markers_myl)
utils::write.csv(markers_myl, "RAW_unfiltered_DE_Tumor_vs_Normal_allmyeloidcells_SCT.csv", row.names = FALSE)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.05 
#logFC already applied
markers_myl_de_sig <- markers_myl[markers_myl$p_val_adj <= 0.05 & abs(markers_myl$avg_log2FC) > 0.58,
]
markers_myl_de_sig

# Save final filtered and significant DE for Myeloid cells
markers_myl_de_sig$gene <- rownames(markers_myl_de_sig)
utils::write.csv(markers_myl_de_sig, "SIG_filtered_DE_Tumor_vs_Normal_allmyeloidcells_SCT.csv", row.names = FALSE)

#Visualization by dot plot or violin plot of top 20 up or downregulated genes
#Split/subset by direction
up_myl <- subset(markers_myl_de_sig,
                   p_val_adj <= 0.05 & avg_log2FC > 0.58)

down_myl <- subset(markers_myl_de_sig,
                     p_val_adj <= 0.05 & avg_log2FC < -0.58)

#Pick top 20 genes in each set
top_up_genes_myl   <- head(rownames(up_myl),   min(20, nrow(up_myl)))
top_down_genes_myl <- head(rownames(down_myl), min(20, nrow(down_myl)))

#DotPlots
p_up_myl <- DotPlot(
  objs_myl,
  features = rev(top_up_genes_myl),   # reversed so top hits appear on top of plot
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Upregulated DEGs in Myeloid cells-Tumor vs. Normal")
p_up_myl
ggsave("7_harmonyintegrated_myeloidcells_top20upregulatedDEG_dotplot.jpg", p_up_myl,
       units = "in", width = 12, height = 8, dpi = 300)

p_down_myl <- DotPlot(
  objs_myl,
  features = rev(top_down_genes_myl),
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Downregulated DEGs in Myeloid cells-Tumor vs. Normal")
p_down_myl
ggsave("7_harmonyintegrated_myeloidcells_top20downregulatedDEG_dotplot.jpg", p_down_myl,
       units = "in", width = 12, height = 8, dpi = 300)

#Violin plot
#Upregulated in Tumor (top 20)
p_up_vln_myl <- VlnPlot(
  objs_myl,
  features = top_up_genes_myl,
  group.by = "Condition",   # x-axis: Normal vs Tumor
  assay    = "SCT",
  pt.size  = 0.1
)
p_up_vln_myl
ggsave("7_harmonyintegrated_myeloidcells_top20upregulatedDEG_violinplot.jpg", p_up_vln_myl,
       units = "in", width = 12, height = 12, dpi = 300)

#Downregulated in Tumor (i.e., Normal-up; top 20)
p_down_vln_myl <- VlnPlot(
  objs_myl,
  features = top_down_genes_myl,
  group.by = "Condition",
  assay    = "SCT",
  pt.size  = 0.1
)
p_down_vln_myl
ggsave("7_harmonyintegrated_myeloidcells_top20downpregulatedDEG_violinplot.jpg", p_down_vln_myl,
       units = "in", width = 12, height = 12, dpi = 300)

######################################################################################################################

# 4) Subset to the target annotation: NK cells
objs_NK <- subset(objs, subset = SingleR_ENCODE == "NK cells")
Idents(objs_NK) <- objs_NK$Condition
table(Idents(objs_NK))

markers_NK <- FindMarkers(
  object          = objs_NK,
  ident.1         = "Tumor",
  ident.2         = "Normal", 
  assay = "SCT",
  slot ="data", #default
  fc.slot = "data", #default
  only.pos        = FALSE, #default # returns pos and neg markers
  min.pct         = 0.25, 
  logfc.threshold = 0, #to get all genes
  verbose = TRUE, #default
  test.use        = "wilcox",   #default  
  recorrect_umi = FALSE, #default is TRUE remove to not get error message
  base = 2, #default  
  max.cells.per.ident = Inf, #default
)
markers_NK
# Save unfiltered DE for NK cells
markers_NK$gene <- rownames(markers_NK)
utils::write.csv(markers_NK, "RAW_unfiltered_DE_Tumor_vs_Normal_allNKcells_SCT.csv", row.names = FALSE)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.05 
#logFC already applied
markers_NK_de_sig <- markers_NK[markers_NK$p_val_adj <= 0.05 & abs(markers_NK$avg_log2FC) > 0.58,
]
markers_NK_de_sig

# Save final filtered and significant DE for NK cells
markers_NK_de_sig$gene <- rownames(markers_NK_de_sig)
utils::write.csv(markers_NK_de_sig, "SIG_filtered_DE_Tumor_vs_Normal_allNKcells_SCT.csv", row.names = FALSE)

#Visualization by dot plot or violin plot of top 20 up or downregulated genes
#Split/subset by direction
up_NK <- subset(markers_NK_de_sig,
                 p_val_adj <= 0.05 & avg_log2FC > 0.58)

down_NK <- subset(markers_NK_de_sig,
                   p_val_adj <= 0.05 & avg_log2FC < -0.58)

#Pick top 20 genes in each set
top_up_genes_NK   <- head(rownames(up_NK),   min(20, nrow(up_NK)))
top_down_genes_NK <- head(rownames(down_NK), min(20, nrow(down_NK)))

#DotPlots
p_up_NK <- DotPlot(
  objs_NK,
  features = rev(top_up_genes_NK),   # reversed so top hits appear on top of plot
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Upregulated DEGs in NK cells-Tumor vs. Normal")
p_up_NK
ggsave("7_harmonyintegrated_NKcells_top20upregulatedDEG_dotplot.jpg", p_up_NK,
       units = "in", width = 12, height = 8, dpi = 300)

p_down_NK <- DotPlot(
  objs_NK,
  features = rev(top_down_genes_NK),
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Downregulated DEGs in NK cells-Tumor vs. Normal")
p_down_NK
ggsave("7_harmonyintegrated_NKcells_top20downregulatedDEG_dotplot.jpg", p_down_NK,
       units = "in", width = 12, height = 8, dpi = 300)

#Violin plot
#Upregulated in Tumor (top 20)
p_up_vln_NK <- VlnPlot(
  objs_NK,
  features = top_up_genes_NK,
  group.by = "Condition",   # x-axis: Normal vs Tumor
  assay    = "SCT",
  pt.size  = 0.1
)
p_up_vln_NK
ggsave("7_harmonyintegrated_NKcells_top20upregulatedDEG_violinplot.jpg", p_up_vln_NK,
       units = "in", width = 12, height = 12, dpi = 300)

#Downregulated in Tumor (i.e., Normal-up; top 20)
p_down_vln_NK <- VlnPlot(
  objs_NK,
  features = top_down_genes_NK,
  group.by = "Condition",
  assay    = "SCT",
  pt.size  = 0.1
)
p_down_vln_NK
ggsave("7_harmonyintegrated_NKcells_top20downpregulatedDEG_violinplot.jpg", p_down_vln_NK,
       units = "in", width = 12, height = 12, dpi = 300)

######################################################################################################################

# 5) Subset to the target annotation: B-cells
objs_b <- subset(objs, subset = SingleR_ENCODE == "B-cells")
Idents(objs_b) <- objs_b$Condition
table(Idents(objs_b))

markers_b <- FindMarkers(
  object          = objs_b,
  ident.1         = "Tumor",
  ident.2         = "Normal", 
  assay = "SCT",
  slot ="data", #default
  fc.slot = "data", #default
  only.pos        = FALSE, #default # returns pos and neg markers
  min.pct         = 0.25, 
  logfc.threshold = 0, #to get all genes
  verbose = TRUE, #default
  test.use        = "wilcox",   #default  
  recorrect_umi = FALSE, #default is TRUE remove to not get error message
  base = 2, #default  
  max.cells.per.ident = Inf, #default
)
markers_b
# Save unfiltered DE for B-cells
markers_b$gene <- rownames(markers_b)
utils::write.csv(markers_b, "RAW_unfiltered_DE_Tumor_vs_Normal_allBcells_SCT.csv", row.names = FALSE)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.05 
#logFC already applied
markers_b_de_sig <- markers_b[markers_b$p_val_adj <= 0.05 & abs(markers_b$avg_log2FC) > 0.58,
]
markers_b_de_sig

# Save final filtered and significant DE for B-cells
markers_b_de_sig$gene <- rownames(markers_b_de_sig)
utils::write.csv(markers_b_de_sig, "SIG_filtered_DE_Tumor_vs_Normal_allBcells_SCT.csv", row.names = FALSE)

#Visualization by dot plot or violin plot of top 20 up or downregulated genes
#Split/subset by direction
up_b <- subset(markers_b_de_sig,
                p_val_adj <= 0.05 & avg_log2FC > 0.58)

down_b <- subset(markers_b_de_sig,
                  p_val_adj <= 0.05 & avg_log2FC < -0.58)

#Pick top 20 genes in each set
top_up_genes_b   <- head(rownames(up_b),   min(20, nrow(up_b)))
top_down_genes_b <- head(rownames(down_b), min(20, nrow(down_b)))

#DotPlots
p_up_b <- DotPlot(
  objs_b,
  features = rev(top_up_genes_b),   # reversed so top hits appear on top of plot
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Upregulated DEGs in B-cells-Tumor vs. Normal")
p_up_b
ggsave("7_harmonyintegrated_Bcells_top20upregulatedDEG_dotplot.jpg", p_up_b,
       units = "in", width = 12, height = 8, dpi = 300)

p_down_b <- DotPlot(
  objs_b,
  features = rev(top_down_genes_b),
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Downregulated DEGs in B-cells-Tumor vs. Normal")
p_down_b
ggsave("7_harmonyintegrated_Bcells_top20downregulatedDEG_dotplot.jpg", p_down_b,
       units = "in", width = 12, height = 8, dpi = 300)

#Violin plot
#Upregulated in Tumor (top 20)
p_up_vln_b <- VlnPlot(
  objs_b,
  features = top_up_genes_b,
  group.by = "Condition",   # x-axis: Normal vs Tumor
  assay    = "SCT",
  pt.size  = 0.1
)
p_up_vln_b
ggsave("7_harmonyintegrated_Bcells_top20upregulatedDEG_violinplot.jpg", p_up_vln_b,
       units = "in", width = 12, height = 12, dpi = 300)

#Downregulated in Tumor (i.e., Normal-up; top 20)
p_down_vln_b <- VlnPlot(
  objs_b,
  features = top_down_genes_b,
  group.by = "Condition",
  assay    = "SCT",
  pt.size  = 0.1
)
p_down_vln_b
ggsave("7_harmonyintegrated_Bcells_top20downpregulatedDEG_violinplot.jpg", p_down_vln_b,
       units = "in", width = 12, height = 12, dpi = 300)

######################################################################################################################

# 6) Subset to the target annotation: T-cells
objs_t <- subset(objs, subset = SingleR_ENCODE == "T-cells")
Idents(objs_t) <- objs_t$Condition
table(Idents(objs_t))

markers_t <- FindMarkers(
  object          = objs_t,
  ident.1         = "Tumor",
  ident.2         = "Normal", 
  assay = "SCT",
  slot ="data", #default
  fc.slot = "data", #default
  only.pos        = FALSE, #default # returns pos and neg markers
  min.pct         = 0.25, 
  logfc.threshold = 0, #to get all genes
  verbose = TRUE, #default
  test.use        = "wilcox",   #default  
  recorrect_umi = FALSE, #default is TRUE remove to not get error message
  base = 2, #default  
  max.cells.per.ident = Inf, #default
)
markers_t
# Save unfiltered DE for T-cells
markers_t$gene <- rownames(markers_t)
utils::write.csv(markers_t, "RAW_unfiltered_DE_Tumor_vs_Normal_allTcells_SCT.csv", row.names = FALSE)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.05 
#logFC already applied
markers_t_de_sig <- markers_t[markers_t$p_val_adj <= 0.05 & abs(markers_t$avg_log2FC) > 0.58,
]
markers_t_de_sig

# Save final filtered and significant DE for T-cells
markers_t_de_sig$gene <- rownames(markers_t_de_sig)
utils::write.csv(markers_t_de_sig, "SIG_filtered_DE_Tumor_vs_Normal_allTcells_SCT.csv", row.names = FALSE)

#Visualization by dot plot or violin plot of top 20 up or downregulated genes
#Split/subset by direction
up_t <- subset(markers_t_de_sig,
               p_val_adj <= 0.05 & avg_log2FC > 0.58)

down_t <- subset(markers_t_de_sig,
                 p_val_adj <= 0.05 & avg_log2FC < -0.58)

#Pick top 20 genes in each set
top_up_genes_t   <- head(rownames(up_t),   min(20, nrow(up_t)))
top_down_genes_t <- head(rownames(down_t), min(20, nrow(down_t)))

#DotPlots
p_up_t <- DotPlot(
  objs_t,
  features = rev(top_up_genes_t),   # reversed so top hits appear on top of plot
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Upregulated DEGs in T-cells-Tumor vs. Normal")
p_up_t
ggsave("7_harmonyintegrated_Tcells_top20upregulatedDEG_dotplot.jpg", p_up_t,
       units = "in", width = 12, height = 8, dpi = 300)

p_down_t <- DotPlot(
  objs_t,
  features = rev(top_down_genes_t),
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Downregulated DEGs in T-cells-Tumor vs. Normal")
p_down_t
ggsave("7_harmonyintegrated_Tcells_top20downregulatedDEG_dotplot.jpg", p_down_t,
       units = "in", width = 12, height = 8, dpi = 300)

#Violin plot
#Upregulated in Tumor (top 20)
p_up_vln_t <- VlnPlot(
  objs_t,
  features = top_up_genes_t,
  group.by = "Condition",   # x-axis: Normal vs Tumor
  assay    = "SCT",
  pt.size  = 0.1
)
p_up_vln_t
ggsave("7_harmonyintegrated_Tcells_top20upregulatedDEG_violinplot.jpg", p_up_vln_t,
       units = "in", width = 12, height = 12, dpi = 300)

#Downregulated in Tumor (i.e., Normal-up; top 20)
p_down_vln_t <- VlnPlot(
  objs_t,
  features = top_down_genes_t,
  group.by = "Condition",
  assay    = "SCT",
  pt.size  = 0.1
)
p_down_vln_t
ggsave("7_harmonyintegrated_Tcells_top20downpregulatedDEG_violinplot.jpg", p_down_vln_t,
       units = "in", width = 12, height = 12, dpi = 300)

######################################################################################################################

# 7) Subset to the target annotation: Endothelial cells
objs_endo <- subset(objs, subset = SingleR_ENCODE == "Endothelial cells")
Idents(objs_endo) <- objs_endo$Condition
table(Idents(objs_endo))

markers_endo <- FindMarkers(
  object          = objs_endo,
  ident.1         = "Tumor",
  ident.2         = "Normal", 
  assay = "SCT",
  slot ="data", #default
  fc.slot = "data", #default
  only.pos        = FALSE, #default # returns pos and neg markers
  min.pct         = 0.25, 
  logfc.threshold = 0, #to get all genes
  verbose = TRUE, #default
  test.use        = "wilcox",   #default  
  recorrect_umi = FALSE, #default is TRUE remove to not get error message
  base = 2, #default  
  max.cells.per.ident = Inf, #default
)
markers_endo
# Save unfiltered DE for Endothelial cells
markers_endo$gene <- rownames(markers_endo)
utils::write.csv(markers_endo, "RAW_unfiltered_DE_Tumor_vs_Normal_allendothelialcells_SCT.csv", row.names = FALSE)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.05 
#logFC already applied
markers_endo_de_sig <- markers_endo[markers_endo$p_val_adj <= 0.05 & abs(markers_endo$avg_log2FC) > 0.58,
]
markers_endo_de_sig

# Save final filtered and significant DE for Endothelial cells
markers_endo_de_sig$gene <- rownames(markers_endo_de_sig)
utils::write.csv(markers_endo_de_sig, "SIG_filtered_DE_Tumor_vs_Normal_allendothelialcells_SCT.csv", row.names = FALSE)

#Visualization by dot plot or violin plot of top 20 up or downregulated genes
#Split/subset by direction
up_endo <- subset(markers_endo_de_sig,
               p_val_adj <= 0.05 & avg_log2FC > 0.58)

down_endo <- subset(markers_endo_de_sig,
                 p_val_adj <= 0.05 & avg_log2FC < -0.58)

#Pick top 20 genes in each set
top_up_genes_endo   <- head(rownames(up_endo),   min(20, nrow(up_endo)))
top_down_genes_endo <- head(rownames(down_endo), min(20, nrow(down_endo)))

#DotPlots
p_up_endo <- DotPlot(
  objs_endo,
  features = rev(top_up_genes_endo),   # reversed so top hits appear on top of plot
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Upregulated DEGs in Endothelial cells-Tumor vs. Normal")
p_up_endo
ggsave("7_harmonyintegrated_endothelialcells_top20upregulatedDEG_dotplot.jpg", p_up_endo,
       units = "in", width = 12, height = 8, dpi = 300)

p_down_endo <- DotPlot(
  objs_endo,
  features = rev(top_down_genes_endo),
  group.by = "Condition",
  assay    = "SCT"
) + RotatedAxis() +
  ggtitle("Top 20 Downregulated DEGs in Endothelial cells-Tumor vs. Normal")
p_down_endo
ggsave("7_harmonyintegrated_endothelialcells_top20downregulatedDEG_dotplot.jpg", p_down_endo,
       units = "in", width = 12, height = 8, dpi = 300)

#Violin plot
#Upregulated in Tumor (top 20)
p_up_vln_endo <- VlnPlot(
  objs_endo,
  features = top_up_genes_endo,
  group.by = "Condition",   # x-axis: Normal vs Tumor
  assay    = "SCT",
  pt.size  = 0.1
)
p_up_vln_endo
ggsave("7_harmonyintegrated_endothelialcells_top20upregulatedDEG_violinplot.jpg", p_up_vln_endo,
       units = "in", width = 12, height = 12, dpi = 300)

#Downregulated in Tumor (i.e., Normal-up; top 20)
p_down_vln_endo <- VlnPlot(
  objs_endo,
  features = top_down_genes_endo,
  group.by = "Condition",
  assay    = "SCT",
  pt.size  = 0.1
)
p_down_vln_endo
ggsave("7_harmonyintegrated_endothelialcells_top20downpregulatedDEG_violinplot.jpg", p_down_vln_endo,
       units = "in", width = 12, height = 12, dpi = 300)
