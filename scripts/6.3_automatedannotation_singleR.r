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

#blue encode first
blue_encode <- celldex::BlueprintEncodeData() #blueprint and ENCODE. #use this as it as stroma and immune cells
unique(blue_encode$label.main)
table(blue_encode$label.main)

unique(blue_encode$label.fine)
table(blue_encode$label.fine)

#For many of the celldex datasets, there are two separate categories of labels. The main labels capture the family of cell types available (e.g. T cells), while the fine category captures the specific subtypes (e.g. Cd4 effector T cells). The references in celldex are species-specific, which has to be accounted for in any reference-based cell type annotation.

#cell-level cell type annotation using main level or fine level if needed
annot_ENCODE= SingleR(test=objs.sce, ref=blue_encode, labels=blue_encode$label.main)
head(annot_ENCODE)

#https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/Seurat_DifferentialExpression_Classification/#4-cell-annotation

table(annot_ENCODE$pruned.labels, useNA="ifany") #useNA can be used turned on in the `table` function

# Inspect quality of the predictions
# Save as high-resolution JPEG
jpeg("6_harmonyintegrated_singleR_ENCODEpredictions_quality.jpg",
     units = "in", width = 15, height = 12, res = 300)
plotScoreHeatmap(annot_ENCODE)
dev.off()

#plotDeltaDistribution(annot_ENCODE, ncol = 4, dots.on.top = FALSE)

# Add to Seurat object
objs <- AddMetaData(objs, annot_ENCODE$pruned.labels, col.name = 'SingleR_ENCODE')

# Visualize them on the tSNE
DimPlot(objs, reduction= "tSNE.harmony", group.by = "SingleR_ENCODE", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedtSNE_singleRautomatedcellannotations_ENCODE.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedtSNE_singleRautomatedcellannotations_ENCODE.jpg", units="in", width=12, height=8, res=300)
dev.off()

# Visualize them on the UMAP
DimPlot(objs, reduction= "umap.harmony", group.by = "SingleR_ENCODE", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedUMAP_singleRautomatedcellannotations_ENCODE.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedUMAP_singleRautomatedcellannotations_ENCODE.jpg", units="in", width=12, height=8, res=300)
dev.off()
############################################################################################################################
#DICE (not apt)
dice <- celldex::DatabaseImmuneCellExpressionData() 
unique(dice$label.main)
table(dice$label.main)
unique(dice$label.fine)
table(dice$label.fine)

#cell-level cell type annotation using main level or fine level if needed
annot_dice= SingleR(test=objs.sce, ref=dice, labels=dice$label.main)
head(annot_dice)
table(annot_dice$pruned.labels, useNA="ifany") #useNA can be used turned on in the `table` function

# Inspect quality of the predictions
plotScoreHeatmap(annot_dice)
plotDeltaDistribution(annot_dice, ncol = 4, dots.on.top = FALSE)

# Add to seurat object
objs <- AddMetaData(objs, annot_dice$pruned.labels, col.name = 'SingleR_DICE')

# Visualise them on the t-SNE
DimPlot(objs, reduction= "tSNE.harmony", group.by = "SingleR_DICE", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedtSNE_singleRautomatedcellannotations_DICE.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedtSNE_singleRautomatedcellannotations_DICE.jpg", units="in", width=12, height=8, res=300)
dev.off()

# Visualise them on the UMAP
DimPlot(objs, reduction= "umap.harmony", group.by = "SingleR_DICE", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedUMAP_singleRautomatedcellannotations_DICE.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedUMAP_singleRautomatedcellannotations_DICE.jpg", units="in", width=12, height=8, res=300)
dev.off()
############################################################################################################################
#Monacco (not apt)
monacco <- celldex::MonacoImmuneData()
unique(monacco$label.main)
table(monacco$label.main)
unique(monacco$label.fine)
table(monacco$label.fine)

#cell-level cell type annotation using main level or fine level if needed
annot_monacco= SingleR(test=objs.sce, ref=monacco, labels=monacco$label.main)
head(annot_monacco)

table(annot_monacco$pruned.labels, useNA="ifany") #useNA can be used turned on in the `table` function

# Inspect quality of the predictions
plotScoreHeatmap(annot_monacco)
plotDeltaDistribution(annot_monacco, ncol = 4, dots.on.top = FALSE)

# Add to seurat object
objs <- AddMetaData(objs, annot_monacco$pruned.labels, col.name = 'SingleR_monacco')

# Visualise them on the UMAP
DimPlot(objs, reduction= "umap.harmony", group.by = "SingleR_monacco", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedUMAP_singleRautomatedcellannotations_monacco.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedUMAP_singleRautomatedcellannotations_monacco.jpg", units="in", width=12, height=8, res=300)
dev.off()

# Visualise them on the tSNE
DimPlot(objs, reduction= "tSNE.harmony", group.by = "SingleR_monacco", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedtSNE_singleRautomatedcellannotations_monacco.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedtSNE_singleRautomatedcellannotations_monacco.jpg", units="in", width=12, height=8, res=300)
dev.off()
############################################################################################################################
#Human primary cell atlas (HPCA) -works
HPCA <- celldex::HumanPrimaryCellAtlasData()
unique(HPCA$label.main)
table(HPCA$label.main)
unique(HPCA$label.fine)
table(HPCA$label.fine)

#cell-level cell type annotation using main level or fine level
annot_HPCA= SingleR(test=objs.sce, ref=HPCA, labels=HPCA$label.main)
head(annot_HPCA)

table(annot_HPCA$pruned.labels, useNA="ifany") #useNA can be used turned on in the `table` function

# Inspect quality of the predictions
# Save as high-resolution JPEG
jpeg("6_harmonyintegrated_singleR_HPCApredictions_quality.jpg",
     units = "in", width = 15, height = 12, res = 300)
plotScoreHeatmap(annot_ENCODE)
dev.off()

#plotDeltaDistribution(annot_HPCA, ncol = 4, dots.on.top = FALSE)

# Add to seurat object
objs <- AddMetaData(objs, annot_HPCA$pruned.labels, col.name = 'SingleR_HPCA')

# Visualise them on the UMAP
DimPlot(objs, reduction= "umap.harmony", group.by = "SingleR_HPCA", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedUMAP_singleRautomatedcellannotations_HPCA.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedUMAP_singleRautomatedcellannotations_HPCA.jpg", units="in", width=12, height=8, res=300)
dev.off()

# Visualise them on the tSNE
DimPlot(objs, reduction= "tSNE.harmony", group.by = "SingleR_HPCA", label = T , repel = T, label.size = 3) 

#TO save image 
ggsave("6_harmonyintegratedtSNE_singleRautomatedcellannotations_HPCA.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedtSNE_singleRautomatedcellannotations_HPCA.jpg", units="in", width=12, height=8, res=300)
dev.off()