#Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

#Run for the first time
#BiocManager::install("glmGamPoi") #needed for sctransform faster run

#To fix errors
#install.packages("devtools")
#devtools::install_github("satijalab/seurat", ref = "fix/v.5.3.1")

#For the first time only
install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)

#install UMAP if not already installed
#reticulate::py_install(packages='umap-learn')  

# load the libraries if required
library(Seurat)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(DoubletFinder)

# Read Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
patient_3_T<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_3_T_scTransformed.RDS")
#Check default assay
DefaultAssay(object = patient_3_T) 

#SCT is default which we need for LDR/clustering before doubletfinder.
#if not change using DefaultAssay(object = adp_filt) <- "SCT"
#Run Linear dimensionality reduction or PCA for each sample and determine dimensionality using elbow plot
#Run Clustering and non-linear dimensionality reduction for each samples using findneighbours, findclusters and UMAP/tSNE

patient_3_T  <- RunPCA(patient_3_T, features = VariableFeatures(object = patient_3_T))
# Examine and visualize PCA results a few different ways
print(patient_3_T [["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(patient_3_T, dims = 1:2, reduction = "pca") 
#TO save image 
ggsave("4_patient_3_T_PCA_VizDimLoadings.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_PCA_VizDimLoadings.jpg", units="in", width=10, height=8, res=300)
dev.off()

DimPlot(patient_3_T, reduction = "pca") # + NoLegend() will remove patient annotations on the side
#TO save image 
ggsave("4_patient_3_T_PC1_PC2_pcaplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_PC1_PC2_pcaplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#TO save image 
jpeg("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_3_T_PC1_to_15_dimheatmap.jpg",
     width=10,
     height=8,
     units="in",
     res=300)
#par(mar = c(10, 10, 3, 3))

DimHeatmap(patient_3_T, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#Determine the ‘dimensionality’ of the dataset
#top Principal components represent robust data. So to decide which statistical significant PCs to include for downstream analysis.
#After elbow, PCs do not vary much. 
#Also, good to be liberal with PC inclusion for downstream analysis and not to be conservative with number of PCs chosen for the analysis. So around first 15- 20 PCs is good.
ElbowPlot(patient_3_T) #more faster and straightforward than the jackstraw plot that Seurat previously used.
#this by default runs only the first 20 PCs
#TO save image 
ggsave("4_patient_3_T_elbowplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_elbowplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

ElbowPlot(patient_3_T, ndims = 50, reduction = "pca")
#we can observe an ‘elbow’ around PC9-12, suggesting that the majority of true signal is captured in the first 12 PCs.
#But we will choose 15-20 PCs for downstream analysis.
#TO save image 
ggsave("4_patient_3_T_elbowplot_50dimensions.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_elbowplot_50dimensions.jpg", units="in", width=10, height=8, res=300)
dev.off()

#clustering and non-linear dimensionality reduction
patient_3_T <- FindNeighbors(patient_3_T, dims = 1:20) # USE SAME PCs in UMAP/ tSNE plot
patient_3_T <- FindClusters(patient_3_T)
patient_3_T  <- RunUMAP(patient_3_T , dims = 1:20)
DimPlot(patient_3_T , reduction = "umap", label = TRUE, repel = TRUE)
#TO save image 
ggsave("4_patient_3_T_UMAP.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_UMAP.jpg", units="in", width=10, height=8, res=300)
dev.off()

patient_3_T <- RunTSNE(object = patient_3_T)
DimPlot(object = patient_3_T, reduction = "tsne")
#TO save image 
ggsave("4_patient_3_T_tSNE.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_tSNE.jpg", units="in", width=10, height=8, res=300)
dev.off()

#DoubletFinder
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_patient_3_T <- paramSweep(patient_3_T, PCs = 1:20, sct = TRUE)
sweep.stats_patient_3_T <- summarizeSweep(sweep.res.list_patient_3_T, GT = FALSE)
bcmvn_patient_3_T <- find.pK(sweep.stats_patient_3_T)

ggplot(bcmvn_patient_3_T, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()
#TO save image 
ggsave("4_patient_3_T_DoubletFinder_bcvmnplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_DoubletFinder_bcvmnplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

pK <- bcmvn_patient_3_T %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- patient_3_T@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(patient_3_T@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
patient_3_T <- doubletFinder(patient_3_T, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = NULL, sct = TRUE)
patient_3_T

# visualize doublets
DimPlot(patient_3_T, reduction = 'umap', group.by = "DF.classifications_0.25_0.09_616")

# number of singlets and doublets
count_singlet_doublet<- data.frame(table(patient_3_T@meta.data$DF.classifications_0.25_0.09_616))
count_singlet_doublet

#export as CSV
write.csv(count_singlet_doublet, "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_3_T_DoubletFinder_singletORdoubletnumbers.csv", row.names = TRUE)

# grab the correct DF classification column dynamically
df_col <- grep("^DF.classifications", colnames(patient_3_T@meta.data), value = TRUE)[1]

# visualize doublets vs singlets on UMAP
DimPlot(
  patient_3_T,
  reduction = "umap",
  group.by  = df_col,
  pt.size   = 0.4
) + ggtitle("patient_3_T: DoubletFinder classifications")

#TO save image 
ggsave("4_patient_3_T_DoubletFinder_singletORdoubletclassifications_UMAP.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_3_T_DoubletFinder_singletORdoubletclassifications_UMAP.jpg", units="in", width=10, height=8, res=300)
dev.off()

#pick singlet cell barcodes
cells_keep <- rownames(patient_3_T@meta.data)[ patient_3_T@meta.data[[df_col]] == "Singlet" ]
# subset by explicit cell list
patient_3_T_singlets <- subset(patient_3_T, cells = cells_keep)

#quick sanity check
table(patient_3_T_singlets@meta.data[[df_col]])

# Save Seurat objects 
saveRDS(patient_3_T_singlets, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_3_T_singletsonly.RDS")
