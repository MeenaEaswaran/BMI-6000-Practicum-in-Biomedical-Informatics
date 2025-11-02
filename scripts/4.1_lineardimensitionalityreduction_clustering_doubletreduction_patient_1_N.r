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
patient_1_N<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3_patient_1_N_scTransformed.RDS")
#Check default assay
DefaultAssay(object = patient_1_N) 

#SCT is default which we need for LDR/clustering before doubletfinder.
#if not change using DefaultAssay(object = adp_filt) <- "SCT"
#Run Linear dimensionality reduction or PCA for each sample and determine dimensionality using elbow plot
#Run Clustering and non-linear dimensionality reduction for each samples using findneighbours, findclusters and UMAP/tSNE

patient_1_N  <- RunPCA(patient_1_N, features = VariableFeatures(object = patient_1_N))
# Examine and visualize PCA results a few different ways
print(patient_1_N [["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(patient_1_N, dims = 1:2, reduction = "pca") 
#TO save image 
ggsave("4_patient_1_N_PCA_VizDimLoadings.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_PCA_VizDimLoadings.jpg", units="in", width=10, height=8, res=300)
dev.off()

DimPlot(patient_1_N, reduction = "pca") # + NoLegend() will remove patient annotations on the side
#TO save image 
ggsave("4_patient_1_N_PC1_PC2_pcaplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_PC1_PC2_pcaplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#TO save image 
jpeg("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_1_N_PC1_to_15_dimheatmap.jpg",
     width=10,
     height=8,
     units="in",
     res=300)
#par(mar = c(10, 10, 3, 3))

DimHeatmap(patient_1_N, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#Determine the ‘dimensionality’ of the dataset
#top Principal components represent robust data. So to decide which statistical significant PCs to include for downstream analysis.
#After elbow, PCs do not vary much. 
#Also, good to be liberal with PC inclusion for downstream analysis and not to be conservative with number of PCs chosen for the analysis. So around first 15- 20 PCs is good.
ElbowPlot(patient_1_N) #more faster and straightforward than the jackstraw plot that Seurat previously used.
#this by default runs only the first 20 PCs
#TO save image 
ggsave("4_patient_1_N_elbowplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_elbowplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

ElbowPlot(patient_1_N, ndims = 50, reduction = "pca")
#we can observe an ‘elbow’ around PC9-12, suggesting that the majority of true signal is captured in the first 12 PCs.
#But we will choose 15-20 PCs for downstream analysis.
#TO save image 
ggsave("4_patient_1_N_elbowplot_50dimensions.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_elbowplot_50dimensions.jpg", units="in", width=10, height=8, res=300)
dev.off()


#clustering and non-linear dimensionality reduction
patient_1_N <- FindNeighbors(patient_1_N, dims = 1:20) # USE SAME PCs in UMAP/ tSNE plot
patient_1_N <- FindClusters(patient_1_N)
patient_1_N  <- RunUMAP(patient_1_N , dims = 1:20)
DimPlot(patient_1_N , reduction = "umap", label = TRUE, repel = TRUE)
#TO save image 
ggsave("4_patient_1_N_UMAP.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_UMAP.jpg", units="in", width=10, height=8, res=300)
dev.off()

patient_1_N <- RunTSNE(object = patient_1_N)
DimPlot(object = patient_1_N, reduction = "tsne")
#TO save image 
ggsave("4_patient_1_N_tSNE.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_tSNE.jpg", units="in", width=10, height=8, res=300)
dev.off()

#DoubletFinder
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_patient_1_N <- paramSweep(patient_1_N, PCs = 1:20, sct = TRUE)
sweep.stats_patient_1_N <- summarizeSweep(sweep.res.list_patient_1_N, GT = FALSE)
bcmvn_patient_1_N <- find.pK(sweep.stats_patient_1_N)

ggplot(bcmvn_patient_1_N, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()
#TO save image 
ggsave("4_patient_1_N_DoubletFinder_bcvmnplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_DoubletFinder_bcvmnplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

pK <- bcmvn_patient_1_N %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- patient_1_N@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(patient_1_N@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
patient_1_N <- doubletFinder(patient_1_N, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = NULL, sct = TRUE)
patient_1_N

# visualize doublets
DimPlot(patient_1_N, reduction = 'umap', group.by = "DF.classifications_0.25_0.17_98")

# number of singlets and doublets
count_singlet_doublet<- data.frame(table(patient_1_N@meta.data$DF.classifications_0.25_0.17_98))

#export as CSV
write.csv(count_singlet_doublet, "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_1_N_DoubletFinder_singletORdoubletnumbers.csv", row.names = TRUE)

# grab the correct DF classification column dynamically
df_col <- grep("^DF.classifications", colnames(patient_1_N@meta.data), value = TRUE)[1]

# visualize doublets vs singlets on UMAP
DimPlot(
  patient_1_N,
  reduction = "umap",
  group.by  = df_col,
  pt.size   = 0.4
) + ggtitle("patient_1_N: DoubletFinder classifications")

#TO save image 
ggsave("4_patient_1_N_DoubletFinder_singletORdoubletclassifications_UMAP.jpg", units="in", width=10, height=8, dpi=300)
jpeg("4_patient_1_N_DoubletFinder_singletORdoubletclassifications_UMAP.jpg", units="in", width=10, height=8, res=300)
dev.off()

#pick singlet cell barcodes
cells_keep <- rownames(patient_1_N@meta.data)[ patient_1_N@meta.data[[df_col]] == "Singlet" ]
# subset by explicit cell list
patient_1_N_singlets <- subset(patient_1_N, cells = cells_keep)

#quick sanity check
table(patient_1_N_singlets@meta.data[[df_col]])

# Save Seurat objects 
saveRDS(patient_1_N_singlets, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_1_N_singletsonly.RDS")
