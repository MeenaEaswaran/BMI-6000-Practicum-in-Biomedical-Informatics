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
options(future.globals.maxSize = 1e9)
#library(DoubletFinder)

#Read Harmony integrated Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
objs<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_harmony_integrated_seurat_30dims.rds")

#Check default assay
DefaultAssay(object = objs) 
DefaultAssay(objs) <- "SCT" #set if not

# Rename all identities based on manual observations
objs <- RenameIdents(object = objs, 
                                  "0" = "Epithelial cells",
                                  "1" = "T-cells",
                                  "2" = "Epithelial cells",
                                  "3" = "Epithelial cells",
                                  "4" = "Epithelial cells",
                                  "5" = "Myeloid cells and NK cells",
                                  "6" = "Epithelial cells",
                                  "7" = "Epithelial cells",
                                  "8" = "Epithelial cells",
                                  "9" = "Epithelial cells",
                                  "10" = "Epithelial cells",
                                  "11" = "T-cells",
                                  "12" = "Epithelial cells",
                                  "13" = "Epithelial cells",
                                  "14" = "Fibroblasts",
                                  "15" = "T-cells",
                                  "16" = "Myeloid cells and NK cells",
                                  "17" = "B-cells",  
                                  "18" = "B-cells", 
                                  "19" = "Endothelial cells", 
                                  "20" = "Endothelial cells",
                                  "21" = "Epithelial cells",
                                  "22" = "T-cells",
                                  "23" = "T-cells")

# Plot the UMAP
plot1 <- DimPlot(object = objs, 
        reduction = "umap.harmony", 
        label = TRUE,
        label.size = 5,
        #label.box = TRUE,
        repel = TRUE) + ggtitle("UMAP—Manually Annotated Cell Clusters on Harmony Integrated data")

plot1

#TO save image 
ggsave("6_harmonyintegratedUMAP_manualcellannotations.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedUMAP_manualcellannotations.jpg", units="in", width=12, height=8, res=300)
dev.off()

# Plot the tSNE
# First run the tSNE analysis
objs <- RunTSNE(
  objs,
  reduction      = "harmony",
  dims           = 1:30,
  reduction.name = "tSNE.harmony"   
)

plot2 <- DimPlot(object = objs, 
                 reduction = "tSNE.harmony", 
                 label = FALSE,
                 #label.size = 5,
                 #label.box = TRUE,
                 #repel = TRUE
                 ) + ggtitle("t-SNE—Manually Annotated Cell Clusters on Harmony Integrated data")

plot2

#TO save image 
ggsave("6_harmonyintegratedtSNE_manualcellannotations.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedtSNE_manualcellannotations.jpg", units="in", width=12, height=8, res=300)
dev.off()

###########
# make sure Idents are your manual labels
# Idents(objs) are already set by your RenameIdents() step

library(patchwork)

cell_types <- c(
  "Epithelial cells","T-cells","B-cells",
  "Myeloid cells and NK cells","Fibroblasts","Endothelial cells"
)

plots_by_type <- lapply(cell_types, function(ct) {
  DimPlot(
    objs,
    reduction = "tSNE.harmony",
    cells.highlight = WhichCells(objs, idents = ct),
    cols = c("grey85", "red"),
    raster = TRUE, pt.size = 0.2
  ) + ggtitle(ct) + NoLegend()
})

# arrange like your figure (3 columns × 2 rows)
p <- wrap_plots(plots_by_type, ncol = 3)
p

#TO save image 
ggsave("6_harmonyintegratedtSNE_manualcellannotations_onecelltypehighlight.jpg", units="in", width=12, height=8, dpi=300)
jpeg("6_harmonyintegratedtSNE_manualcellannotations_onecelltypehighlight.jpg", units="in", width=12, height=8, res=300)
dev.off()

# Save harmony integrated Seurat object with tSNE reduction now
saveRDS(objs, "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithtSNE_seurat_30dims.rds")
