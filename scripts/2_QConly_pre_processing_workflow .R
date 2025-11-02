#Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

##### Standard pre-processing workflow (includes 4 steps)
# 1. Quality control and selecting cells for further analysis
  # Check QC on individual samples and then do QC on joint samples
# 2. Data normalization 
# 3. Identification of highly variable features (feature selection)
# 4. Data scaling

# load the libraries if required
library(Seurat)
library(tidyverse)

# Read Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes-SKIP IF NOT REQUIRED)
#All controls
patient_1_N <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_1_N.RDS")
patient_2_N <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_2_N.RDS")
patient_3_N <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_3_N.RDS")
#All cases
patient_1_T <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_1_T.RDS")
patient_2_T <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_2_T.RDS")
patient_3_T <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_3_T.RDS")

# 1. Quality control and selecting cells for further analysis #
##### QC metrics: "nFeature_RNA", "nCount_RNA", "percent.mt" 
# Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics commonly used by the community include

#The number of unique genes detected in each cell.
#Low-quality cells or empty droplets will often have very few genes (low nFeature_RNA & nCount_RNA)
# Cell doublets or multiplets may exhibit an aberrantly high gene count ( high values of nFeature_RNA & nCount_RNA)

# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)

#The percentage of reads that map to the mitochondrial genome
#Low-quality / dying cells often exhibit extensive mitochondrial contamination (high percent.mt)
#We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
#We use the set of all genes starting with MT- as a set of mitochondrial genes

#From this link: https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/Seurat_QC_to_Clustering/
#nCount_RNA - the absolute number of RNA molecules (UMIs) per cell (i.e., count depth). Each unique RNA molecule (non-PCR duplicates) will have its own Unique Molecular Identifier (UMI).
#high total count - potential doublets or multiplets
#low total count - potential ambient mRNA (not real cells)
#Cell Ranger threshold set at 500 UMIs

#nFeature_RNA - number of genes expressed (detected) per cell.
#high number of detected genes - potential doublets or multiplets
#low number of detected genes - potential ambient mRNA (not real cells)

#Percent mitochondrial(percent.mt) - the fraction of reads from mitochondrial genes per cell
#high mtDNA - cellular degradation

#In general, low quality cells and empty droplets will have few genes (low n_feature_RNA, low nCount_RNA), whereas doublets and multiplets will exhibit high n_feature_RNA and high nCount_RNA. Low quality or dying cells will also have high percent_mt.

# View individual sample QC metrics
#For controls: patient_1_N
range(patient_1_N$nFeature_RNA)
range(patient_1_N$nCount_RNA)

# store mitochondrial percentage in object meta data
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
patient_1_N[["percent.mt"]] <- PercentageFeatureSet(patient_1_N, pattern = "^MT-")
view(patient_1_N@meta.data)
range(patient_1_N$percent.mt) 

#include Ribosomal content
patient_1_N[["percent.ribo"]] <- PercentageFeatureSet(patient_1_N, pattern = "^RP[SL]")

#include hemoglobin genes (blodd contamination risk during sample isolation)
patient_1_N[["percent.hb"]] <- PercentageFeatureSet(patient_1_N, pattern = "^HB[^(P|E|S)]")

#we have to get rid of genes with very high mt content based on this range

##### Selecting cells for further analysis
# Visualize QC metrics as a violin plot for all the cells
VlnPlot(patient_1_N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)
patient_1_N

#TO save image 
ggsave("2_raw_patient_1_N_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_1_N_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used if required
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(patient_1_N, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plot2 <- FeatureScatter(patient_1_N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plot3 <- FeatureScatter(patient_1_N, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plot4 <- FeatureScatter(patient_1_N, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plot1 + plot2 + plot3 + plot4

#TO save image 
ggsave("2_raw_patient_1_N_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_1_N_featurescatterplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#straight line is added to the feature scatter plot using geom line
#good quality data set follows the straight line typically.
# Especially in plot2, we should not expect to see any cells in the top left corner (experiment captured only higher number of genes which are deeply sequenced enough) or lower right corner (experiment captured only lower number of genes repeatedly which results in higher transcript or RNA counts). These are lower quality cells. 
------------------------------------------------------------------------------------------------------------------------------------
#For controls: patient_2_N
range(patient_2_N$nFeature_RNA)
range(patient_2_N$nCount_RNA)
patient_2_N[["percent.mt"]] <- PercentageFeatureSet(patient_2_N, pattern = "^MT-")
view(patient_2_N@meta.data)
range(patient_2_N$percent.mt) 
patient_2_N[["percent.ribo"]] <- PercentageFeatureSet(patient_2_N, pattern = "^RP[SL]")
patient_2_N[["percent.hb"]] <- PercentageFeatureSet(patient_2_N, pattern = "^HB[^(P|E|S)]")

#we have to get rid of genes with very high mt content based on this range

VlnPlot(patient_2_N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)
patient_2_N
#TO save image 
ggsave("2_raw_patient_2_N_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_2_N_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

plot5 <- FeatureScatter(patient_2_N, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plot6 <- FeatureScatter(patient_2_N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plot7 <- FeatureScatter(patient_2_N, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plot8 <- FeatureScatter(patient_2_N, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plot5 + plot6+ plot7 + plot8
#TO save image 
ggsave("2_raw_patient_2_N_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_2_N_featurescatterplot.jpg", units="in", width=10, height=8, res=300)
dev.off()
---------------------------------------------------------------------------------------------------------------------------------------
#For controls: patient_3_N
range(patient_3_N$nFeature_RNA)
range(patient_3_N$nCount_RNA)
patient_3_N[["percent.mt"]] <- PercentageFeatureSet(patient_3_N, pattern = "^MT-")
view(patient_3_N@meta.data)
range(patient_3_N$percent.mt) 
patient_3_N[["percent.ribo"]] <- PercentageFeatureSet(patient_3_N, pattern = "^RP[SL]")
patient_3_N[["percent.hb"]] <- PercentageFeatureSet(patient_3_N, pattern = "^HB[^(P|E|S)]")

#we have to get rid of genes with very high mt content based on this range

VlnPlot(patient_3_N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)
patient_3_N
#TO save image 
ggsave("2_raw_patient_3_N_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_3_N_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

plot9 <- FeatureScatter(patient_3_N, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plot10 <- FeatureScatter(patient_3_N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plot11 <- FeatureScatter(patient_3_N, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plot12 <- FeatureScatter(patient_3_N, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plot9 + plot10+ plot11 + plot12
#TO save image 
ggsave("2_raw_patient_3_N_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_3_N_featurescatterplot.jpg", units="in", width=10, height=8, res=300)
dev.off()
---------------------------------------------------------------------------------------------------------------------------------------
#For cases :patient_1_T
range(patient_1_T$nFeature_RNA)
range(patient_1_T$nCount_RNA)
patient_1_T[["percent.mt"]] <- PercentageFeatureSet(patient_1_T, pattern = "^MT-")
view(patient_1_T@meta.data)
range(patient_1_T$percent.mt) 

patient_1_T[["percent.ribo"]] <- PercentageFeatureSet(patient_1_T, pattern = "^RP[SL]")
patient_1_T[["percent.hb"]] <- PercentageFeatureSet(patient_1_T, pattern = "^HB[^(P|E|S)]")

#we have to get rid of genes with very high mt content based on this range

VlnPlot(patient_1_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)
patient_1_T
#TO save image 
ggsave("2_raw_patient_1_T_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_1_T_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

plot13 <- FeatureScatter(patient_1_T, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plot14 <- FeatureScatter(patient_1_T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plot15 <- FeatureScatter(patient_1_T, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plot16 <- FeatureScatter(patient_1_T, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plot13 + plot14+ plot15 + plot16
#TO save image 
ggsave("2_raw_patient_1_T_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_1_T_featurescatterplot.jpg", units="in", width=10, height=8, res=300)
dev.off()
---------------------------------------------------------------------------------------------------------------------------------------
#For cases :patient_2_T
range(patient_2_T$nFeature_RNA)
range(patient_2_T$nCount_RNA)
patient_2_T[["percent.mt"]] <- PercentageFeatureSet(patient_2_T, pattern = "^MT-")
view(patient_2_T@meta.data)
range(patient_2_T$percent.mt) 

patient_2_T[["percent.ribo"]] <- PercentageFeatureSet(patient_2_T, pattern = "^RP[SL]")
patient_2_T[["percent.hb"]] <- PercentageFeatureSet(patient_2_T, pattern = "^HB[^(P|E|S)]")

#we have to get rid of genes with very high mt content based on this range

VlnPlot(patient_2_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)
patient_2_T
#TO save image 
ggsave("2_raw_patient_2_T_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_2_T_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

plot17 <- FeatureScatter(patient_2_T, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plot18 <- FeatureScatter(patient_2_T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plot19 <- FeatureScatter(patient_2_T, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plot20 <- FeatureScatter(patient_2_T, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plot17 + plot18+ plot19 + plot20
#TO save image 
ggsave("2_raw_patient_2_T_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_2_T_featurescatterplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

---------------------------------------------------------------------------------------------------------------------------------------
#For cases :patient_3_T
range(patient_3_T$nFeature_RNA)
range(patient_3_T$nCount_RNA)
patient_3_T[["percent.mt"]] <- PercentageFeatureSet(patient_3_T, pattern = "^MT-")
view(patient_3_T@meta.data)
range(patient_3_T$percent.mt) 
patient_3_T[["percent.ribo"]] <- PercentageFeatureSet(patient_3_T, pattern = "^RP[SL]")
patient_3_T[["percent.hb"]] <- PercentageFeatureSet(patient_3_T, pattern = "^HB[^(P|E|S)]")

#we have to get rid of genes with very high mt content based on this range

VlnPlot(patient_3_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)
patient_3_T
#TO save image 
ggsave("2_raw_patient_3_T_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_3_T_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

plot21 <- FeatureScatter(patient_3_T, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plot22 <- FeatureScatter(patient_3_T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plot23 <- FeatureScatter(patient_3_T, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plot24 <- FeatureScatter(patient_3_T, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plot21 + plot22+ plot23 + plot24
#TO save image 
ggsave("2_raw_patient_3_T_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_raw_patient_3_T_featurescatterplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

---------------------------------------------------------------------------------------------------------------------------------------

#Filter further to get the best quality of cells for downstream analysis (filter out low quality cells)
#Controls:patient_1_N
patient_1_N_filtered_1 <- subset(patient_1_N, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA <= 20000 & percent.mt <= 20 & percent.hb<=1)
VlnPlot(patient_1_N_filtered_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)

#TO save image 
ggsave("2_filtered_patient_1_N_filtered_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_1_N_filtered_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Feature scatter plots to recheck quality IF REQUIRED
plota <- FeatureScatter(patient_1_N_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plotb <- FeatureScatter(patient_1_N_filtered_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plotc <- FeatureScatter(patient_1_N_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plotd <- FeatureScatter(patient_1_N_filtered_1, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plota + plotb + plotc + plotd
#TO save image 
ggsave("2_filtered_patient_1_N_filtered_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_1_N_filtered_featurscattereplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#No need to adjust/QC since ribosomal content look s good after mito adjustment, but can still regress during normalization if required as a precaution.
---------------------------------------------------------------------------------------------------------------------------------------
#Controls:patient_2_N
patient_2_N_filtered_1 <- subset(patient_2_N, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA <= 20000 & percent.mt <= 20 & percent.hb<=1)
VlnPlot(patient_2_N_filtered_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)

#TO save image 
ggsave("2_filtered_patient_2_N_filtered_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_2_N_filtered_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Feature scatter plots to recheck quality IF REQUIRED
plote <- FeatureScatter(patient_2_N_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plotf <- FeatureScatter(patient_2_N_filtered_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plotg <- FeatureScatter(patient_2_N_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
ploth <- FeatureScatter(patient_2_N_filtered_1, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plote + plotf + plotg + ploth
#TO save image 
ggsave("2_filtered_patient_2_N_filtered_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_2_N_filtered_featurscattereplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#No need to adjust/QC since ribosomal content look s good after mito adjustment, but can still regress during normalization if required as a precaution.
---------------------------------------------------------------------------------------------------------------------------------------
#Controls:patient_3_N
patient_3_N_filtered_1 <- subset(patient_3_N, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA <= 20000 & percent.mt <= 20 & percent.hb<=1)
VlnPlot(patient_3_N_filtered_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)

#TO save image 
ggsave("2_filtered_patient_3_N_filtered_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_3_N_filtered_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Feature scatter plots to recheck quality IF REQUIRED
ploti <- FeatureScatter(patient_3_N_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plotj<- FeatureScatter(patient_3_N_filtered_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plotk <- FeatureScatter(patient_3_N_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plotl <- FeatureScatter(patient_3_N_filtered_1, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
ploti + plotj + plotk + plotl
#TO save image 
ggsave("2_filtered_patient_3_N_filtered_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_3_N_filtered_featurscattereplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

---------------------------------------------------------------------------------------------------------------------------------------
#Cases:patient_1_T
patient_1_T_filtered_1 <- subset(patient_1_T, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA <= 20000 & percent.mt <= 20 & percent.hb<=1)
VlnPlot(patient_1_T_filtered_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)

#TO save image 
ggsave("2_filtered_patient_1_T_filtered_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_1_T_filtered_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Feature scatter plots to recheck quality IF REQUIRED
plotm <- FeatureScatter(patient_1_T_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plotn<- FeatureScatter(patient_1_T_filtered_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
ploto <- FeatureScatter(patient_1_T_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plotp <- FeatureScatter(patient_1_T_filtered_1, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plotm + plotn + ploto + plotp
#TO save image 
ggsave("2_filtered_patient_1_T_filtered_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_1_T_filtered_featurscattereplot.jpg", units="in", width=10, height=8, res=300)
dev.off()


---------------------------------------------------------------------------------------------------------------------------------------
  
#Cases:patient_2_T
patient_2_T_filtered_1 <- subset(patient_2_T, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA <= 20000 & percent.mt <= 20 & percent.hb<=1)
VlnPlot(patient_2_T_filtered_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)

#TO save image 
ggsave("2_filtered_patient_2_T_filtered_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_2_T_filtered_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Feature scatter plots to recheck quality IF REQUIRED
plotq <- FeatureScatter(patient_2_T_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plotr<- FeatureScatter(patient_2_T_filtered_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
plots <- FeatureScatter(patient_2_T_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plott <- FeatureScatter(patient_2_T_filtered_1, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plotq + plotr + plots + plott
#TO save image 
ggsave("2_filtered_patient_2_T_filtered_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_2_T_filtered_featurscattereplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

---------------------------------------------------------------------------------------------------------------------------------------  
#Cases:patient_3_T
patient_3_T_filtered_1 <- subset(patient_3_T, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA <= 20000 & percent.mt <= 20 & percent.hb<=1)
VlnPlot(patient_3_T_filtered_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)

#TO save image 
ggsave("2_filtered_patient_3_T_filtered_vinplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_3_T_filtered_vinplot.jpg", units="in", width=10, height=8, res=300)
dev.off()

#Feature scatter plots to recheck quality IF REQUIRED
plotw <- FeatureScatter(patient_3_T_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method='lm')
plotx<- FeatureScatter(patient_3_T_filtered_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')
ploty <- FeatureScatter(patient_3_T_filtered_1, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_smooth(method='lm')
plotz <- FeatureScatter(patient_3_T_filtered_1, feature1 = "percent.ribo", feature2 = "percent.mt") + geom_smooth(method='lm')
plotw + plotx + ploty + plotz
#TO save image 
ggsave("2_filtered_patient_3_T_filtered_featurscattereplot.jpg", units="in", width=10, height=8, dpi=300)
jpeg("2_filtered_patient_3_T_filtered_featurscattereplot.jpg", units="in", width=10, height=8, res=300)
dev.off()
  

# Save Seurat objects 
saveRDS(patient_1_N_filtered_1, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_1_N_filtered.RDS")
saveRDS(patient_2_N_filtered_1, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_2_N_filtered.RDS")
saveRDS(patient_3_N_filtered_1, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_3_N_filtered.RDS")
saveRDS(patient_1_T_filtered_1, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_1_T_filtered.RDS")
saveRDS(patient_2_T_filtered_1, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_2_T_filtered.RDS")
saveRDS(patient_3_T_filtered_1, file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2_patient_3_T_filtered.RDS")