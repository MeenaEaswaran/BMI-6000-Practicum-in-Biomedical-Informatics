# Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

# SCENIC (Single cell regulatory network interference and clustering)

#Introduction to SCENIC
# SCENIC is a tool to simultaneously reconstruct gene regulatory networks and identify stable cell states from single-cell RNA-seq data. The gene regulatory network is inferred based on co-expression and DNA motif analysis, and then the network activity is analyzed in each cell to identify the recurrent cellular states.

#Required, SCENIC R is based on three R packages (for first time only)
BiocManager::install(c("AUCell", "RcisTarget"), force=TRUE)
BiocManager::install(c("GENIE3"), force = TRUE) # Can be replaced by GRNBoost (Linux & MacOS)


#Optional (but highly recommended):
#To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools"), force = TRUE)
#For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"), force = TRUE)
#To support parallel execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"), force = TRUE)

remotes::install_github("bokeh/rbokeh", force = TRUE)

#SCENIC
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC", force = TRUE) 
packageVersion("SCENIC")

devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE, force = TRUE)


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

#install cellchatDB and dependencies for the first time
#devtools::install_github("jinworks/CellChat")
#install.packages('NMF')
#devtools::install_github("jokergoo/circlize")
#devtools::install_github("jokergoo/ComplexHeatmap")
#UMAP-learn needed but already installed
#install.packages("future.callr")
# load the libraries if required
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
#library(sctransform)
library(glmGamPoi)
library(harmony)
library(metap)
library(multtest)
library(presto)
library(SingleR)
library(celldex)
#library(scRNAseq)
#library(SingleCellExperiment)
library(pheatmap)
#options(future.globals.maxSize = 1e9)
#library(CellChat)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(NMF)
library(compiler)
library(parallel)
#library(future.callr)
#library(DoubletFinder)
options(stringsAsFactors = FALSE)
library(reticulate)
library(SCopeLoomR)

#Create a directory to keep all input output files in one place
#This is where the cisTarget_databases folder should also be stored
dir.create("SCENIC_LSCC_2percent")
setwd("SCENIC_LSCC_2percent")

### Load data
loom_path <- "~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/10_scenic_input.loom"

loom <- open_loom(loom_path)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)

#Cell info/phenodata
# cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)

cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$SingleR_ENCODE))

dir.create("int")

saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(SingleR_ENCODE=c("Epithelial cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Myeloid cells"="magenta4", 
                           "B-cells"="hotpink", 
                           "NK cells"="red3", 
                           "T-cells"="skyblue", 
                           "Endothelial cells"="darkblue"))
colVars$SingleR_ENCODE <- colVars$SingleR_ENCODE[intersect(names(colVars$SingleR_ENCODE), cellInfo$SingleR_ENCODE)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$SingleR_ENCODE, legend=names(colVars$SingleR_ENCODE))


#Initialize settings
#In order to keep consistent settings across the multiple steps of SCENIC, most functions in SCENIC package use a common object where the options for the current run are stored. This object replaces the “arguments” for most functions, and should be created at the begining of a SCENIC run with the function initializeScenic().

#The default settings should be valid for most analyses. The parameters that need to be specified in all runs is the organism (mgi for mouse, hgnc for human, or dmel for fly), and the directory where the RcisTarget databases are stored (you may create a link in the current directory to avoid duplicating them, e.g. in linux: system("ln -s ~/path/to/dbs databases")).

#For details on the options that can be modified check the help of ?initializeScenic or of the specific function that takes it as input.

#Species-specific databases: human, mouse and Drosophila (use human)
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs
#These should be downloaded from the old folders in cisTarget databases as the SCENIC R version precedes the latest pySCENIC version


library(SCENIC)
org <- "hgnc" # or mgi, or dmel
dbDir <- "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC LSCC_2percent" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

#This error message will pop-up but if you see the environment, there will be motifAnnotations datatable that would have autopopulated that you can set as motifaddnotations_hgnc.

motifAnnotations_hgnc <- motifAnnotations

#RERUN SCENIC OPTIONS NOW. The error message would have gone and the below will pop-up.

#Motif databases selected: 
#hg19-500bp-upstream-7species.mc9nr.feather 
#hg19-tss-centered-10kb-7species.mc9nr.feather
#Using the column 'features' as feature index for the ranking database.
#Using the column 'features' as feature index for the ranking database.


#Add cellInfo and colVars
scenicOptions@inputDatasetInfo$cellInfo <- 
  "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/cellInfo.Rds"

scenicOptions@inputDatasetInfo$colVars <- 
  "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/colVars.Rds"

#Save to use at a later time
saveRDS(scenicOptions, 
        file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/scenicOptions.Rds") 

########################################################################################################################

scenicOptions <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/scenicOptions.Rds")

#Co-expression network
#1: Gene filter/selection
#(Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#adjusted with stricter calculations (moderate strict (2%)) - used this 
genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions,
                           minCountsPerGene = 3 * .02 * ncol(exprMat),
                           minSamples = ncol(exprMat) * .02)

#We have 8494 genes after filtering with adjusted features.

#check interesting genes
interestingGenes <- c("SOX9", "SOX10", "DLX5", "SOX2", "SOX4") #human nomenclature
#what is missing?
interestingGenes[which(!interestingGenes %in% genesKept)] 
#  "SOX10" "DLX5"  are missing. 

#filter the expression matrix to contain only these 10333 genes. 
#This matrix is now ready for the co-expression analysis.
#How many genes and cells are left after filtering your expression matrix.
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
#8494 genes across 30117 cells.

#2: Correlation (+ or - correlations translates to activation or suppression by TF's on downstream genes)
#In order to distinguish potential activation from repression, we will split the targets into positive- and negative-correlated targets (i.e. Spearman correlation between the TF and the potential target).
runCorrelation(exprMat_filtered, scenicOptions)

#If correlation gives vector memory limit issue, do the below, IF REQUIRED.

# SINCE THE VECTOR MEMORY EXHAUST LIMIT REACHES WHEN RUNNING PCA FOR TUMOR SAMPLES, RSTUDIO LIMIT HAS TO BE RESET USING TERMINAL IN RSTUDIO AND THEN APPLICATION HAS TO BE RESTARTED. FOLLOW THE CODE/INSTRUCTIONS ON STACKOVERFLOW AS SHOWN BELOW

#https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos

library(usethis) 
usethis::edit_r_environ()
#when the tab opens up in R studio, add this to the 1st line: R_MAX_VSIZE=100Gb (or whatever memory you wish to allocate).
#Re-start R and/or restart computer and run the R command again that gave you the memory error.
#RERUN ABOVE COMMANDS after loading relevant libraries.

#3: Run GENIE3 to infer potential transcription factor targets
#Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)
#In runGenie3(exprMat_filtered, scenicOptions) :
#Only 611 (38%) of the 1605 TFs in the database were found in the dataset with moderate strict (2%) filtering options.



#IV: Build the gene regulatory network & Identify cell states:
#Build the gene regulatory network: 
#1: Get co-expression modules
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

#Number of links between TFs and targets (weight>=0.001): 2860201
#[,1]
#nTFs          611
#nTargets     8494
#nGeneSets    3663
#nLinks    4009516

#2: Get regulons (with RcisTarget: TF motif analysis)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) #** Only for toy run!!
#IGNORE warning messages 

#Identify cell states: 
#3: Score GRN (regulons) in the cells (with AUCell) 
#Optional: log expression (for TF expression plot, it does not affect any other calculation)
exprMat_log <- log2(exprMat+1) #normalized the original expression matrix 

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

#Save the analysis
saveRDS(scenicOptions, 
        file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/scenicOptions.Rds")

##########################################################################################################################


BiocManager::install("AUCell", force = TRUE)
install.packages(c("shiny","DT","plotly","NMF","rhandsontable")) 
library(shiny)
library(AUCell)
library(SCENIC)
library(tidyverse)
getAnywhere("AUCell_createViewerApp")      # should now find it

#4: Optional steps

#Visually inspect the 3.5_AUCellThresholds_Info.tsy in the int folder.
#Highest threshold is 0.916 and the lowest is 0.00931.
#4.1: Look the threshold results in shinyAPP to adjust the threshold
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)

#This will not work anymore. SKIP. CANNOT ADJUST THE THRESHOLD. (REPORT THAT DEFAULT THRESHOLD LEVELS WERE USED)
#rbokeh was removed from BioConductor and AUCell was not allowed to depend on it anymore, so the plotTsne_AUCellApp viewer was removed too:aertslab/AUCell#42
#https://github.com/aertslab/SCENIC/issues/448
#https://github.com/aertslab/SCENIC/issues/455
#https://github.com/aertslab/AUCell/issues/42


#4.2: Binarize the network activity (regulon on/off)
#Building the GRN and scoring its activity in AUCell is often enough for datasets with very clear cell types. However, in many cases it is also useful to binarize the activity score into “on/off”; either for easier interpretation, or for maximizing the differences across cell types. It is possible to binarize only specific regulons for exploring/interpreting key TFs. However, binarizing the activity of all the regulons in the dataset allows to create the “Binarized regulon activity matrix”, which can be used for upstream analysis (e.g. clustering). The binarized activity is specially useful to reduce technical biases (e.g. number of detected genes, batch effects), the grouping by patient of origin in cancer datasets, or even cross-species comparisons (see Aibar et al. (2017)).

#To determine in which cells each regulon is active, we will use an AUC threshold. AUCell automatically calculates possible thresholds for the binarization, but they are often too conservative. We recommend to check these thresholds manually before proceeding to the binarization. This can be a iterative process, where the thresholds can be re-adjusted after an initial exploration. Once the final thresholds are selected, the cell-regulon activity will be summarized into a binary activity matrix in which the columns represent the cells and the rows the regulons. The coordinates of the matrix that correspond to active regulons in a given cell will contain a “1” value, and “0” all the others.
# scenicOptions@settings$devType="png"
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions) #USEFUL TO REDUCE TECHNICAL EFFECTS (GROUP DATASETS BY THE ORIGIN OF PATIENTS OR SAMPLES OR ACROSS SPECIES)

#will get an error after the completion of above step as below. 
#Warning message:
#In AUCell_plotTSNE(tSNE = tSNE, exprMat = exprMat, cellsAUC = regulonAUC,  :
                     #Expression plot was requested, but no expression matrix provided.

#Clustering / dimensionality reduction on the regulon activity
#4.3: Cluster cells according to the GRN activity (Optional)
#4.3.1: set number of PCs
nPcs <- c(5) # For demo dataset, # nPcs <- c(5,15,50)

#4.3.2: Calculates the t-SNE based on the regulon activity
scenicOptions@settings$seed <- 123 # same seed for all of them
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames

par(mfrow=c(2,3))
plotTsne_compareSettings("int/tSNE_AUC_05pcs_05perpl.Rds", scenicOptions, 
                         showLegend=FALSE, varName="seuratCluster", cex=.5)

plotTsne_compareSettings("int/tSNE_AUC_05pcs_15perpl.Rds", scenicOptions, 
                         showLegend=FALSE, varName="CellType", cex=.5)

plotTsne_compareSettings("int/tSNE_AUC_05pcs_50perpl.Rds", scenicOptions, 
                         showLegend=FALSE, varName="CellType", cex=.5)

# 4.3.3 Calculates the t-SNE using only "high-confidence" regulons
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50),
                     onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames

plotTsne_compareSettings("int/tSNE_oHC_AUC_05pcs_05perpl.Rds", scenicOptions, 
                         showLegend=FALSE, varName="CellType", cex=.5)

plotTsne_compareSettings("int/tSNE_oHC_AUC_05pcs_15perpl.Rds", scenicOptions, 
                         showLegend=FALSE, varName="CellType", cex=.5)

plotTsne_compareSettings("int/tSNE_oHC_AUC_05pcs_50perpl.Rds", scenicOptions, 
                         showLegend=FALSE, varName="CellType", cex=.5)

# 4.3.4 Save defalt t-SNE
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, 
        file="/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/scenicOptions.Rds")

print(tsneFileName(scenicOptions))

dev.off()

#Regulators for known cell types or clusters
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$SingleR_ENCODE),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")



