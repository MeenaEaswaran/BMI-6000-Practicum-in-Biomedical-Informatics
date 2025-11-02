# Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

#Run for the first time
#BiocManager::install("glmGamPoi") #needed for sctransform faster run

#To fix errors
#install.packages("devtools")
#devtools::install_github("satijalab/seurat", ref = "fix/v.5.3.1")

#install UMAP if not already installed
reticulate::py_install(packages='umap-learn')  

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
library(sctransform)
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
options(future.globals.maxSize = 1e9)
library(CellChat)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(NMF)
library(compiler)
library(parallel)
library(future.callr)
#library(DoubletFinder)
options(stringsAsFactors = FALSE)
library(reticulate)

#To split original harmony integrated dataset by condition-normal and tumor- SKIP IF NOT REQUIRED

#Load integrated object 
objs <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithfinalcelltypeannotations_seurat_30dims.rds")

DefaultAssay(objs) # if needed: DefaultAssay(objs) <- "SCT"

#To run cellchat on tumor samples and normal samples separately
objs$Condition <- factor(objs$Condition, levels = c("Normal","Tumor"))
seu_list <- SplitObject(objs, split.by = "Condition")

saveRDS(seu_list$Normal, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/9_cellchat_seurat_Normal_SCT.rds")
saveRDS(seu_list$Tumor, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/9_cellchat_seurat_Tumor_SCT.rds")


# Read Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
objs_normal<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/9_cellchat_seurat_Normal_SCT.rds")

DefaultAssay(objs_normal) #check if SCT

# Extract the normalized data
data_input <- GetAssayData(objs_normal, assay = "SCT", layer = "data")

#Create CellChat object
cellchat <- createCellChat(object = data_input, meta = objs_normal@meta.data, group.by = "SingleR_ENCODE") #we are grouping by cell types similar to the original paper

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

#Before users can use CellChat to infer cell-to-cell communication, they need to set up a ligand-receptor interaction database and identify overexpressed ligands or receptors.
#CellChatDB is a hand-curated database of human and mouse literature supporting ligand-receptor interactions. CellChatDB v2 contains ~3,300 validated molecular interactions, including ~40% secretory autocrine/paracrine signaling interactions, ~17% extracellular matrix (ECM)-receptor interactions, ~13% cell- Cell contact interactions and approximately 30% of non-protein signaling. Compared with CellChatDB v1, CellChatDB v2 adds more than 1000 protein and non-protein interactions, such as metabolism and synaptic signaling. It should be noted that for gene-related molecules that are not directly measured in single-cell RNA sequencing, CellChat v2 estimates the expression of ligands and receptors, using key mediators or enzymes of these molecules, in order to potentially pass through non-protein mediated communication.
#Users can update CellChatDB by adding their own curated ligand-receptor pairs.

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
ggsave("9_cellchatDB_human_database_normalsamples.jpg", units="in", width=10, height=8, dpi=300)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat

#Preprocessing the expression data for cell-cell communication analysis
#To infer the cell state-specific communications, CellChat identifies over-expressed ligands or receptors in one cell group and then identifies over-expressed ligand-receptor interactions if either ligand or receptor are over-expressed.

#We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. One might be concerned about the possible artifact introduced by this diffusion process, however, it will only introduce very weak communications. By default CellChat uses the raw data (i.e., object@data.signaling) instead of the projected data. To use the projected data, users should run the function projectData before running computeCommunProb, and then set raw.use = FALSE when running computeCommunProb.

#projectdata not supported in seurat V5. Its smoothData as below.also future is not supported.

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
# do parallel
#future::plan("multiprocess", workers = 4) #not supported anymore in V5 Seurat
#https://github.com/satijalab/seurat/issues/9385

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)

#Part II: Inference of cell-cell communication network
#CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value and peforming a permutation test. CellChat models the probability of cell-cell communication by integrating gene expression with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.

#The number of inferred ligand-receptor pairs clearly depends on the method for calculating the average gene expression per cell group. By default, CellChat uses a statistically robust mean method called ‘trimean’, which produces fewer interactions than other methods. However, we find that CellChat performs well at predicting stronger interactions, which is very helpful for narrowing down on interactions for further experimental validations.

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

#CellChat can also visualize aggregated intercellular communication networks. For example, use a circle plot to show the number of interactions or the total interaction strength (weight) between any two groups of cells.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions:Normal samples")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength:Normal samples")


#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
mat <- cellchat@net$weight

# Layout: 3 rows x 4 cols; expand outer & inner margins; allow drawing beyond plot area
par(mfrow = c(3, 4),
    mar   = c(1.2, 1.2, 3.2, 1.2),  # plot margins (b, l, t, r) per panel
    oma   = c(1, 1, 1, 1),          # outer margins around the full grid
    xpd   = NA,                     # allow labels to extend past plot region
    cex   = 1.1)                    # global size multiplier for base plotting

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight     = groupSize,
    weight.scale      = TRUE,
    edge.weight.max   = max(mat),
    label.edge        = FALSE,
    title.name        = rownames(mat)[i],
    vertex.label.cex  = 1.1,   # bigger node labels
    arrow.size        = 0.5    # slightly larger arrows (optional)
  )
}

#Part III: Visualization of cell-cell communication network
#Upon infering the cell-cell communication network, CellChat provides various functionality for further data exploration, analysis, and visualization.

#It provides several ways for visualizing cell-cell communication network, including hierarchical plot, circle plot, Chord diagram, and bubble plot.

#It provides an easy-to-use tool for extracting and visualizing high-order information of the inferred networks. For example, it allows ready prediction of major signaling inputs and outputs for cell populations and how these populations and signals coordinate together for functions.

#It can quantitatively characterize and compare the inferred cell-cell communication networks using an integrated approach by combining social network analysis, pattern recognition, and manifold learning approaches.

#Here we can take the input of a signal path as an example or all pathways. All signaling pathways showing significant communication can be accessed via cellchat@netP$pathways.
##########
cellchat@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)

#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#CellChat can also show all the significant interactions mediated by L-R pairs and signaling pathways, and interactions provided by users from some cell groups to other cell groups using the function netVisual_bubble (option A) and netVisual_chord_gene (option B).

## (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = levels(cellchat@idents), remove.isolate = TRUE, thresh = 0.01 )

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pathways.show.all <- cellchat@netP$pathways
pathways.show.all #58 pathways here

#show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show.all)


#use for reference in the plots below to be saved
#> levels(cellchat@idents)
#[1] "B-cells"           "Endothelial cells" "Epithelial cells"  "Fibroblasts"       "Myeloid cells"     "NK cells"         
#[7] "T-cells" 

#B cells 
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:7), signaling = pathways.show.all, remove.isolate = TRUE, thresh = 0.01)

#Endothelial cells 
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1, 3:7), signaling = pathways.show.all, remove.isolate = TRUE, thresh = 0.01)

#epithelial cells 
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2, 4:7), signaling = pathways.show.all, remove.isolate = TRUE, thresh = 0.01)

#fibroblasts 
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3, 5:7), signaling = pathways.show.all, remove.isolate = TRUE, thresh = 0.01)

#myeloid cells 
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4, 6:7), signaling = pathways.show.all, remove.isolate = TRUE, thresh = 0.01)

#NK cells 
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:5, 7), signaling = pathways.show.all, remove.isolate = TRUE, thresh = 0.01)

#T cells 
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6), signaling = pathways.show.all, remove.isolate = TRUE, thresh = 0.01)


#Systems analysis of cell-cell communication network
#CellChat allows ready identification of dominant senders, receivers, mediators and influencers in the intercellular communication network by computing several network centrality measures for each cell group. 
#Compute and visualize the network centrality scores
ptm = Sys.time()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all, width = 8, height = 2.5, font.size = 10)

#visualize dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1

#Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
#We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups. In this heatmap, colobar represents the relative signaling strength of a signaling pathway across cell groups (NB: values are row-scaled). The top colored bar plot shows the total signaling strength of a cell group by summarizing all signaling pathways displayed in the heatmap. The right grey bar plot shows the total signaling strength of a signaling pathway by summarizing all cell groups displayed in the heatmap.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# draw with larger height
draw(ht1 + ht2, height = unit(20, "cm")) #SAVE

#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together

#Identify and visualize outgoing communication pattern of secreting cells

library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing") #SAVE
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 4.
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns) #SAVE
 
# river plot
netAnalysis_river(cellchat, pattern = "outgoing") #SAVE
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing") # dot plot #SAVE


#Identify and visualize incoming communication pattern of target cells
#Shows how target cells (i.e. cells that are signal receivers) coordinate with each other and how they coordinate with certain signaling pathways to respond to received signals.
selectK(cellchat, pattern = "incoming") #SAVE
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 5.
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns) #SAVE

# river plot
netAnalysis_river(cellchat, pattern = "incoming")  #SAVE

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# **Manifold and classification learning analysis of signaling networks**
#Further, CellChat is able to quantify the similarity between all significant signaling pathways and then group them based on their cellular communication network similarity. Grouping can be done either based on the functional or structural similarity.

#Functional similarity: High degree of functional similarity indicates major senders and receivers are similar, and it can be interpreted as the two signaling pathways or two ligand-receptor pairs exhibit similar and/or redundant roles. The functional similarity analysis requires the same cell population composition between two datasets.

#Structural similarity: A structural similarity was used to compare their signaling network structure, without considering the similarity of senders and receivers.


#do below to avoid pythonb UMAP error
library(reticulate)
use_condaenv("base", conda = "/opt/anaconda3/bin/conda", required = TRUE) 

py_config()

reticulate::py_install(packages = 'umap-learn', envname = "base", pip = TRUE)

py_module_available("umap")

#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


##Identification of signal groups based on structural similarities.
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")

#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")

#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


#Save the CellChat object
saveRDS(cellchat, file = "9_cellchat_normalsamples.rds")
