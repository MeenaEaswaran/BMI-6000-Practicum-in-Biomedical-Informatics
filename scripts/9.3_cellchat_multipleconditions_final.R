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

#Load CellChat object of each dataset and merge them together
#Read Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
cellchat.normal<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/9_cellchatouput_normalsamplesonly.rds")
cellchat.tumor<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/9_cellchatouput_tumorsamplesonly.rds")

object.list <- list(N = cellchat.normal, T = cellchat.tumor)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

#Export the merged CellChat object and the list of the two separate objects for later use
save(object.list, file = "9_cellchat_object.list_LSCC_normal_tumor.RData")
save(cellchat, file = "9_cellchat_merged__LSCC_normal_tumor.RData")
########################################################################################################

#Part I: Identify altered interactions and cell populations

#CellChat employs a top-down approach, i.e., starting with the big picture and then refining it in a greater detail on the signaling mechanisms, to identify signaling changes at different levels, including altered interactions, cell populations, signaling pathways and ligand-receptor pairs. First, CellChat starts with a big picture to answer the following questions:
#Whether the cell-cell communication is enhanced or not
#The interaction between which cell types is significantly changed
#How the major sources and targets change from one condition to another

#Compare the total number of interactions and interaction strength
#To answer the question on whether the cell-cell communication is enhanced or not, CellChat compares the total number of interactions and interaction strength of the inferred cell-cell communication networks from different biological conditions.
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

ggsave("9_cellchatDB_TvsN_interactions_strength_barplot.jpg", units="in", width=10, height=8, dpi=300)

#Compare the number of interactions and interaction strength among different cell populations
#To identify the interaction between which cell populations showing significant changes, CellChat compares the number of interactions and interaction strength among different cell populations using a circle plot with differential interactions (option A), a heatmap with differential interactions (option B) and two circle plots with the number of interactions or interaction strength per dataset (option C). Alternatively, users can examine the differential number of interactions or interaction strength among coarse cell types by aggregating the cell-cell communication based on the defined cell groups (option D).

#(A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#(B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets.
#CellChat can also show differential number of interactions or interaction strength in greater details using a heatmap. The top colored bar plot represents the sum of each column of the absolute values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of each row of the absolute values (outgoing signaling). Therefore, the bar height indicates the degree of change in terms of the number of interactions or interaction strength between the two conditions. In the colorbar, red (or blue ) represents increased (or decreased) signaling in the second dataset compared to the first one.

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#(C) Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets.

#The above differential network analysis only works for pairwise datasets. If there are more datasets for comparison, CellChat can directly show the number of interactions or interaction strength between any two cell populations in each dataset.

# To better control the node size and edge weights of the inferred networks across different datasets, CellChat computes the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#Compare the major sources and targets in a 2D space
#(A) Identify cell populations with significant changes in sending or receiving signals

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
ggsave("9_cellchatDB_TvsN_Identifycellpopulationswithsignificantchangesinsendingorreceivingsignals.jpg", units="in", width=10, height=8, dpi=300)

#we see that all cell types except b-cells have some  changes
#(B) Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epithelial cells")
#> Visualizing differential outgoing and incoming signaling changes from N to T
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Fibroblasts")
#> Visualizing differential outgoing and incoming signaling changes from N to T
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Myeloid cells")
#> Visualizing differential outgoing and incoming signaling changes from N to T
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NK cells")
#> Visualizing differential outgoing and incoming signaling changes from N to T
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "B-cells")
#> Visualizing differential outgoing and incoming signaling changes from N to T
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T-cells")
#> Visualizing differential outgoing and incoming signaling changes from N to T
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg7 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endothelial cells")
#> Visualizing differential outgoing and incoming signaling changes from N to T
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2,gg3, gg4, gg5, gg6, gg7))

#########################################################################################

#Part II: Identify altered signaling with distinct network architecture and interaction strength
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)

#########################################################################################

#Identify altered signaling with distinct interaction strength

#By comparing the information flow/interaction strength of each signaling pathway, CellChat identifies signaling pathways that: (i) turn off, (ii) decrease, (iii) turn on, or (iv) increase, by changing their information flow at one condition as compared to another condition. Identify the altered signaling pathways or ligand-receptor pairs based on the overall information flow by following option A, and based on the outgoing (or incoming) signaling patterns by following option B.

#(A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2
ggsave("9_cellchatDB_TvsN_Compare the overall information flow of each signaling pathway or ligand-receptor pair.jpg", units="in", width=10, height=10, dpi=300)


#(B) Compare outgoing (or incoming) signaling patterns associated with each cell population
#The above analysis summarize the information from the outgoing and incoming signaling together. CellChat can also compare the outgoing (or incoming) signaling pattern between two datasets, allowing to identify signaling pathways/ligand-receptors that exhibit different signaling patterns. We can combine all the identified signaling pathways from different datasets and thus compare them side by side, including outgoing signaling, incoming signaling and overall signaling by aggregating outgoing and incoming signaling together.

#CellChat uses heatmap plot to show the contribution of signals (signaling pathways or ligand-receptor pairs) to cell groups in terms of outgoing or incoming signaling. In this heatmap, color bar represents the relative signaling strength of a signaling pathway across cell groups (Note that values are row-scaled). The top colored bar plot shows the total signaling strength of a cell group by summarizing all signaling pathways displayed in the heatmap. The right grey bar plot shows the total signaling strength of a signaling pathway by summarizing all cell groups displayed in the heatmap.

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(1, "cm"), height = unit(20, "cm")) #SAVE

#incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(1, "cm"), height = unit(20, "cm")) #SAVE

########################################################################################################
#Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities
#use for reference in the plots below to be saved
#> levels(cellchat@idents)
#[1] "B-cells"           "Endothelial cells" "Epithelial cells"  "Fibroblasts"       "Myeloid cells"     "NK cells"         
#[7] "T-cells" 

#B cells 
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:7),  comparison = c(1, 2), angle.x = 45, thresh = 0.01, remove.isolate = TRUE)

#Endothelial cells 
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1, 3:7),  comparison = c(1, 2), angle.x = 45, thresh = 0.01, remove.isolate = TRUE)

#epithelial cells 
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2, 4:7),  comparison = c(1, 2), angle.x = 45, thresh = 0.01, remove.isolate = TRUE)

#fibroblasts 
jpeg("9_cellchat_Identify dysfunctional signaling by comparing the communication probabities_tvsN_fibroblasts.jpeg",
    width = 2800,          # width in pixels
    height = 4500,         # increase height (in pixels)
    res = 300)             # resolution (300 dpi is print-quality)

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3, 5:7),  comparison = c(1, 2), angle.x = 45, thresh = 0.01, remove.isolate = TRUE)

dev.off()

#myeloid cells 
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4, 6:7),  comparison = c(1, 2), angle.x = 45, thresh = 0.01, remove.isolate = TRUE)

#NK cells 
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:5, 7),  comparison = c(1, 2), angle.x = 45, thresh = 0.01, remove.isolate = TRUE)

#T cells 
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45, thresh = 0.01, remove.isolate = TRUE)

################################################################################################################

#Identify dysfunctional signaling by using differential expression analysis
#The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups. Alternative, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential expression analysis (DEA). Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, and then obtain the upgulated and down-regulated signaling based on the fold change of ligands in the sender cells and receptors in the receiver cells.

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "T"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in Tumor
net.up <- subsetCommunication(cellchat, net = net, datasets = "T",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in Normal, i.e.,downregulated in Tumor
net.down <- subsetCommunication(cellchat, net = net, datasets = "N",ligand.logFC = -0.05, receptor.logFC = NULL)
################################################################################################################

#Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs
#to show significant LR pairs only up/downregulated in tumor samples via bubble plot or chord plot
#we are showing only bubble plot

#B-cells
#bubble plot 
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(2:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(2:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2 #save

######################################################################################
#Endothelial cells 
#bubble plot 
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 2, targets.use = c(1, 3:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 2, targets.use = c(1, 3:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2 #save
######################################################################################
#epithelial cells 
#bubble plot 
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(1:2, 4:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 3, targets.use = c(1:2, 4:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2 #save
######################################################################################
#fibroblasts 
#bubble plot 
jpeg("9_cellchat_significant up-regulated and down-regulated signaling ligand-receptor pairs_bubbleplot_fibroblasts.jpeg",
    width = 2800,          # width in pixels
    height = 4500,         # increase height (in pixels)
    res = 300)    

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(1:3, 5:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(1:3, 5:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2 #save

dev.off()

######################################################################################
#myeloid cells 
#bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 5, targets.use = c(1:4, 6:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 5, targets.use = c(1:4, 6:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2 #save
######################################################################################
#NK cells 
#bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 6, targets.use = c(1:5, 7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 6, targets.use = c(1:5, 7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2 #save

######################################################################################
#T cells 
#bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 7, targets.use = c(1:6), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 7, targets.use = c(1:6), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2 #save


