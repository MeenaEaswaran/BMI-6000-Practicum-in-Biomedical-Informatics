# Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent")

library(shiny)
library(AUCell)
library(SCENIC)
library(tidyverse)
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


#Read SCENIC related RDS
scenicOptions<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/scenicOptions.Rds")

cellInfo <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/SCENIC_LSCC_2percent/int/cellInfo.Rds")


#Regulators for known cell types or clusters in SCENIC tutorial 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$SingleR_ENCODE),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")



auc <- getAUC(regulonAUC)                     # regulons x cells
stopifnot(all(colnames(auc) %in% rownames(cellInfo)))
cellInfo <- cellInfo[colnames(auc), , drop=FALSE]  # align


#Regulon activity by Condition (e.g., Normal vs Tumor)
regulonActivity_byCond <- sapply(split(colnames(auc), cellInfo$Condition),
                                 function(cells) rowMeans(auc[, cells, drop=FALSE]))
regulonActivity_byCond_scaled <- t(scale(t(regulonActivity_byCond), center=TRUE, scale=TRUE))

ComplexHeatmap::Heatmap(regulonActivity_byCond_scaled,
                        name = "Regulon z-score",
                        column_title = "Condition")

#2) Regulon activity by Condition × Cell type 
grp <- interaction(cellInfo$Condition, cellInfo$SingleR_ENCODE, drop=TRUE)  # e.g., Tumor.T-cells
regulonActivity_byCondCell <- sapply(split(colnames(auc), grp),
                                     function(cells) rowMeans(auc[, cells, drop=FALSE]))
regulonActivity_byCondCell_scaled <- t(scale(t(regulonActivity_byCondCell), center=TRUE, scale=TRUE))

ComplexHeatmap::Heatmap(regulonActivity_byCondCell_scaled,
                        name = "Regulon z-score",
                        column_title = "Condition × Cell type",
                        cluster_columns = FALSE)

# what we want :) 

library(ComplexHeatmap)
library(grid)
library(circlize)
library(RColorBrewer)

# 1) Build the matrix as you already do
grp <- interaction(cellInfo$Condition, cellInfo$SingleR_ENCODE, drop = TRUE)  # e.g., "Tumor.T-cells"
regulonActivity_byCondCell <- sapply(split(colnames(auc), grp),
                                     function(cells) rowMeans(auc[, cells, drop = FALSE]))
regulonActivity_byCondCell_scaled <- t(scale(t(regulonActivity_byCondCell), center = TRUE, scale = TRUE))

# 2) Parse "Condition.Celltype" from column names and set desired orders
cell_order <- c("Epithelial cells", "Fibroblasts", "Myeloid cells",
                "NK cells", "T-cells", "B-cells", "Endothelial cells")  # <- adjust to match your labels exactly

cn <- colnames(regulonActivity_byCondCell_scaled)
parts <- do.call(rbind, strsplit(cn, "\\.", fixed = FALSE))
cond  <- factor(parts[, 1], levels = c("Normal", "Tumor"))
ctype <- factor(parts[, 2], levels = cell_order)

ord <- order(cond, ctype, na.last = TRUE)
mat <- regulonActivity_byCondCell_scaled[, ord, drop = FALSE]
cond_ord  <- cond[ord]
ctype_ord <- ctype[ord]

## --- Colors for Cell type in the NEW order you requested ---
ctype_cols <- setNames(
  c("#332288", "#1CA71C", "#DDCC77", "#999933", "#882255", "#117733", "#654522"),
  cell_order
)

## Top annotation (Condition + Cell type)
cond_cols <- c(Normal = "grey", Tumor = "red")
ta <- HeatmapAnnotation(
  Condition = cond_ord,
  `Cell type` = ctype_ord,
  col = list(Condition = cond_cols, `Cell type` = ctype_cols),
  annotation_name_side = "left"
)

## Draw (unchanged)
Heatmap(
  mat,
  name = "Regulon z-score",
  #column_title = "Condition × Cell type",
  cluster_columns = FALSE,
  column_split   = cond_ord,
  column_gap     = unit(4, "mm"),
  top_annotation = ta,
  show_column_dend = FALSE
) |> draw(heatmap_legend_side = "right", annotation_legend_side = "right")


#Another way
library(dplyr)
library(tidyr)
library(ggplot2)

col_meta <- data.frame(
  colname   = colnames(mat),
  Condition = factor(as.character(cond_ord), levels = c("Normal","Tumor")),
  cell_type = factor(as.character(ctype_ord), levels = cell_order),
  stringsAsFactors = FALSE
)

df_plot <- as.data.frame(mat) |>
  tibble::rownames_to_column("regulon") |>
  tidyr::pivot_longer(-regulon, names_to = "colname", values_to = "Z") |>
  dplyr::left_join(col_meta, by = "colname") |>
  dplyr::select(Condition, cell_type, regulon, Z)


regulon_order <- df_plot |>
  group_by(regulon) |>
  summarise(meanZ = mean(Z, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(meanZ)) |>
  pull(regulon)

df_plot <- df_plot |>
  mutate(cell_type = factor(cell_type, levels = rev(cell_order)))

z_label_threshold <- 0   # set to 2.5 to mimic seaborn style

df_plot <- df_plot |>
  mutate(label = ifelse(!is.na(Z) & abs(Z) >= z_label_threshold,
                        sprintf("%.1f", Z), ""))


p <- ggplot(df_plot, aes(x = regulon, y = cell_type, fill = Z)) +
  geom_tile(color = "grey80", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, na.value = "white",
    name = "Z-score"   # legend title
  ) +
  facet_wrap(~ Condition, ncol = 2) +
  coord_fixed() +
  labs(x = "Regulons", y = "Cell types") +  # axis labels
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 12),        # Normal/Tumor facet labels
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(face = "bold", size = 10),
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, colour = NA)
  )

p

ggsave("10_anotherheatmap_condition_celltype_regulon_2percent.png", p, width = 20, height = 6, dpi = 300)


#To view on SCope
### Load data
loom_path <- "~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/10_scenic_input.loom"

loom <- open_loom(loom_path)
exprMat <- get_dgem(loom)
dgem <- exprMat
head(colnames(dgem))


scenicOptions@fileNames$output["loomFile",] <- "~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/10_scenic_ouput.loom"
export2loom(scenicOptions, exprMat)
