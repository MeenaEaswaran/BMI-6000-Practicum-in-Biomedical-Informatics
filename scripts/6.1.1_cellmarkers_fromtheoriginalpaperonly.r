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

#Quick type check with canonical markers in the dataset first and also other markers mined from literature
#Approach: markers similar to the one used in the original paper
DefaultAssay(objs) <- "SCT"

#Map aliases → canonical symbols
marker_panels <- list(
    Epithelial = c("KRT15","KRT18","KRT19","EPCAM", "KRT5"),
    Fibroblast = c("alpha-SMA","FAP","S100A4"), # alpha-SMA(ACTA2)
    Myeloid    = c("CD33","CD68","CD1E","LYZ","LAMP3"), # CD33 = SIGLEC3
    NK         = c("CD56","CD16","NKP46","NKP30"),  # CD56=NCAM1; CD16=FCGR3A; NKP46=NCR1; NKP30=NCR3
    T_cells    = c("CD2", "CD3D","CD3E","CD3G"), 
    B_cells    = c("CD19","CD79A","CD79B"),
    Endothelial= c("VEGFR", "TEK","CD54")  # VEGFR -> KDR/FLT1; CD54 -> ICAM1
     )

#Keep only genes present in the object & report any missing
present <- rownames(objs)
missing_report <- lapply(marker_panels, function(v) setdiff(v, intersect(v, present)))
marker_panels <- lapply(marker_panels, function(v) intersect(v, present))
print(missing_report)  # look here to see which aliases weren't found

#now check with aliases
#Map aliases → canonical symbols
marker_panels <- list(
   Epithelial = c("KRT15","KRT18","KRT19","EPCAM", "KRT5"),
   Fibroblast = c("ACTA2","FAP","S100A4"), # alpha-SMA(ACTA2)
   Myeloid    = c("CD33","CD68","CD1E","LYZ","LAMP3"), # CD33 = SIGLEC3
   NK         = c("NCAM1","FCGR3A","NCR1","NCR3"),   # CD56=NCAM1; CD16=FCGR3A; NKP46=NCR1; NKP30=NCR3
   T_cells    = c("CD2","CD3D","CD3E","CD3G"),
   B_cells    = c("CD19","CD79A","CD79B"),
   Endothelial= c("KDR","FLT1","TEK","ICAM1")       # VEGFR -> KDR/FLT1; CD54 -> ICAM1
  )

#Keep only genes present in the object & report any missing
present <- rownames(objs)
missing_report <- lapply(marker_panels, function(v) setdiff(v, intersect(v, present)))
marker_panels <- lapply(marker_panels, function(v) intersect(v, present))
print(missing_report)  # look here to see which aliases weren't found
#all not present 
#proceed


#DotPlot: quick multi-panel overview by cluster 
panel_features <- unique(unlist(marker_panels))
DotPlot(objs, features = panel_features, group.by = "seurat_clusters", assay = "SCT") +
  RotatedAxis() + ggtitle("Original paper marker overview by Harmony integrated clusters")

#Save image
ggsave(
  filename = "6_canonicalmarkers_originalpaperonly_harmonyintegratedUMAPclusters_rawimage.jpg",
  units = "in", width = 18, height = 8, dpi = 300
)

#for a more detailed plot ordered by clusters and marker types
library(Seurat)
library(ggplot2)
library(scales)

#Prepare ordering & spans 
panel_order <- names(marker_panels)

feat_order <- unique(unlist(marker_panels[panel_order], use.names = FALSE))
feat_order <- feat_order[feat_order %in% rownames(objs[["SCT"]])]

panel_map <- data.frame(
  feature = unlist(marker_panels[panel_order], use.names = FALSE),
  panel   = rep(panel_order, lengths(marker_panels[panel_order])),
  stringsAsFactors = FALSE
)
panel_map <- subset(panel_map, feature %in% feat_order)
panel_map <- panel_map[!duplicated(panel_map$feature), ]

xpos <- match(panel_map$feature, feat_order)
xmin <- tapply(xpos, panel_map$panel, min) - 0.5
xmax <- tapply(xpos, panel_map$panel, max) + 0.5
panel_spans <- data.frame(panel = names(xmin),
                          xmin  = as.numeric(xmin),
                          xmax  = as.numeric(xmax))
panel_spans <- panel_spans[order(panel_spans$xmin), ]

#DISTINCT colors for top bars (HCL space) 
cb_palette <- function(n) {
  base <- c(
    # Okabe–Ito (7)
        "#332288", "#1CA71C", "#DDCC77", "#999933","#882255","#117733","#654522"
  )
  if (n <= length(base)) base[seq_len(n)] else
    grDevices::hcl(h = seq(15, 375, length.out = n + 1)[-1], c = 90, l = 60)
}
panel_cols <- setNames(cb_palette(nrow(panel_spans)), panel_spans$panel)

# Customize panel names
panel_rename <- c(
    Epithelial = "Epithelial cells",
    Fibroblast = "Fibroblast",
    NK = "NK cells",
   Myeloid = "Myeloid cells",
   T_cells = "T cells",
  B_cells = "B cells",
  Endothelial = "Endothelial cells"
  )

#Fallback: if a panel isn't listed above, use a prettified default
pretty_fallback <- function(x) gsub("_", " ", x)
panel_labels <- ifelse(is.na(panel_rename[panel_spans$panel]),
                       pretty_fallback(panel_spans$panel),
                       panel_rename[panel_spans$panel])

#Lines/strips geometry 
vlines_x_inner <- head(panel_spans$xmax, -1)
vline_left     <- head(panel_spans$xmin, 1)
vline_right    <- tail(panel_spans$xmax, 1)
sep_col        <- "grey40"

p <- DotPlot(objs, features = feat_order, group.by = "seurat_clusters", assay = "SCT") +
  RotatedAxis() +
  ggtitle("Original paper marker overview by Harmony integrated clusters") +
  labs(x = "Genes (Features)", y = "Cluster identity")

df <- p$data
n_rows   <- length(levels(df$id))
bar_ymin <- n_rows + 0.70
bar_ymax <- n_rows + 1.30

bars_df <- transform(panel_spans, ymin = bar_ymin, ymax = bar_ymax, fill = panel)

vseg_df <- data.frame(
  x    = c(vline_left, vlines_x_inner, vline_right),
  y    = 0.5,
  yend = bar_ymin       # touch the strip (no gap)
)
hseg_df <- data.frame(x = vline_left, xend = vline_right, y = 0.5)

#Set symmetric red–white–blue scale for Average Expression 
L <- max(abs(range(df$avg.exp.scaled, na.rm = TRUE)))  # symmetric limits

#plot the final dot plot
p2 <- p +
  geom_rect(data = bars_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE) +
  scale_fill_manual(
    values = panel_cols,
    breaks = panel_spans$panel,     # legend order = panel order
    labels = panel_labels,          # custom legend labels
    name   = "Marker panel"
  ) +
  geom_segment(data = vseg_df,
               aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, color = sep_col, linewidth = 0.4) +
  geom_segment(data = hseg_df,
               aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, color = sep_col, linewidth = 0.4) +
  scale_color_gradient2(
    low = "#2C7FB8", mid = "#FFFFFF", high = "#D7301F",  # blue–white–red
    midpoint = 0, limits = c(-L, L),
    breaks = pretty(c(-L, L), n = 5), name = "Average Expression"
  ) +
  scale_y_discrete(expand = expansion(add = c(0.60, 1.40))) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin  = margin(t = 26, r = 6, b = 6, l = 6),
    axis.line    = element_blank(),
    panel.border = element_blank()
  )
print(p2)

# Optional: hide the legend for the panel colors
# p2 + guides(fill = "none")

# Save image
ggsave(
  filename = "6_canonicalmarkers_originalpaper_harmonyintegratedUMAPclusters_groupedbymarkers.jpg",
  units = "in", width = 18, height = 8, dpi = 300
)

#by Condition instead of cluster (preliminary)
DotPlot(objs, features = panel_features, group.by = "Condition", assay = "SCT") + RotatedAxis()

# Save image
ggsave(
  filename = "6_canonicalmarkers_originalpaper_harmonyintegratedUMAPclusters_rawimagebycondition.jpg",
  units = "in", width = 18, height = 8, dpi = 300
)

#more detailed dot plot
#Prepare ordering & spans 
panel_order <- names(marker_panels)

feat_order <- unique(unlist(marker_panels[panel_order], use.names = FALSE))
feat_order <- feat_order[feat_order %in% rownames(objs[["SCT"]])]

panel_map <- data.frame(
  feature = unlist(marker_panels[panel_order], use.names = FALSE),
  panel   = rep(panel_order, lengths(marker_panels[panel_order])),
  stringsAsFactors = FALSE
)
panel_map <- subset(panel_map, feature %in% feat_order)
panel_map <- panel_map[!duplicated(panel_map$feature), ]

xpos <- match(panel_map$feature, feat_order)
xmin <- tapply(xpos, panel_map$panel, min) - 0.5
xmax <- tapply(xpos, panel_map$panel, max) + 0.5
panel_spans <- data.frame(panel = names(xmin),
                          xmin  = as.numeric(xmin),
                          xmax  = as.numeric(xmax))
panel_spans <- panel_spans[order(panel_spans$xmin), ]

#DISTINCT colors for top bars (HCL space) 
cb_palette <- function(n) {
  base <- c(
    # Okabe–Ito (7)
    "#332288", "#1CA71C", "#DDCC77", "#999933","#882255","#117733","#654522"
  )
  if (n <= length(base)) base[seq_len(n)] else
    grDevices::hcl(h = seq(15, 375, length.out = n + 1)[-1], c = 90, l = 60)
}
panel_cols <- setNames(cb_palette(nrow(panel_spans)), panel_spans$panel)

# Customize panel names
panel_rename <- c(
  Epithelial = "Epithelial cells",
  Fibroblast = "Fibroblast",
  NK = "NK cells",
  Myeloid = "Myeloid cells",
  T_cells = "T cells",
  B_cells = "B cells",
  Endothelial = "Endothelial cells"
)

#Fallback: if a panel isn't listed above, use a prettified default
pretty_fallback <- function(x) gsub("_", " ", x)
panel_labels <- ifelse(is.na(panel_rename[panel_spans$panel]),
                       pretty_fallback(panel_spans$panel),
                       panel_rename[panel_spans$panel])

#Lines/strips geometry 
vlines_x_inner <- head(panel_spans$xmax, -1)
vline_left     <- head(panel_spans$xmin, 1)
vline_right    <- tail(panel_spans$xmax, 1)
sep_col        <- "grey40"

#Base DotPlot: group by Condition
p3 <- DotPlot(objs, features = feat_order, group.by = "Condition", assay = "SCT") +
  RotatedAxis() +
  ggtitle("Original paper marker overview by Condition") +
  labs(x = "Genes (Features)", y = "Tissue type")

df <- p3$data
n_rows <- length(levels(df$id))

# ---- Slim strip with zero bottom padding on y-scale ----
bar_gap    <- 0.46         # distance from top row center to bottom of strip
bar_height <- 0.16         # thickness of the colored strip
bar_ymin   <- n_rows + bar_gap
bar_ymax   <- bar_ymin + bar_height

# just enough top padding so the strip is fully visible; bottom padding = 0
top_add <- max(0, bar_ymax - (n_rows + 0.5)) + 0.02

bars_df <- transform(panel_spans, ymin = bar_ymin, ymax = bar_ymax, fill = panel)

# make verticals TOUCH the strip (no gap)
vseg_df <- data.frame(
  x    = c(vline_left, vlines_x_inner, vline_right),
  y    = 0.5,
  yend = bar_ymin
)

#bottom horizontal axis line exactly on panel boundary
hseg_df <- data.frame(x = vline_left, xend = vline_right, y = 0.5)

#symmetric red–white–blue scale for Average Expression 
L <- max(abs(range(df$avg.exp.scaled, na.rm = TRUE)))  

#plot the final dot plot
p4 <- p3 +
  # Colored bar per marker group
  geom_rect(data = bars_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE) +
  scale_fill_manual(
    values = panel_cols,
    breaks = panel_spans$panel,     # legend order = panel order
    labels = panel_labels,          # custom legend labels
    name   = "Marker panel"
  ) +
    # Boxed verticals + bottom horizontal line (all same color)
  geom_segment(data = vseg_df,
               aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, color = sep_col, linewidth = 0.4) +
  geom_segment(data = hseg_df,
               aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, color = sep_col, linewidth = 0.4) +
    scale_color_gradient2(
    low = "#2C7FB8", mid = "#FFFFFF", high = "#D7301F",  # blue–white–red
    midpoint = 0, limits = c(-L, L),
    breaks = pretty(c(-L, L), n = 5), name = "Average Expression"
  ) +
      # ZERO bottom padding; only the minimal top padding required for the strip
  scale_y_discrete(expand = expansion(mult = c(0, 0), add = c(0, top_add))) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin  = margin(t = 16, r = 6, b = 6, l = 6),
    axis.line    = element_blank(),
    panel.border = element_blank()
  )

print(p4)

# Optional: remove legend
# p4 + guides(fill = "none")

# Save image
ggsave(
  filename = "6_canonicalmarkers_originalpaper_harmonyintegrated_basedoncondition.jpg",
  units = "in", width = 18, height = 8, dpi = 300
)
##################################################################################
##################################################################################
##################################################################################

#The images below are to overlay individual clusters on the harmony integrated umap with cluster number
#draw the original UMAP first
p_umap <- DimPlot(
  objs,
  reduction = "umap.harmony",
  group.by  = "seurat_clusters",
  label     = TRUE, 
  repel = TRUE,
  pt.size   = 0.3
) + ggtitle("Harmony UMAP — clusters")
p_umap
# Save image
ggsave(
  filename = "6_harmonyintegrated_UMAP_with cluster names.jpg",
  units = "in", width = 8, height = 8, dpi = 300
)
##################################################################################
# flatten to unique gene list :Epithelial (LSCC originates from epithelial cells)
feats_epithelial <- intersect(marker_panels$Epithelial, rownames(objs))
p_epi <- FeaturePlot(
  objs,
  features  = feats_epithelial,
  reduction = "umap.harmony",
  order     = TRUE,
  pt.size   = 0.3,  
  label = TRUE,
  repel=TRUE
)
p_epi

# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplot_withepithelialmarkers.jpg",
  units = "in", width = 10, height = 8, dpi = 300
)

#split by condition
p_epi_split <-FeaturePlot(
  objs,
  features  = feats_epithelial,
  reduction = "umap.harmony",
  split.by = "Condition",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE,
  keep.scale = "feature", 
  cols = c("grey",
           "red")
) 

p_epi_split 

# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplotBYCONDITION_withepithelialmarkers.jpg",
  units = "in", width = 12, height = 12, dpi = 300
)
##################################################################################
# flatten to unique gene list: fibroblast
feats_fibo <- intersect(marker_panels$Fibroblast, rownames(objs))
pfibo <- FeaturePlot(
  objs,
  features  = feats_fibo,
  reduction = "umap.harmony",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE
)
pfibo
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplot_withfibroblast_markers.jpg",
  units = "in", width = 10, height = 8, dpi = 300
)
#split by condition
pfibo_split <-FeaturePlot(
  objs,
  features  = feats_fibo,
  reduction = "umap.harmony",
  split.by = "Condition",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE,
  keep.scale = "feature", 
  cols = c("grey",
           "red")
) 
pfibo_split 
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplotBYCONDITION_withfibroblast_markers.jpg",
  units = "in", width = 12, height = 12, dpi = 300
)
##################################################################################
# flatten to unique gene list: Myeloid
feats_myl <- intersect(marker_panels$Myeloid, rownames(objs))
pmyl<- FeaturePlot(
  objs,
  features  = feats_myl,
  reduction = "umap.harmony",
  order     = TRUE,
  pt.size   = 0.3,  
  label = TRUE,
  repel=TRUE
)
pmyl
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplot_withmyeloid_markers.jpg",
  units = "in", width = 10, height = 8, dpi = 300
)
#split by condition
pmyl_split <-FeaturePlot(
  objs,
  features  = feats_myl,
  reduction = "umap.harmony",
  split.by = "Condition",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE,
  keep.scale = "feature", 
  cols = c("grey",
           "red")
) 
pmyl_split 
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplotBYCONDITION_withmyeloid_markers.jpg",
  units = "in", width = 12, height = 12, dpi = 300
)
##################################################################################
# flatten to unique gene list: NK cells
feats_nk <- intersect(marker_panels$NK, rownames(objs))
pnk <- FeaturePlot(
  objs,
  features  = feats_nk,
  reduction = "umap.harmony",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE
)
pnk
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplot_withNK_markers.jpg",
  units = "in", width = 10, height = 9, dpi = 300
)
#split by condition
pnk_split <-FeaturePlot(
  objs,
  features  = feats_nk,
  reduction = "umap.harmony",
  split.by = "Condition",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE,
  keep.scale = "feature", 
  cols = c("grey",
           "red")
) 
pnk_split 
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplotBYCONDITION_withNK_markers.jpg",
  units = "in", width = 10, height = 10, dpi = 300
)
##################################################################################
# flatten to unique gene list: T cells
feats_tcell <- intersect(marker_panels$T_cells, rownames(objs))
ptcells <- FeaturePlot(
  objs,
  features  = feats_tcell,
  reduction = "umap.harmony",
  order     = TRUE,
  pt.size   = 0.3,  
  label = TRUE,
  repel=TRUE
)
ptcells
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplot_withTcell_markers.jpg",
  units = "in", width = 10, height = 9, dpi = 300
)
#split by condition
ptcells_split <-FeaturePlot(
  objs,
  features  = feats_tcell,
  reduction = "umap.harmony",
  split.by = "Condition",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE,
  keep.scale = "feature", 
  cols = c("grey",
           "red")
) 
ptcells_split 
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplotBYCONDITION_withTcell_markers.jpg",
  units = "in", width = 12, height = 12, dpi = 300
)
##################################################################################
# flatten to unique gene list: Bcells
feats_bcell <- intersect(marker_panels$B_cells, rownames(objs))
pbcell <- FeaturePlot(
  objs,
  features  = feats_bcell,
  reduction = "umap.harmony",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE
)
pbcell
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplot_withBcells_markers.jpg",
  units = "in", width = 10, height = 9, dpi = 300
)
#split by condition
pbcell_split <-FeaturePlot(
  objs,
  features  = feats_bcell,
  reduction = "umap.harmony",
  split.by = "Condition",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE,
  keep.scale = "feature", 
  cols = c("grey",
           "red")
) 
pbcell_split 
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplotBYCONDITION_withBcells_markers.jpg",
  units = "in", width = 10, height = 10, dpi = 300
)
##################################################################################
# flatten to unique gene list: endothelial
feats_endo <- intersect(marker_panels$Endothelial, rownames(objs))
pendo <- FeaturePlot(
  objs,
  features  = feats_endo,
  reduction = "umap.harmony",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE
)
pendo
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplot_withendothelial_markers.jpg",
  units = "in", width = 10, height = 9, dpi = 300
)
#split by condition
pendo_split <-FeaturePlot(
  objs,
  features  = feats_endo,
  reduction = "umap.harmony",
  split.by = "Condition",
  order     = TRUE,
  pt.size   = 0.3,
  label = TRUE,
  repel=TRUE,
  keep.scale = "feature", 
  cols = c("grey",
           "red")
) 
pendo_split 
# Save image
ggsave(
  filename = "6_harmonyintegrated_featureplotBYCONDITION_withendothelial_markers.jpg",
  units = "in", width = 12, height = 12, dpi = 300
)
#################################################################################