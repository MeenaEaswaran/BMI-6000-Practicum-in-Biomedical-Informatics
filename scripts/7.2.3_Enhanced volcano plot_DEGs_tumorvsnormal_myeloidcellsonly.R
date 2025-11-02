# Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")
#install.packages("readr")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("BiocManager")
#BiocManager::install('EnhancedVolcano')

#import data
library(readr)
myl <- read.csv("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/imagesORcsv/7_Differentialexpressionforspecificcelltypes_Normalvscontrol/RAW_unfiltered_DE_Tumor_vs_Normal_allmyeloidcells_SCT.csv")
head(myl)
library(ggplot2)
library(EnhancedVolcano)

#For the most basic volcano plot, only a single data-frame, data-matrix, or tibble of test results is required, containing point labels, log2FC, and adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
EnhancedVolcano(myl,
                lab = rownames(myl),
                x = 'avg_log2FC',
                y = 'p_val_adj')

#Modify cut-offs for log2FC and P value; specify title; adjust point and label size; remove grid lines; legend positioning
EnhancedVolcano(myl,
                lab = myl$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Tumor vs. normal - Myeloid cells',
                xlab = expression(bold(log[2] ~ fold ~ change)),
                ylab = expression (bold(''-''~log[10] ~ padj)),
                pCutoff = 5.00E-02,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 4.0,
                labCol = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                colAlpha = 1,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0)
#TOsave images without changing up and down colors
ggsave("7_myeloidcells_volcanoplot_defaultcolors_TvsN.jpg", units="in", width=10, height=8, dpi=300)

#Modify to unicolor DEGs all in red
EnhancedVolcano(myl,
                lab = myl$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Tumor vs. normal - Myeloid cells',
                xlab = expression(bold(log[2] ~ fold ~ change)),
                ylab = expression (bold(''-''~log[10] ~ padj)),
                pCutoff = 5.00E-02,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                col=c('black', 'black', 'black', '#DC3220'),
                colAlpha = 1,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0)
#Modify colors
keyvals <- ifelse(
  myl$avg_log2FC < -0.58 & myl$p_val_adj <= 0.05 ,'#0072b2', 
  ifelse(myl$avg_log2FC > 0.58 & myl$p_val_adj <= 0.05 , '#DC3220', '#D6D1D1')
)
keyvals[is.na(keyvals)] <- '#D6D1D1'
names(keyvals)[keyvals == '#DC3220'] <- 'Upregulated'
names(keyvals)[keyvals == '#D6D1D1'] <- 'Not significant'
names(keyvals)[keyvals == '#0072b2'] <- 'Downregulated'

# Calculate counts
upregulated_count <- sum(myl$avg_log2FC > 0.58 & myl$p_val_adj <= 0.05, na.rm = TRUE)
downregulated_count <- sum(myl$avg_log2FC < -0.58 & myl$p_val_adj <= 0.05, na.rm = TRUE)

# Modify the EnhancedVolcano plot
EnhancedVolcano(myl,
                lab = myl$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Tumor vs. normal - Myeloid cells',
                selectLab = rownames(myl)[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))],
                xlab = expression(bold(log[2] ~ fold ~ change)),
                ylab = expression (bold(''-''~log[10] ~ padj)),
                pCutoff = 5.00E-02,
                FCcutoff = 0.58,
                pointSize = 1.0,
                labSize = 4.0,
                #xlim = c(-5, 5),
                #ylim = c(0, 30),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                colCustom = keyvals,
                colAlpha = 4/5) +
  # Annotate downregulated count on the left side top
  annotate("text", x = 1, y = 260, label = paste0("#Downregulated genes: ", downregulated_count),
           hjust = 0, size = 4, color = "#0072b2") +
  # Annotate upregulated count on the right side top
  annotate("text", x = 1, y = 275, label = paste0("#Upregulated genes: ", upregulated_count),
           hjust = 0, size = 4, color = "#DC3220")

ggsave("7_myeloidcells_volcanoplot_finalcolors_TvsN.jpg", units="in", width=10, height=8, dpi=300)

#plot with select genes from feature plots for Tumor vs. normal 
# Define the specific genes you want to label from feature plots (initial) 
#Define subtypes as needed
#additional_genes <- c("CD68") # CD33 = SIGLEC3 #all myeloid cells
#CD33, CD1E, LYZ, LAMP3 are not significantly differentially expressed. So did not show in plot

#additional_genes <- c("ITGAM", "CD66b", "CXCR2") #ITGAM-CD11b, CEACAM8-CD66b #neutrophils
#ITGAM", "CD66b", "CXCR2" are not significantly differentially expressed. So did not show in plot .

#additional_genes <- c("CD274", "CD163")  # CD-274 (PD-L1-found in other cell types ) #macrophages
#CD274 is not significantly differentially expressed.

additional_genes <- c("CD68", "CD163")  # CD-274  #macrophages #final

# Check if the additional genes are in the dataset
valid_additional_genes <- additional_genes[additional_genes %in% myl$gene]
# Combine the additional genes with the existing selected labels
select_labels <- unique(c(rownames(myl)[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))], valid_additional_genes))

#plotnow with annotated hub genes
EnhancedVolcano(myl,
                lab = myl$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Tumor vs. normal - Myeloid cells (macrophages)', #change as neutrophil and macrophage subtypes when required
                selectLab = select_labels,
                xlab = expression(bold(log[2] ~ fold ~ change)),
                ylab = expression (bold(''-''~log[10] ~ padj)),
                pCutoff = 5.00E-02,
                FCcutoff = 0.58,
                pointSize = 1, #1
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'top',
                legendLabSize = 12, #14
                legendIconSize = 4.0, #2
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                colCustom = keyvals,
                colAlpha = 4/5)

#TOsave images without changing up and down colors
ggsave("7_myeloidcells_volcanoplot_finalcolors_TvsN_allcellsannotated(macrophagesmarkersonlyexpressed).jpg", units="in", width=10, height=8, dpi=300)
########################################################################################

