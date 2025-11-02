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
library(scRNAseq)
library(SingleCellExperiment)
library(pheatmap)
options(future.globals.maxSize = 1e9)
#library(DoubletFinder)

#Load integrated object 
objs <- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/6_harmony_integratedwithfinalcelltypeannotations_seurat_30dims.rds")

DefaultAssay(objs) # if needed: DefaultAssay(objs) <- "SCT"

# =========================
# COUNTS & PROPORTIONS with Sample / Condition / Patient
# =========================

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# --- Safety: create NA-removed copy for cell-type summaries if missing ---
if (!exists("plot_obj")) {
  plot_obj <- subset(objs, subset = !is.na(SingleR_ENCODE))
  plot_obj$SingleR_ENCODE <- droplevels(plot_obj$SingleR_ENCODE)
}

# --- 0) Quick checks (all should exist) ---
stopifnot(all(c("SingleR_ENCODE","Condition","Patient","Sample") %in% colnames(plot_obj@meta.data)))
stopifnot(all(c("Condition","Patient","Sample","seurat_clusters") %in% colnames(objs@meta.data)))

# Convert metadata to data.frames
meta_lbl <- as.data.frame(plot_obj@meta.data)  # NA-removed for SingleR_ENCODE
meta_all <- as.data.frame(objs@meta.data)      # all cells for clusters/conditions

## (Optional) if plyr got attached, detach it to avoid conflicts
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

## NA-removed metadata for cell-type summaries
meta_lbl <- tibble::as_tibble(plot_obj@meta.data, rownames = "cell_id")

## All-cell metadata for cluster/condition totals
meta_all <- tibble::as_tibble(objs@meta.data, rownames = "cell_id")

# ---------- 1) Cell-type summaries (NA-removed) ----------
#celltype_overall <- meta_lbl %>%
#dplyr::count(SingleR_ENCODE, name = "n") %>%
#dplyr::arrange(dplyr::desc(n)) %>%
#dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1))

celltype_by_condition <- meta_lbl %>%
  dplyr::filter(!is.na(Condition)) %>%
  dplyr::count(Condition, SingleR_ENCODE, name = "n") %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1)) %>%
  dplyr::ungroup()

#celltype_by_patient <- meta_lbl %>%
# dplyr::filter(!is.na(Patient)) %>%
#dplyr::count(Patient, SingleR_ENCODE, name = "n") %>%
#dplyr::group_by(Patient) %>%
#dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1)) %>%
#dplyr::ungroup()

celltype_by_sample <- meta_lbl %>%
  dplyr::filter(!is.na(Sample)) %>%
  dplyr::count(Sample, SingleR_ENCODE, name = "n") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1)) %>%
  dplyr::ungroup()

# ---------- 2) Cluster summaries (all cells) ----------
cluster_factor <- function(x) {
  xs <- as.character(x)
  if (all(grepl("^\\d+$", xs))) factor(xs, levels = sort(unique(as.numeric(xs))), ordered = TRUE)
  else forcats::fct_inorder(xs)
}

#cluster_counts <- meta_all %>%
#dplyr::mutate(seurat_clusters = cluster_factor(seurat_clusters)) %>%
#dplyr::count(seurat_clusters, name = "n") %>%
# dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1))

cluster_by_condition <- meta_all %>%
  dplyr::filter(!is.na(Condition)) %>%
  dplyr::mutate(seurat_clusters = cluster_factor(seurat_clusters)) %>%
  dplyr::count(Condition, seurat_clusters, name = "n") %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1)) %>%
  dplyr::ungroup()

# ---------- 3) Condition / Patient / Sample totals (all cells) ----------
condition_counts <- meta_all %>%
  dplyr::filter(!is.na(Condition)) %>%
  dplyr::count(Condition, name = "n") %>%
  dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1))

#patient_counts <- meta_all %>%
#dplyr::filter(!is.na(Patient)) %>%
# dplyr::count(Patient, name = "n") %>%
#dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1))

sample_counts <- meta_all %>%
  dplyr::filter(!is.na(Sample)) %>%
  dplyr::count(Sample, name = "n") %>%
  dplyr::mutate(prop = n/sum(n), prop_pct = round(100*prop, 1))


library(ggplot2)
library(forcats)
library(scales)

# ==== CELL TYPES (from plot_obj; NA-removed) ====

## 1) By Condition
# 100% stacked (proportions)

p_ct_cond_prop <- ggplot(celltype_by_condition,
                         aes(x = Condition, y = prop, fill = SingleR_ENCODE)) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.01),   # 0.00–1.00
    breaks = seq(0, 1, 0.25),                 # 0.00, 0.25, 0.50, 0.75, 1.00
    limits = c(0, 1)
  ) +
  labs(title = "Cell-type composition by Condition",
       x = NULL, y = "Proportion") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL))

p_ct_cond_prop

ggsave("6_celltypeproprotion_bycondition.jpg", p_ct_cond_prop, width = 8, height = 5, dpi = 300)

# Stacked counts (IGNORE)
p_ct_cond_counts <- ggplot(celltype_by_condition,
                           aes(x = Condition, y = n, fill = SingleR_ENCODE)) +
  geom_col(width = 0.8) +
  labs(title = "Cell-type counts by Condition",
       x = NULL, y = "Cells") +
  theme_classic()
p_ct_cond_counts 
ggsave("celltype_by_condition_counts.jpg", p_ct_cond_counts, width = 8, height = 5, dpi = 300)


## 2) By Patient (IGNORE)
# Order patients as they appear
celltype_by_patient$Patient <- fct_inorder(celltype_by_patient$Patient)

p_ct_patient_prop <- ggplot(celltype_by_patient,
                            aes(x = Patient, y = prop, fill = SingleR_ENCODE)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  labs(title = "Cell-type composition by Patient (100%)",
       x = "Patient", y = "Proportion") +
  theme_classic()
p_ct_patient_prop

p_ct_patient_counts <- ggplot(celltype_by_patient,
                              aes(x = Patient, y = n, fill = SingleR_ENCODE)) +
  geom_col() +
  labs(title = "Cell-type counts by Patient",
       x = "Patient", y = "Cells") +
  theme_classic()
p_ct_patient_counts
ggsave("celltype_by_patient_100pct.jpg", p_ct_patient_prop, width = 10, height = 5, dpi = 300)
ggsave("celltype_by_patient_counts.jpg", p_ct_patient_counts, width = 10, height = 5, dpi = 300)

## 3) By Sample 
celltype_by_sample$Sample <- fct_inorder(celltype_by_sample$Sample)

p_ct_sample_prop <- ggplot(celltype_by_sample,
                           aes(x = Sample, y = prop, fill = SingleR_ENCODE)) +
  geom_col() +
  scale_y_continuous(
    labels = label_number(accuracy = 0.01),   # 0.00–1.00
    breaks = seq(0, 1, 0.25),                 # 0.00, 0.25, 0.50, 0.75, 1.00
    limits = c(0, 1)
  )  +
  labs(title = "Cell-type composition by Patient Sample",
       x = "Sample", y = "Proportion") +
  #coord_flip() +
  theme_classic() + guides(fill = guide_legend(title = NULL))

p_ct_sample_prop 
ggsave("6_celltypeproprotion_bypatientsample.jpg", p_ct_sample_prop, width = 8, height = 5, dpi = 300)

p_ct_cond_prop <- ggplot(celltype_by_condition,
                         aes(x = Condition, y = prop, fill = SingleR_ENCODE)) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.01),   # 0.00–1.00
    breaks = seq(0, 1, 0.25),                 # 0.00, 0.25, 0.50, 0.75, 1.00
    limits = c(0, 1)
  ) +
  labs(title = "Cell-type composition by Condition",
       x = NULL, y = "Proportion") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL))

#IGNORE BELOW
p_ct_sample_counts <- ggplot(celltype_by_sample,
                             aes(x = Sample, y = n, fill = SingleR_ENCODE)) +
  geom_col() +
  labs(title = "Cell-type counts by Sample",
       x = "Sample", y = "Cells") +
  coord_flip() +
  theme_classic()
p_ct_sample_counts
