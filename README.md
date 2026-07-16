# BMI 6000: Practicum in Biomedical Informatics

This repository documents a **bioinformatics project** for the BMI 6000: Practicum in Biomedical Informatics course. The project, titled **"Single-cell analysis of the laryngeal squamous cell carcinoma (LSCC) and healthy adjacent tissue"**, aimed to characterize the cellular heterogeneity and regulatory landscape of LSCC at a single-cell resolution. Using single-cell transcriptomic methods in R and public datasets from NCBI GEO, the study identified alterations in cellular composition, gene regulatory networks, and intracellular communication that distinguish carcinoma in situ LSCC from neighboring healthy epithelial tissue. 

---
# Single-cell analysis of the laryngeal squamous cell carcinoma (LSCC) and healthy adjacent tissue

## Project Overview

#### Bioinformatic Workflow: 

The analysis was performed primarily in **R** using **Seurat version 5.0**, with additional tools for data normalization, doublet detection, automated cell-type annotation, functional enrichment, gene regulatory network analysis, and cell-cell communication analysis.

<img src="assets/Figure%201_Bioinformatic%20workflow.png" alt="Single-cell RNA-seq workflow" width="600">

---

### 1. Data Retrieval

- The LSCC single-cell RNA-sequencing dataset was retrieved from **NCBI GEO** under accession ID **[GSE206332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206332)**.
- **Study Cohort:**
  - The dataset consists of samples collected from **3 individual patients**. For each patient, exactly two tissue types were obtained:
  - **Tumor:** Carcinoma in situ lesion (3 samples in total).
  - **Normal:** Healthy/adjacent non-cancerous normal laryngeal epithelial mucosa (3 samples in total).

---

### 2. Seurat Object Creation, Quality Control, and Data Preprocessing

- **Quality-control metrics:**
  - `nFeature_RNA`
  - `nCount_RNA`
  - `percent.mt`
  - `percent.ribo`
  - `percent.hb`

- **Filtering criteria:**
  - `nFeature_RNA`: 200–6000
  - `nCount_RNA`: <20000
  - `percent.mt`: <20%
  - `percent.hb`: <1%

Cells with low transcript complexity, unusually high transcript counts, elevated mitochondrial RNA, or excessive hemoglobin-gene expression were removed.

- R scripts with "1_ and 2_" as prefixes can be found in this [folder](scripts).

## Normal samples initial QC

<img src="assets/Figure%202_Filtering%20_normal%20samples.png" alt="Normal samples Inital QC" width="600">

## Tumor samples initial QC

<img src="assets/Figure%203_Filtering%20tumor%20samples.png" alt="Tumor samples initial QC" width="600">

---

### 3. Normalization and Cell-Cycle Assessment

- Each sample was normalized independently using **sctransform**, in place of the standard Seurat data normalization workflow.
-  Mitochondrial, ribosomal, and hemoglobin-expression percentages were evaluated as technical covariates.
- `percent.mt`, `percent.ribo`, and `percent.hb` were regressed during normalization.
-  `G1, S, and G2/M-phase` cell-cycle scores were assessed at the individual-sample level and were also regressed during normalization.

Cell cycle effects had minimal influence on the overall transcriptional heterogeneity in both normal and tumor samples during SCTransform normalization.

- R scripts with "3_" as prefix can be found in this [folder](scripts).

## Normal samples: Effect of regressing cell cycle scores during SCTransform normalization

<img src="assets/Figure%204_sctransform%20with%20cell%20cycle%20score%20regression_normal.png" alt="Normal samples: Effect of regressing cell cycle scores during SCTransform normalization" width="600">

## Tumor samples: Effect of regressing cell cycle scores during SCTransform normalization

<img src="assets/Figure%205_sctransform%20with%20cell%20cycle%20score%20regression_tumor.png" alt="Tumor samples: Effect of regressing cell cycle scores during SCTransform normalization" width="600">

---

### 4. Sample-Level Dimensionality Reduction, Clustering, Doublet Detection and Removal

- The analyses below were used to evaluate sample quality, identify major cell populations, examine technical effects, and support doublet detection.
  - Principal Component Analysis (PCA) 
  - Uniform Manifold Approximation and Projection (UMAP)
  - t-distributed Stochastic Neighbor Embedding (t-SNE)

- Doublets were identified independently within each sample using the **DoubletFinder** R package.
- Only cells classified as singlets were retained for downstream integration and analysis.
- Doublet removal was performed before integration to reduce the influence of artificial cell profiles on clustering and cell-type annotation.

- R scripts with "4.1-4.6_" as prefixes can be found in this [folder](scripts).

##  Elbow plots representing the variance across the first 50 principal components across the normal and tumor samples after normalization

<img src="assets/Figure%206_%20Elbowplots_normalandtumor.png" alt="Elbow plots representing the variance across the first 50 principal components across the normal and tumor samples after normalization" width="600">

##  Initial UMAP projections of individual normal and tumor samples before doublet removal and data integration

<img src="assets/Figure%207_%20initial_LDA_UMAPs_normal%20and%20tumor.png" alt="Initial UMAP projections of individual normal and tumor samples before doublet removal and data integration" width="600">

##  Initial t-SNE projections of individual normal and tumor samples before doublet removal and data integration

<img src="assets/Figure%208_%20initial_LDA_t-SNEs_normal%20and%20tumor.png" alt="Initial t-SNE projections of individual normal and tumor samples before doublet removal and data integration" width="600">

##  Identification and removal of doublets from the normal samples

<img src="assets/Figure%209_%20Doubletfinder_normal.png" alt="Identification and removal of doublets from the normal samples" width="600">

##  Identification and removal of doublets from the tumor samples

<img src="assets/Figure%2010_%20Doubletfinder_tumor.png" alt="Identification and removal of doublets from the tumor samples" width="600">

---

### 5. Data Integration and Batch-Effect Correction

- Multiple integration approaches were evaluated to reduce patient- and sample-specific batch effects while preserving biologically meaningful variation.
- **Integration methods compared:**
  - Canonical Correlation Analysis (CCA)
  - Reciprocal PCA (RPCA)
  - Harmony

- CCA, Harmony, and RPCA integration greatly improved shared cell alignment, yielding highly consistent clustering with 23, 24, and 24 clusters, respectively.
- **Harmony-integrated** data was selected for downstream analyses based on its proven accuracy and scalability in aligning diverse single-cell datasets.

- Cell cycle scores showed no major source of variation in the integrated dataset.

- R scripts with "5.1-5.6_" as prefixes can be found in this [folder](scripts).

##  Comparison of data integration approaches across samples

<img src="assets/Figure%2011_%20Integration%20analysis.png" alt="Comparison of data integration approaches across samples" width="600">

##  UMAP visualizations of integrated datasets colored by cell cycle phases

<img src="assets/Figure%2012_%20Integration%20analysis_cell%20cycle%20score.png" alt="UMAP visualizations of integrated datasets colored by cell cycle phases" width="600">

---

### 6. Cell-Type Annotation and Proportion Analysis

- Cell identities were assigned using a combination of:
  - Literature-guided canonical marker genes
  - Cluster-specific gene-expression patterns
  - Automated reference-based annotation using the **SingleR** R package
    - Blueprint/ENCODE, HPCA, DICE, and MonaccoImmuneData reference datasets
  - Manual review of epithelial, immune, stromal, and endothelial marker expression

- Major annotated cell populations included:
  - Epithelial cells
  - Fibroblasts
  - Myeloid cells
  - NK cells
  - B-cells
  - T-cells
  - Endothelial cells

- Cell-type proportions were calculated by sample and experimental condition to compare cellular composition between tumor and normal tissues.

- R scripts with "6.1-6.5_" as prefixes can be found in this [folder](scripts).

##  Literature-guided marker gene expression across Harmony-integrated clusters and experimental conditions

<img src="assets/Figure%2014_%20dotplot_originalstudyandliteraturemarkers_cluster_NandT.png" alt="Literature-guided marker gene expression across Harmony-integrated clusters and experimental conditions" width="600">

##  Visualization of manually annotated cell clusters in the Harmony-integrated dataset

<img src="assets/Figure%2015_%20Manual%20cell%20type%20annotation_UMAP%20and%20t-SNE.png" alt="Visualization of manually annotated cell clusters in the Harmony-integrated dataset" width="600">

##  SingleR automated cell type annotations using multiple reference datasets 

<img src="assets/Figure%2016_SingleRAutomationonly.png" alt="SingleR automated cell type annotations using multiple reference datasets" width="600">

##  Manual + SingleR-automated cell type annotation and composition analysis for the Harmony-integrated dataset

<img src="assets/Figure%2017_finalcellannotation_cellproportion.png" alt="Manual + SingleR-automated cell type annotation and composition analysis for the Harmony-integrated dataset" width="600">

---

### 7. Differential Expression and Functional Enrichment Analysis

- Differential expression analysis was performed between **tumor and normal cells within each annotated cell type**.
  - **Comparison:** Tumor vs. Normal
  - **Differential-expression criteria:**
    - Adjusted *p*-value < 0.05
    - |Average log2 fold change| > 0.58

-Differentially expressed genes were evaluated using **[Metascape](https://metascape.org/gp/index.html#/main/step1)** functional enrichment tool.
  - GO Biological Process (BP), Cellular Compartment (CC), and Molecular Function (MF), and KEGG pathway enrichments.

- Cell-type-specific biological interpretation

- The analysis identified transcriptional alterations associated with extracellular-matrix  remodeling, immune activation, and changes in stress responses.

- R scripts with "7.1_-7.2.7_" as prefixes can be found in this [folder](scripts).

##  Differential expression and functional enrichment analysis of epithelial cells 

<img src="assets/Figure%2018_Differentialexpression_functionalenrichment_epithelial%20cells.png" alt="Differential expression and functional enrichment analysis of epithelial cells" width="600">

##  Differential expression and functional enrichment analysis of fibroblasts

<img src="assets/Figure%2019_Differentialexpression_functionalenrichment_fibroblasts.png" alt="Differential expression and functional enrichment analysis of fibroblasts" width="600">

##  Differential expression and functional enrichment analysis of myeloid cells 

<img src="assets/Figure%2020_Differentialexpression_functionalenrichment_myeloidcells.png" alt="Differential expression and functional enrichment analysis of myeloid cells" width="600">

##  Differential expression and functional enrichment analysis of NK cells 

<img src="assets/Figure%2021_Differentialexpression_functionalenrichment_NKcells.png" alt="Differential expression and functional enrichment analysis of NK cells" width="600">

##  Differential expression and functional enrichment analysis of B-cells 

<img src="assets/Figure%2022_Differentialexpression_functionalenrichment_Bcells.png" alt="Differential expression and functional enrichment analysis of B-cells" width="600">

##  Differential expression and functional enrichment analysis of T-cells 

<img src="assets/Figure%2023_Differentialexpression_functionalenrichment_Tcells.png" alt="Differential expression and functional enrichment analysis of T-cells" width="600">

##  Differential expression and functional enrichment analysis of endothelial cells 

<img src="assets/Figure%2024_Differentialexpression_functionalenrichment_Endothelialcells.png" alt="Differential expression and functional enrichment analysis of endothelial cells" width="600">

---

### 8. Gene Regulatory Network Analysis

- Gene regulatory network analysis was performed using **Single-cell regulatory network inference and clustering (SCENIC; R version)**.

- SCENIC was used to:
  - Infer regulon activity between tumor and normal conditions
  - Examine cell-type-specific transcriptional regulation

- Distinct transcriptional regulatory patterns between tumor and normal conditions across each cell type were observed.

- R scripts with "10.1_-10.3_" as prefixes can be found in this [folder](scripts).

##  SCENIC regulon activity heatmap

<img src="assets/Figure%2025_SCENIC_heatmap.png" alt="SCENIC regulon activity heatmap" width="600">

---

### 9. Cell-Cell Communication Analysis

- Intercellular signaling networks were evaluated using **CellChat** R package, to identify ligand-receptor interactions across the tumor and normal datasets.
  
- This analysis was used to investigate how tumor-associated epithelial, immune, endothelial, and stromal populations may reshape the local tissue microenvironment.

- R scripts with "9.1_-9.3_" as prefixes can be found in this [folder](scripts).

##  Cell-cell communication analysis of the normal and tumor samples

<img src="assets/Figure%2027_Intercellularcommunication_cellchat.png" alt="Cell-cell communication analysis of the normal and tumor samples" width="600">

##  Epithelial cells: Significant L-R interactions

<img src="assets/Figure%2028_Sig_L-R_epithelialcells.png" alt="Epithelial cells: Significant L-R interactions" width="600">

##  Fibroblasts: Significant L-R interactions

<img src="assets/Figure%2029_Sig_L-R_fibroblasts.png" alt="Fibroblasts: Significant L-R interactions" width="600">

##  Myeloid cells: Significant L-R interactions

<img src="assets/Figure%2030_Sig_L-R_myeloidcells.png" alt="Myeloid cells: Significant L-R interactions" width="600">

##  NK cells: Significant L-R interactions 

<img src="assets/Figure%2031_Sig_L-R_NKcells.png" alt="NK cells: Significant L-R interactions " width="600">

##  B-cells: Significant L-R interactions 

<img src="assets/Figure%2032_Sig_L-R_Bcells.png" alt="B-cells: Significant L-R interactions" width="600">

##  T-cells: Significant L-R interactions 

<img src="assets/Figure%2033_Sig_L-R_Tcells.png" alt="T-cells: Significant L-R interactions " width="600">

##  Endothelial cells: Significant L-R interactions 

<img src="assets/Figure%2034_Sig_L-R_Endothelialcells.png" alt="Endothelial cells: Significant L-R interactions " width="600">

---

For questions or issues, please contact the repository maintainer. Refer to the [presentation slides](assets/Easwaran_Meena_Practicum_RESULTS%20SUMMARY_COMPILED-compressed.pdf) and/or the [final practicum report](assets/Final_Practicum_report_Fall%202025_Easwaran_Meena.pdf)for detailed information and results.

This repository is **solely for educational purposes and serves as a backup** for my graduate school capstone paper and Master's degree requirements related to the **BMI 6000: Practicum in Biomedical Informatics** course at McWilliams School of Biomedical Informatics at UTHealth Houston.
