# BMI 6000: Practicum in Biomedical Informatics

This repository documents a **bioinformatics project** for the BMI 6000: Practicum in Biomedical Informatics course. The project, titled **"Single-cell analysis of the laryngeal squamous cell carcinoma (LSCC) and healthy adjacent tissue"**, aimed to characterize the cellular heterogeneity and regulatory landscape of LSCC at a single-cell resolution. Using single-cell transcriptomic methods in R and public datasets from NCBI GEO, the study identified alterations in cellular composition, gene regulatory networks, and intracellular communication that distinguish carcinoma in situ LSCC from neighboring healthy epithelial tissue. 

---
# Single-cell analysis of the laryngeal squamous cell carcinoma (LSCC) and healthy adjacent tissue

## Project Overview

#### Bioinformatic Workflow 

![Single-cell RNA-seq workflow](assets/Figure%201_Bioinformatic%20workflow.png)

---

### 1. Data Retrieval

- The LSCC single-cell RNA-sequencing dataset was retrieved from **NCBI GEO** under accession ID **[GSE206332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206332)**.
- **Study Cohort:**
  - The dataset consists of samples collected from **3 individual patients**. For each patient, exactly two tissue types were obtained:
  - **Tumor:** Carcinoma in situ lesion (3 samples in total).
  - **Normal:** Healthy/adjacent non-cancerous normal laryngeal epithelial mucosa (3 samples in total).

---

### 2. Seurat Object Creation, Quality Control, and Data Preprocessing
- **Toolkit used:** Seurat version 5.0 R package for creating a Seurat object, initial quality control performed per sample, and data preprocessing.

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

- R scripts can be found in this [folder](scripts).

**Figures:**
## Normal samples initial QC

![Normal samples Inital QC](assets/Figure%202_Filtering%20_normal%20samples.png)

## Tumor samples initial QC

![Tumor samples initial QC](assets/Figure%203_Filtering%20tumor%20samples.png)

---

### 3. Normalization and Cell-Cycle Assessment
**- **Toolkit used:** the sctransform R package for independent-sample normalization, in place of the standard Seurat datapreprocessing workflow.

-  Mitochondrial, ribosomal, and hemoglobin-expression percentages were evaluated as technical covariates.
- `percent.mt`, `percent.ribo`, and `percent.hb` were regressed during normalization.
-  `G1, S, and G2/M-phase` cell-cycle scores were assessed at the individual-sample level and were also regressed during normalization.

Cell cycle effects had minimal influence on the overall transcriptional heterogeneity in both normal and tumor samples during SCTransform normalization.

- R scripts can be found in this [folder](scripts).

**Figures:**
## Normal samples: Effect of regressing cell cycle scores during SCTransform normalization

![Normal samples: Effect of regressing cell cycle scores during SCTransform normalization](assets/Figure%204_sctransform%20with%20cell%20cycle%20score%20regression_normal.png)

## Tumor samples: Effect of regressing cell cycle scores during SCTransform normalization

![Tumor samples: Effect of regressing cell cycle scores during SCTransform normalization](assets/Figure%205_sctransform%20with%20cell%20cycle%20score%20regression_tumor.png)

---

### 4. Sample-Level Dimensionality Reduction, Clustering, Doublet Detection and Removal

Initial dimensionality reduction and clustering were performed independently for each sample.

- Principal Component Analysis (PCA)
- Nearest-neighbor graph construction
- Initial graph-based clustering
- Uniform Manifold Approximation and Projection (UMAP)
- t-distributed Stochastic Neighbor Embedding (t-SNE)

These analyses were used to evaluate sample quality, identify major cell populations, examine technical effects, and support doublet detection.

- Doublets were identified independently within each sample using the **DoubletFinder** R package.
- Only cells classified as singlets were retained for downstream integration and analysis.
- Doublet removal was performed before integration to reduce the influence of artificial cell profiles on clustering and cell-type annotation.

---

### 5. Data Integration and Batch-Effect Correction

Multiple integration approaches were evaluated to reduce patient- and sample-specific batch effects while preserving biologically meaningful variation.

- **Integration methods compared:**
  - Canonical Correlation Analysis (CCA)
  - Reciprocal PCA (RPCA)
  - Harmony

The resulting embeddings and cluster structures were compared across methods. **Harmony** was selected for downstream analyses based on sample mixing and preservation of biologically interpretable cell populations.

---

### 6. Integrated Dimensionality Reduction and Clustering

Following data integration:

- PCA was used to summarize major sources of variation.
- Graph-based clustering was performed on the integrated dataset.
- UMAP and t-SNE embeddings were generated for visualization.
- Cell-cycle scores were reevaluated on the integrated dataset.
- Cluster structure and sample distribution were examined before cell-type annotation.

---

### 7. Cell-Type Annotation and Proportion Analysis

Cell identities were assigned using a combination of:

- Literature-guided canonical marker genes
- Cluster-specific gene-expression patterns
- Automated reference-based annotation using **SingleR**
- Manual review of epithelial, immune, stromal, and endothelial marker expression

Major annotated cell populations included:

- Epithelial cells
- T cells
- B cells
- Myeloid cells
- Natural killer cells
- Fibroblasts
- Endothelial cells

Cell-type proportions were calculated by sample and experimental condition to compare cellular composition between tumor and normal tissues.

---

### 8. Differential Expression and Functional Enrichment Analysis

Differential expression analysis was performed between **tumor and normal cells within each annotated cell type**.

- **Comparison:** Tumor vs. Normal
- **Differential-expression criteria:**
  - Adjusted *p*-value < 0.05
  - |Average log2 fold change| > 0.58

Differentially expressed genes were evaluated using:

- Gene Ontology
- KEGG pathways
- Metascape functional enrichment
- Cell-type-specific biological interpretation

The analysis identified transcriptional alterations associated with epithelial remodeling, immune activation, extracellular-matrix organization, tumor-associated signaling, and changes in tissue homeostasis.

---

### 9. Gene Regulatory Network Analysis

Gene regulatory network analysis was performed using **SCENIC**.

SCENIC was used to:

- Identify transcription factor–target gene modules
- Infer regulon activity
- Compare regulatory programs between tumor and normal conditions
- Examine cell-type-specific transcriptional regulation

Complementary transcription factor enrichment was performed using **Metascape** and the **TRRUST** database.

Regulatory programs involving transcription factors such as **TP63, TP53, HES1, KLF, ETS, STAT, FOS, and JUN family members** were examined in the context of epithelial transformation and tumor-associated signaling.

---

### 10. Cell-Cell Communication Analysis

Intercellular signaling networks were evaluated using **CellChat**.

Tumor and normal datasets were compared based on:

- Number of inferred interactions
- Overall interaction strength
- Sender and receiver cell populations
- Signaling pathway activity
- Ligand-receptor interactions

This analysis was used to investigate how tumor-associated epithelial, immune, endothelial, and stromal populations may reshape the local tissue microenvironment.

---

For questions or issues, please contact the repository maintainer. Refer to the [presentation slides](assets/Easwaran_Meena_Practicum_RESULTS%20SUMMARY_COMPILED-compressed.pdf) and/or the [final practicum report](assets/Final_Practicum_report_Fall%202025_Easwaran_Meena.pdf)for detailed information and results.

This repository is **solely for educational purposes and serves as a backup** for my graduate school capstone paper and Master's degree requirements related to the **BMI 6000: Practicum in Biomedical Informatics** course at McWilliams School of Biomedical Informatics at UTHealth Houston.
