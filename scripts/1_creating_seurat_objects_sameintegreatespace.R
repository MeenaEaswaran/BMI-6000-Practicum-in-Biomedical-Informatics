#Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

#Install Seurat
install.packages('Seurat')

#Load the library
library(Seurat)

#Laryngeal squamous cell carcinoma (LSCC) vs normal tissue 
#Dataset used  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206332

#Downstream analysis:
  #cellular heterogeneity and regulatory landscape 
  #cellular composition, gene regulatory networks, and cell-cell communication 

#KEY TERMS:
#BARCODES (COLUMNS) ARE SHORT DNA BARCODE TAGS TO IDENTIFY READS THAT ARE FROM THE SAME CELL.
#FEATURES ARE GENES (ROWS).
#COUNT MATRIX: A COUNT MATRIX REPRESENTING THE NUMBER OF UNIQUE OBERSERVATIONS OF EACH FEATURE WITHIN EACH CELL BARCODE.
#UMI: UNIQUE MOLECULAR IDENTIFIERS; MOLECULAR TAGS THAT CAN BE APPLIED TO DETECT & QUANTIFY UNIQUE TRANSCRIPTS. 
#DOUBLETS: WHEN TWO CELLS ARE ENCAPSULATED INTO ONE REACTION VOLUME. (IMPORTANT DURING QC)

#Read barcodes, features (genes) and matrix files for each patient's carcinoma in tissue (T) or normal mucosal tissue (N) sample
# All controls
patient_1_N.data <- Read10X(data.dir = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/1N/")
patient_2_N.data <- Read10X(data.dir = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2N/") 
patient_3_N.data <- Read10X(data.dir = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3N/") 
# All cases
patient_1_T.data <- Read10X(data.dir = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/1T/") 
patient_2_T.data <- Read10X(data.dir = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/2T/") 
patient_3_T.data <- Read10X(data.dir = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/3T/") 

#dgCMatrix is comppressed column oriented numeric matrix. This needs to be converted to Seurat objects.
#dots mean no expression (no count);features or gene in rows;  each cell barcode is in the columns;
#The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.

#View data organization for one or all tissue samples IF REQUIRED (SKIP)
# All controls
patient_1_N.data 
patient_2_N.data 
patient_3_N.data 
# All cases
patient_1_T.data 
patient_2_T.data 
patient_3_T.data 

#Load the library
library(tidyverse)

#Create Seurat objects
#All controls
patient_1_N <- CreateSeuratObject(counts = patient_1_N.data, project = "patient_1_N", min.cells = 3, min.features = 200) 
#View class of object if required to confirm that it is a Seurat object
class(patient_1_N) #SKIP IF NOT NEEDED
patient_2_N <- CreateSeuratObject(counts = patient_2_N.data, project = "patient_2_N", min.cells = 3, min.features = 200) 
patient_3_N <- CreateSeuratObject(counts = patient_3_N.data, project = "patient_3_N", min.cells = 3, min.features = 200) 

#All cases
patient_1_T <- CreateSeuratObject(counts = patient_1_T.data, project = "patient_1_T", min.cells = 3, min.features = 200) 
patient_2_T <- CreateSeuratObject(counts = patient_2_T.data, project = "patient_2_T", min.cells = 3, min.features = 200) 
patient_3_T <- CreateSeuratObject(counts = patient_3_T.data, project = "patient_3_T", min.cells = 3, min.features = 200) 

#View Seurat objects IF REQUIRED OR CONFIRM FROM THE ENVIRONMENT ON THE RIGHT SIDE
#All controls
patient_1_N
patient_2_N
patient_3_N
#All cases
patient_1_T
patient_2_T
patient_3_T

#To view column and row names for each sample (SKIP IF NOT DOING THE ABOVE VIEW STEP)
colnames(patient_1_N[]) #will show barcodes
rownames(patient_1_N[]) #will show genes
view(patient_1_N) 
view(patient_1_N@meta.data)

#Save Seurat objects for individual control and case samples
#All controls
saveRDS(patient_1_N, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_1_N.RDS")
saveRDS(patient_2_N, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_2_N.RDS")
saveRDS(patient_3_N, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_3_N.RDS")
#All cases
saveRDS(patient_1_T, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_1_T.RDS")
saveRDS(patient_2_T, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_2_T.RDS")
saveRDS(patient_3_T, file = "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/patient_3_T.RDS")
