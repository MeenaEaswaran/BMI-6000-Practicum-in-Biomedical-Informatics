#Set the working directory
setwd("~/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace")

#Run for the first time
#BiocManager::install("glmGamPoi") #needed for sctransform faster run

#To fix errors
#install.packages("devtools")
#devtools::install_github("satijalab/seurat", ref = "fix/v.5.3.1")

#install UMAP if not already installed
#reticulate::py_install(packages='umap-learn')  

# load the libraries if required
library(Seurat)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
options(future.globals.maxSize = 1e9)
#library(DoubletFinder)

# Read Seurat Objects from previous storage (use this if starting afresh-saves time to not run the previous codes)
patient_1_N<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_1_N_singletsonly.RDS")
patient_2_N<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_2_N_singletsonly.RDS")
patient_3_N<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_3_N_singletsonly.RDS")
patient_1_T<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_1_T_singletsonly.RDS")
patient_2_T<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_2_T_singletsonly.RDS")
patient_3_T<- readRDS("/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/4_patient_3_T_singletsonly.RDS")

#Check default assay
DefaultAssay(object = patient_1_N) 
DefaultAssay(object = patient_2_N) 
DefaultAssay(object = patient_3_N) 
DefaultAssay(object = patient_1_T) 
DefaultAssay(object = patient_2_T) 
DefaultAssay(object = patient_3_T) 

# Add prefixes first so cell names are unique
patient_1_N <- RenameCells(patient_1_N, add.cell.id = "p1N")
patient_2_N <- RenameCells(patient_2_N, add.cell.id = "p2N")
patient_3_N <- RenameCells(patient_3_N, add.cell.id = "p3N")
patient_1_T <- RenameCells(patient_1_T, add.cell.id = "p1T")
patient_2_T <- RenameCells(patient_2_T, add.cell.id = "p2T")
patient_3_T <- RenameCells(patient_3_T, add.cell.id = "p3T")

# Merge all into a single Seurat object for integration
objs <- merge(
  x = patient_1_N,
  y = list(patient_2_N, patient_3_N,
           patient_1_T, patient_2_T, patient_3_T),
  project = "AllSamples"
)

objs

rm(patient_1_N, patient_1_T, patient_2_N, patient_2_T, patient_3_N, patient_3_T)

#save RDs for next steps
saveRDS(objs, "/Users/meena/Desktop/Practicum Fall 2025/Approach_2_case_control_integratedspace/5_all_singlet_objects_combined.rds")