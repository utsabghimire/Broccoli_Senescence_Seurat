# singlecellbroccoli
Seurat Single-Cell RNA-Seq Analysis for Broccoli Senescence 🌱🔬

Overview

This repository contains an R script for single-cell RNA sequencing (scRNA-seq) analysis using Seurat. The script processes and analyzes broccoli inflorescence tissue to study postharvest senescence at the single-cell level.

Files in This Repository

seurat_analysis.R - Main script for scRNA-seq preprocessing, normalization, clustering, and marker identification.
cluster_markers.csv - (Optional) Output file containing marker genes for identified clusters.
Setup and Requirements

1. Install Required R Packages
Ensure you have R (>=4.0) installed. Install necessary libraries in R:

install.packages(c("Matrix", "Seurat"))
Then, load them:

library(Matrix)
library(Seurat)
2. Input Data Requirements
The script expects:

A filtered feature barcode matrix from 10X Genomics (matrix.mtx, barcodes.tsv.gz, features.tsv.gz).
The working directory should be set correctly:
setwd("path_to_your_data_directory")
How to Use This Script

Clone this repository:
git clone https://github.com/your-username/your-repo.git
cd your-repo
Open RStudio and run seurat_analysis.R.
The script performs:
Data preprocessing & QC
Normalization & feature selection
PCA, clustering, and UMAP visualization
Identification of cluster markers
Outputs

Seurat Object: Processed scRNA-seq data (seurat_obj.rds).
Cluster Marker Genes: A .csv file with marker genes per cluster.
Future Enhancements

Integrate with Monocle3 for trajectory inference
Add automated cell-type annotation

Author

Utsab Ghimire
