# --------------------------------------------------------------
# Seurat Single-Cell RNA-Seq Analysis for Broccoli Senescence
# Author: Utsab Ghimire
# Date: 03-02-2025
# Description: This script performs preprocessing, normalization,
#              dimensionality reduction, clustering, and marker
#              gene identification using Seurat.
# --------------------------------------------------------------

# ------------------- 1. Install and Load Packages -------------------
# Check and install required packages
required_packages <- c("Matrix", "Seurat")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(Matrix)
library(Seurat)

# ------------------- 2. Set Working Directory -------------------
# Define the data directory (modify based on your system)
data_dir <- file.path("~", "Downloads", "Pre1", "outs", "filtered_feature_bc_matrix")
setwd(data_dir)

# ------------------- 3. Load Data -------------------
# Define file names
matrix_file <- "matrix.mtx"
barcodes_file <- "barcodes.tsv.gz"
features_file <- "features.tsv.gz"

# Check if files exist before loading
if (!all(file.exists(matrix_file, barcodes_file, features_file))) {
  stop("One or more required files are missing in the directory: ", data_dir)
}

# Read the data
expression_matrix <- readMM(matrix_file)                 # Load matrix
barcodes <- readLines(barcodes_file)                     # Load barcodes (cell identifiers)
features <- read.delim(features_file, header = FALSE)    # Load features (gene information)

# Assign row and column names to the matrix
rownames(expression_matrix) <- features[, 2]  # Use gene names from the features file
colnames(expression_matrix) <- barcodes       # Use barcodes as column names

# Check matrix dimensions
dim(expression_matrix)

# ------------------- 4. Create Seurat Object -------------------
seurat_obj <- CreateSeuratObject(counts = expression_matrix, project = "Broccoli_Senescence")

# Check a summary of the Seurat object
print(seurat_obj)

# ------------------- 5. Quality Control & Filtering -------------------
# View QC metrics (number of genes and RNA counts per cell)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Calculate mitochondrial percentage (modify regex for plant mitochondria)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")  # Adjust as needed

# Filter cells based on quality control thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ------------------- 6. Normalize Data -------------------
seurat_obj <- NormalizeData(seurat_obj)

# Access normalized data for verification
head(GetAssayData(seurat_obj, slot = "data")[, 1:5])

# ------------------- 7. Identify Highly Variable Genes -------------------
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Visualize variable features
top10 <- head(VariableFeatures(seurat_obj), 10)
VariableFeaturePlot(seurat_obj) + LabelPoints(points = top10, repel = TRUE)

# ------------------- 8. Scale Data -------------------
seurat_obj <- ScaleData(seurat_obj)

# ------------------- 9. Perform PCA -------------------
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Visualize PCA results
DimPlot(seurat_obj, reduction = "pca")
ElbowPlot(seurat_obj)

# ------------------- 10. Clustering -------------------
# Find neighbors and cluster cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)  # Adjust number of PCs if necessary
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# ------------------- 11. Run UMAP -------------------
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Visualize UMAP with clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# ------------------- 12. Identify Cluster Markers -------------------
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View the top markers for each cluster
head(cluster_markers)

# Save the markers to a CSV file
write.csv(cluster_markers, "cluster_markers.csv", row.names = FALSE)

# ------------------- 13. Save the Seurat Object -------------------
saveRDS(seurat_obj, file = "seurat_obj.rds")

# ------------------- 14. Load Seurat Object (For Future Use) -------------------
# To reload later, use:
# seurat_obj <- readRDS("seurat_obj.rds")

# ------------------- END OF SCRIPT -------------------
