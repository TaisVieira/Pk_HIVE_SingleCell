library(Seurat)
library(dplyr)
library(plotly)
library(clustree)
library(openxlsx)
library(Matrix)
library(patchwork)

setwd('./Analysis/R analysis/Scripts')
source('./filterScrublet.R')
source('./util_funcs.R')
setwd('./Analysis/P_knowlesi/Pk2/')

# =============================================================================
# LOAD DATA AND CREATE INIITAL OBJECT
# =============================================================================
# Load original data
barcodes <- read.delim("merged/matrix/Pk2_barcodes.tsv", header = FALSE)
features <- read.delim("merged/matrix/Pk2_features.tsv", header = FALSE)
matrix <- readMM("merged/matrix/Pk2_matrix.mtx")

colnames(matrix) <- barcodes[,1]
rownames(matrix) <- features[,2]

# Filter for cells with more than 1000 reads
matrix_filtered <- matrix[,colSums(matrix) > 1000]

# Create object with genes present in at least 10 cells and
# cells that have at least 100 genes
so.pk <- CreateSeuratObject(matrix_filtered, min.cells = 10, min.features = 100)
so.pk

so.pk <- NormalizeData(so.pk)
so.pk <- FindVariableFeatures(so.pk)
so.pk <- ScaleData(so.pk)
so.pk <- RunPCA(so.pk, verbose = FALSE)
ElbowPlot(so.pk, ndims = 50)
so.pk <- RunUMAP(so.pk, dims=1:20, n.components = 3, n.neighbors = 25, 
                 min.dist = 0.4, spread = 2.5, repulsion.strength = 2,
                 metric = 'correlation', verbose = FALSE, return.model = TRUE)
so.pk <- FindNeighbors(so.pk)
so.pk <- FindClusters(so.pk, resolution = 0.3)
DimPlot(so.pk, reduction = 'pca')
DimPlot(so.pk, reduction = 'umap')

# =============================================================================
# EXPORT DATA FOR SCRUBLET
# =============================================================================

# Create output directory
output_dir <- "merged/matrix/scrublet_input"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Extract filtered counts matrix (genes x cells)
counts_matrix <- so.pk@assays$RNA$counts

# Export matrix in Matrix Market format
writeMM(counts_matrix, file = file.path(output_dir, "matrix.mtx"))

# Export barcodes (cell names)
barcodes <- colnames(counts_matrix)
write.table(barcodes, 
            file = file.path(output_dir, "barcodes.tsv"),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# Export features (gene names)
features <- rownames(counts_matrix)
write.table(features, 
            file = file.path(output_dir, "features.tsv"),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# =============================================================================
# CREATE FINAL SCRUBLET OBJECT
# =============================================================================

# Run Scrublet in Python

# Apply Scrublet filtering in R
scrublet_results <- apply_scrublet_filtering(
  matrix_filtered = counts_matrix,
  barcodes_filtered = barcodes,
  features = features,
  scrublet_csv_path = "merged/matrix/scrublet_doublet_scores.csv",
  output_dir = "merged/matrix/matrix_filtered_scrublet"
)

# Create object
matrix_scrublet <- scrublet_results$matrix_final
barcodes_scrublet <- scrublet_results$barcodes_final
features_scrublet <- scrublet_results$features

colnames(matrix_scrublet) <- barcodes_scrublet
rownames(matrix_scrublet) <- features_scrublet

so.scrublet <- CreateSeuratObject(matrix_scrublet)
so.scrublet <- NormalizeData(so.scrublet)
so.scrublet <- FindVariableFeatures(so.scrublet)
so.scrublet <- ScaleData(so.scrublet)
so.scrublet <- RunPCA(so.scrublet, verbose = FALSE)
so.scrublet <- RunUMAP(so.scrublet, dims=1:20, n.components = 3, n.neighbors = 25, 
                       min.dist = 0.4, spread = 2.5, repulsion.strength = 2,
                       metric = 'correlation', verbose = FALSE, return.model = TRUE)
so.scrublet <- FindNeighbors(so.scrublet)
so.scrublet <- FindClusters(so.scrublet, resolution = 0.3)
DimPlot(so.scrublet, reduction = 'pca')
DimPlot(so.scrublet, reduction = 'umap')

# =============================================================================
# CHECK FOR OUTLIERS
# =============================================================================

# Visualize gene and read distribution
Idents(so.scrublet) <- 'orig_idents'
VlnPlot(so.scrublet, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
Idents(so.scrublet) <- 'seurat_clusters'

# Find outlier cells
qcount <- quantile(so.scrublet$nCount_RNA, probs = 0.999)
outlier_cells <- WhichCells(so.scrublet, 
                            expression = nCount_RNA > qcount)
so.scrublet$is_outlier <- colnames(so.scrublet) %in% outlier_cells

# Visualize to confirm they're not driving your structure
DimPlot(so.scrublet, group.by = "is_outlier", split.by = 'is_outlier')

# Filter out outliers
cells_keep <- WhichCells(so.scrublet, 
                         expression = nCount_RNA < qcount)
so.pk.clean <- subset(so.scrublet, cells = cells_keep)

# Visualize filtered distribution
hist(so.pk.clean$nFeature_RNA, breaks = 100, main = "nFeature Distribution")
hist(so.pk.clean$nCount_RNA, breaks = 100, main = "nCount Distribution")

summary(so.pk.clean$nCount_RNA)
summary(so.pk.clean$nFeature_RNA)

# =============================================================================
# CREATE FINAL OBJECT
# =============================================================================

# Normalization is preserved from so.scrublet — recompute variable features,
# scaling, and embeddings on the filtered cell set
so.final <- ScaleData(so.pk.clean)
so.final <- RunPCA(so.final, verbose = FALSE)
so.final <- RunUMAP(so.final, dims=1:20, n.components = 2, n.neighbors = 25, 
                    min.dist = 0.4, spread = 2.5, repulsion.strength = 2,
                    metric = 'correlation', verbose = FALSE, return.model = TRUE)
so.final <- FindNeighbors(so.final, dims = 1:20, reduction = "pca")
so.final <- FindClusters(so.final, resolution = 0.4)
DimPlot(so.final, reduction = 'pca')
DimPlot(so.final, reduction = 'umap')

# =============================================================================
# PLOT COMPARISON
# =============================================================================

wrap_plots(
  DimPlot(so.pk, 
          reduction = "umap") +
    NoAxes() +
    ggtitle("No Filter"),
  
  DimPlot(so.scrublet, 
          reduction = "umap") +
    NoAxes() +
    ggtitle("After Scrublet"),
  
  DimPlot(so.final, 
          reduction = "umap") +
    NoAxes() +
    ggtitle("Without Outliers"),
  
  FeaturePlot(so.pk, 
              features = 'nCount_RNA') +
    NoAxes() +
    ggtitle("No Filter Count"),
  
  FeaturePlot(so.scrublet, 
              features = 'nCount_RNA') +
    NoAxes() +
    ggtitle("After Scrublet Count"),
  
  FeaturePlot(so.final, 
              features = 'nCount_RNA') +
    NoAxes() +
    ggtitle("Without Outliers Count"),
  
  ncol = 3
) + plot_layout(guides = "collect")

VlnPlot(so.final, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, 
        pt.size = 0, group.by = 'orig.ident')

# =============================================================================
# FIND DEG
# =============================================================================

# Find cluster specific markers
deg_clusters <- FindAllMarkers(so.final,  test.use = "wilcox", min.pct = 0.25,
                               min.diff.pct = 0.1, only.pos = TRUE)
deg_clusters_filtered <- deg_clusters[deg_clusters$p_val_adj < 0.05,]
all_markers <- add_gene_descriptions_from_gff_fast(deg_clusters_filtered, 
                                                   "./Analysis/Genomes/Pk-68/PlasmoDB-68_PknowlesiH.gff")
deg_clusters_filtered$gene <- gsub("-", "_", deg_clusters_filtered$gene)
write.xlsx(all_markers, file = "cluster_markers.xlsx", rowNames = FALSE)


# =============================================================================
# SAVE OBJECT
# =============================================================================

saveRDS(so.final, 
        "./Analysis/R analysis/RDS files/S.O.Pk2.rds")




