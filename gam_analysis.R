library(Seurat)
library(patchwork)
library(dplyr)
library(openxlsx)
library(purrr)

setwd('./Analysis/R analysis/Scripts')
source('./util_funcs.R')
source('./pseudotime_funcs.R')
setwd('./Analysis/Paper/pk')

# =============================================================================
# LOAD DATA
# =============================================================================

so.pk <- readRDS("./Analysis/R analysis/RDS files/S.O.Pk.integrated.rds")
so.gam <- subset(so.pk, idents = c(6,7))
so.gam$og_cluster <- so.gam$seurat_clusters

# =============================================================================
# PREPATE DATA
# =============================================================================

so.gam <- NormalizeData(so.gam)
so.gam <- FindVariableFeatures(so.gam)
so.gam <- ScaleData(so.gam)
so.gam <- RunPCA(so.gam, verbose = FALSE)
ElbowPlot(so.gam, ndims = 50)
so.gam <- RunUMAP(so.gam, dims=1:20, metric = 'correlation', verbose = FALSE, return.model = TRUE)
so.gam <- FindNeighbors(so.gam, dims = 1:20)
so.gam <- FindClusters(so.gam, algorithm = 4, random.seed = 42, resolution = 0.3)
DimPlot(so.gam, reduction = 'pca')
DimPlot(so.gam, reduction = 'umap')

# ============================================================================
# RUN PSEUDOTIME
# ============================================================================

Idents(so.gam) <- "seurat_clusters"
sling_data <- slingshot(Embeddings(so.gam, "umap"), clusterLabels = so.gam$seurat_clusters, start.clus = "1")
sling_data <- SlingshotDataSet(sling_data)
cell_id <- as.data.frame(slingBranchID(sling_data))
pseudotime_df <- as.data.frame(slingPseudotime(sling_data))
pseudotime_df$cell <- rownames(pseudotime_df)

# Adds pseudotime metadata to original object
pt_vector <- setNames(pseudotime_df$Lineage1, pseudotime_df$cell)
so.gam$fem_pseudotime <- NA  
cells_to_update <- intersect(colnames(so.gam), names(pt_vector))
so.gam$fem_pseudotime[match(cells_to_update, colnames(so.gam))] <- pt_vector[cells_to_update]

pt_vector <- setNames(pseudotime_df$Lineage2, pseudotime_df$cell)
so.gam$male_pseudotime <- NA  
cells_to_update <- intersect(colnames(so.gam), names(pt_vector))
so.gam$male_pseudotime[match(cells_to_update, colnames(so.gam))] <- pt_vector[cells_to_update]


# Creates trajectory plots
my_cols <- c("#440154FF", "#FFA500")
wrap_plots(
  FeaturePlot(so.gam, features = "fem_pseudotime", reduction = "umap", cols = my_cols, pt.size = 1.3), 
  FeaturePlot(so.gam, features = "male_pseudotime", reduction = "umap", cols = my_cols, pt.size = 1.3),
  ncol = 2
)

saveRDS(so.gam, "./Paper/pk/S.O.Pk.integrated_gametocytes.rds")

# ============================================================================
# FIND DEGs 
# ============================================================================

deg_clusters <- FindAllMarkers(so.gam, test.use = "wilcox", only.pos = TRUE)
deg_clusters_filtered <- deg_clusters[deg_clusters$p_val_adj < 0.05,]
all_markers <- add_gene_descriptions_from_gff_fast(deg_clusters_filtered, 
                                                   "./Analysis/Genomes/Pk-68/PlasmoDB-68_PknowlesiH.gff")
all_markers$gene <- gsub("-", "_", all_markers$gene)
write.xlsx(all_markers, file = "cluster_markers_gam.xlsx", rowNames = FALSE)

# ============================================================================
# GAM SPECIFIC MARKERS
# ============================================================================

gam_markers1 <- c(
  'PKNH-0418700' = 'g27',
  'PKNH-0414600' = 'arom',
  'PKNH-0514800' = 'metiltransferase',
  'PKNH-0314600' = 'nek4',
  'PKNH-1271000' = 'amp',
  'PKNH-1211100' = 'mget',
  'PKNH-0304300' = 's16',
  'PKNH-1213900' = 'l22',
  'PKNH-1254400' = 'p45/48',
  'PKNH-0615700' = 'p25',
  'PKNH-1230500' = 'lap5',
  'PKNH-1310200' = 'zinc finger'
)

gam_markers2 <- c(
  'PKNH-1254400' = 'p45/48',
  'PKNH-0615700' = 'p25',
  'PKNH-1230500' = 'lap5',
  'PKNH-1310200' = 'zinc finger'
)

genes <- names(gam_markers1)[names(gam_markers1) %in% rownames(so.pk)]

plots <- FeaturePlot(so.pk, features = genes, combine = FALSE, dims = c(2,3), order = T)
plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] +
    ggtitle(paste0(gsub("-", "_", genes[i]), " (", gam_markers1[[genes[i]]], ")"))
})
wrap_plots(plots, ncol = 4)

# ============================================================================
# TOP MOLECULAR MARKERS
# ============================================================================

markers <- read.xlsx("cluster_markers_gam.xlsx") %>%
  filter(p_val_adj < 0.05, avg_log2FC >= 0.2)
markers$gene <- gsub('_', '-', markers$gene)

# Asexual cells in the full object
Idents(so.pk) <- "stage"
asex_cells <- colnames(so.pk)[!so.pk@meta.data[["stage"]] %in% c("Gametocyte 1 (putative)", "Gametocyte 2 (putative)")]

# Candidate genes present in both objects
cand <- intersect(unique(markers$gene), rownames(so.pk))
cand <- intersect(cand, rownames(so.gam))

# Expression in asexual stages: % of cells expressing + mean normalized expression
counts_asex <- GetAssayData(so.pk, layer = "data")[cand, asex_cells, drop = FALSE]
asex_stats <- tibble(
  gene      = cand,
  pct_asex  = Matrix::rowMeans(counts_asex > 0),
  mean_asex = Matrix::rowMeans(counts_asex)
)

# Expression in each gametocyte cluster
Idents(so.gam) <- "seurat_clusters"
data_gam <- GetAssayData(so.gam, layer = "data")[cand, , drop = FALSE]
clust    <- so.gam@meta.data[["seurat_clusters"]]

gam_stats <- map_dfr(levels(factor(clust)), function(cl) {
  cells <- colnames(so.gam)[clust == cl]
  tibble(
    cluster  = cl,
    gene     = cand,
    pct_gam  = Matrix::rowMeans(data_gam[, cells, drop = FALSE] > 0),
    mean_gam = Matrix::rowMeans(data_gam[, cells, drop = FALSE])
  )
})

# Combine with the marker table, filter, and rank
res <- markers %>%
  mutate(cluster = as.character(cluster)) %>%
  inner_join(gam_stats, by = c("cluster", "gene")) %>%
  inner_join(asex_stats, by = "gene") %>%
  filter(pct_gam >= 0.4, pct_asex <= 0.15) %>%
  mutate(specificity = mean_gam - mean_asex) %>%          
  arrange(cluster, desc(specificity), desc(avg_log2FC))

top_genes <- res %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup()

gene_list <- split(top_genes$gene, top_genes$cluster)

