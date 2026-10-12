library(Seurat)
library(slingshot)
library(patchwork)
library(ggplot2)
library(openxlsx)
library(dplyr)

setwd('./Analysis/R analysis/Scripts')
source('./util_funcs.R')
source('./analyzeSlingshot.R')
setwd('./Analysis/Paper/pk')

# =============================================================================
# LOAD DATA
# =============================================================================

so.pk <- readRDS("./Analysis/Paper/pk/S.O.Pk.integrated.rds")

# ============================================================================
# FIND DEGs 
# ============================================================================

deg_clusters <- FindAllMarkers(so.pk, test.use = "wilcox", only.pos = TRUE)
deg_clusters_filtered <- deg_clusters[deg_clusters$p_val_adj < 0.05,]
all_markers <- add_gene_descriptions_from_gff_fast(deg_clusters_filtered, 
                                                   "./Analysis/Genomes/Pk-68/PlasmoDB-68_PknowlesiH.gff")
all_markers$gene <- gsub("-", "_", all_markers$gene)
write.xlsx(all_markers, file = "cluster_markers.xlsx", rowNames = FALSE)

# ============================================================================
# RUN PSEUDOTIME
# ============================================================================

Idents(so.pk) <- "seurat_clusters"
sling_data <- slingshot(Embeddings(so.pk, "umap"), clusterLabels = so.pk$seurat_clusters, start.clus = "5")
sling_data <- SlingshotDataSet(sling_data)
cell_id <- as.data.frame(slingBranchID(sling_data))
pseudotime_df <- as.data.frame(slingPseudotime(sling_data))
pseudotime_df$cell <- rownames(pseudotime_df)

# Adds pseudotime metadata to original object
pt_vector <- setNames(pseudotime_df$Lineage1, pseudotime_df$cell)
so.pk$gam_pseudotime <- NA  
cells_to_update <- intersect(colnames(so.pk), names(pt_vector))
so.pk$gam_pseudotime[match(cells_to_update, colnames(so.pk))] <- pt_vector[cells_to_update]

# Subsets original object to just IDC so gametocytes don't influence trajectory
so.idc <- subset(so.pk, idents = c(1:5))
sling_data <- slingshot(Embeddings(so.idc, "umap"), clusterLabels = so.idc$seurat_clusters, start.clus = "5")
sling_data <- SlingshotDataSet(sling_data)
cell_id <- as.data.frame(slingBranchID(sling_data))
pseudotime_df <- as.data.frame(slingPseudotime(sling_data))
pseudotime_df$cell <- rownames(pseudotime_df)

pt_vector <- setNames(pseudotime_df$Lineage1, pseudotime_df$cell)
so.pk$idc_pseudotime <- NA  
cells_to_update <- intersect(colnames(so.pk), names(pt_vector))
so.pk$idc_pseudotime[match(cells_to_update, colnames(so.pk))] <- pt_vector[cells_to_update]


# Creates trajectory plots
my_cols <- c("#440154FF", "#FFA500")
wrap_plots(
  FeaturePlot(so.pk, features = "idc_pseudotime", reduction = "umap", cols = my_cols, pt.size = .5, dims = c(2,3)),
  FeaturePlot(so.pk, features = "gam_pseudotime", reduction = "umap", cols = my_cols, pt.size = .5, dims = c(2,3)),
  ncol = 2
)

# ============================================================================
# ASSIGN STAGES TO CLUSTERS
# ============================================================================

stage_map <- c(
  "1" = "Trophozoite",
  "2" = "Early Schizont",
  "3" = "Late Schizont",
  "4" = "Late Schizont",
  "5" = "Ring",
  "6" = "Gametocyte 1 (putative)",
  "7" = "Gametocyte 2 (putative)"
)

so.pk$stage <- unname(stage_map[as.character(so.pk$seurat_clusters)])

#saveRDS(so.pk, "./Analysis/Paper/pk/S.O.Pk.integrated.rds")

# ============================================================================
# TOP 5 MARKERS
# ============================================================================

all_markers <- read.xlsx('cluster_markers.xlsx')

stage_levels <- c("Ring", "Trophozoite", "Early Schizont", "Late Schizont",
                  "Gametocyte 1 (putative)", "Gametocyte 2 (putative)")
so.pk$stage <- factor(so.pk$stage, levels = stage_levels)
Idents(so.pk) <- "stage"

top5_by_cluster <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  ungroup() %>%
  mutate(stage = stage_map[as.character(cluster)]) %>%
  mutate(stage = factor(stage, levels = stage_levels)) %>%
  arrange(stage, cluster, desc(avg_log2FC))

top_genes_unique <- unique(gsub("_", "-", top5_by_cluster$gene))

DotPlot(so.pk, features = top_genes_unique, group.by = "stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Top 5 DEGs per cluster, grouped by cluster
genes <- read_csv("cluster_markers.csv") %>%
  arrange(cluster) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  pull(gene) %>%
  unique() 

genes <- gsub('_', '-', genes)

# Average expression per cluster, scaled per gene
avg <- AverageExpression(so.pk, features = genes, group.by = "seurat_clusters")$RNA
avg <- t(scale(t(avg[genes, ])))

pheatmap(avg, cluster_rows = FALSE, cluster_cols = FALSE, #border_color = NA, 
         color = colorRampPalette(c("navy", "white", "red"))(100))

