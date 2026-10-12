library(Seurat)
library(patchwork)
library(pheatmap)
library(ggplot2)

setwd('./Analysis/R analysis/Scripts')
source('./util_funcs.R')
source('./pseudotime_funcs.R')
setwd('./Analysis/Paper/Pk_HIVE')

# =============================================================================
# LOAD DATA
# =============================================================================

so.pk <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk.integrated.rds")

# =============================================================================
# AP2-G UP DATASET
# =============================================================================

# Check AP2 expression
ap2 <- c(
  'PKNH-1441700' = 'AP2-G',
  'PKNH-1350200' = 'AP2-G2'
)

genes <- names(ap2)[names(ap2) %in% rownames(so.pk)]
plots <- FeaturePlot(so.pk, features = genes, combine = FALSE, dims = c(2,3), order = T,
                     min.cutoff = "q97")
plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] +
    ggtitle(paste0(gsub("-", "_", genes[i]), " (", ap2[[genes[i]]], ")"))
})
wrap_plots(plots, ncol = 2)

# Create AP2-G Up dataset
genes <- c("PKNH-1441700", "PKNH-1350200")
expr <- FetchData(so.pk, vars = genes, layer = "data") 
cutoffs <- sapply(genes, function(g) quantile(expr[[g]], probs = 0.97, names = FALSE))
names(cutoffs)   
keep <- Reduce(`|`, lapply(genes, function(g) expr[[g]] >= cutoffs[[g]]))
sum(is.na(keep)); sum(keep)   
cells.keep <- rownames(expr)[keep]
ap2_up <- subset(so.pk, cells = cells.keep)
asexual <- subset(so.pk, cells = rownames(expr)[!keep])
asexual <- subset(asexual, idents = c(1:5))

so.pk$AP2_high <- ifelse(colnames(so.pk) %in% colnames(ap2_up), "high", "rest")
VlnPlot(so.pk, features = genes, group.by = "AP2_high", pt.size = 0)
table(subset_obj$stage) / table(so.pk$stage)

stage_map <- c(
  "1" = "Trophozoite",
  "2" = "Early Schizont",
  "3" = "Late Schizont",
  "4" = "Late Schizont",
  "5" = "Ring",
  "6" = "Gametocyte 1 (putative)",
  "7" = "Gametocyte 2 (putative)"
)
ap2_up$stages <- unname(stage_map[as.character(ap2_up$seurat_clusters)])

stage_order <- c(
  "Ring",
  "Trophozoite",
  "Early Schizont",
  "Late Schizont",
  "Gametocyte 1 (putative)",
  "Gametocyte 2 (putative)"
)
ap2_up$stages <- factor(ap2_up$stages, levels = stage_order)

stage_colors <- c(
  "Ring"                  = "#ed7912",
  "Trophozoite"           = "#0c0da0",
  "Early Schizont"        = "#ad0825",
  "Late Schizont"         = "#5e1b83",
  "Gametocyte 1 (putative)" = "#164a11ff",
  "Gametocyte 2 (putative)" = "#e9db30"
)
DimPlot(ap2_up, group.by = "stages", cols = stage_colors, dims = c(2,3))

#saveRDS(ap2_up, "./Analysis/Paper/Pk_HIVE/S.O.Pk_AP2G_Up.rds")
#saveRDS(asexual, "./Analysis/Paper/Pk_HIVE/S.O.Pk_asexual.rds")

# =============================================================================
# COMPARATIVE ANALYSIS
# =============================================================================

so.gam <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk.integrated_gametocytes.rds")
ap2_up <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk_AP2G_Up.rds")
asexual <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk_asexual.rds")

DefaultAssay(so.gam)  <- "RNA"
DefaultAssay(ap2_up)  <- "RNA"
DefaultAssay(asexual) <- "RNA"

get_top_markers <- function(obj, n = 6) {
  FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.25) %>%
    group_by(cluster) %>%
    top_n(n, avg_log2FC) %>%
    pull(gene) %>%
    unique()
}

top_ap2up   <- get_top_markers(ap2_up)
top_asexual <- get_top_markers(asexual)
top_gam     <- get_top_markers(so.gam)

# Classify AP2up genes 
genes_exclusive   <- setdiff(top_ap2up, union(top_asexual, top_gam))
genes_com_asexual <- intersect(top_ap2up, top_asexual)
genes_com_gam     <- intersect(top_ap2up, top_gam)

genes_final <- unique(c(genes_exclusive, genes_com_asexual, genes_com_gam))
genes_final <- genes_final[genes_final %in% rownames(ap2_up) &
                             genes_final %in% rownames(asexual) &
                             genes_final %in% rownames(so.gam)]

# Average expression per cluster 
get_avg_exp <- function(obj, genes, group.by = "seurat_clusters") {
  avg <- AverageExpression(obj, features = genes, group.by = group.by,
                           assays = "RNA")$RNA
  avg[genes, , drop = FALSE]
}

avg_ap2up   <- get_avg_exp(ap2_up,  genes_final)
avg_asexual <- get_avg_exp(asexual, genes_final)
avg_gam     <- get_avg_exp(so.gam,  genes_final)

# Final matrix: mean across clusters of each object
heatmap_matrix <- cbind(
  AP2up     = rowMeans(avg_ap2up),
  Asexual   = rowMeans(avg_asexual),
  Gametocytes = rowMeans(avg_gam)
)

# Row z-score and plot
heatmap_scaled <- t(scale(t(heatmap_matrix)))

mat <- t(heatmap_scaled)[c("Gametocytes", "AP2up", "Asexual"), ]
rownames(mat) <- c("Gametocytes", "AP2g- up", "Asexual")

col_fun <- colorRampPalette(c("black", "#2B0A5E", "#8A4FFF", "#F0669B",
                              "#F58A2E", "#B01E0E", "#7A0A0A"))(100)

pheatmap(mat,
         name           = "Z-score",
         color          = col_fun,
         breaks         = seq(-1, 1, length.out = 101),  
         cluster_rows   = FALSE,
         cluster_cols   = TRUE,
         treeheight_row = 0,
         treeheight_col = 0,
         border_color   = NA,
         angle_col      = "90",
         fontsize_col   = 8,
         fontsize_row   = 10,
         cellheight     = 30,
         main           = NA)

# =============================================================================
# PSEUDOTIME EXPRESSION
# =============================================================================

# genes chosen manually from the heatmap 
genes <- c(
  "PKNH-1353400",
  "PKNH-0817000",
  "PKNH-0615700",
  "PKNH-0615600",
  "PKNH-1245800",
  "PKNH-1411500",
  "PKNH-1324400",
  "PKNH-1300900",
  "PKNH-0600500",
  "PKNH-0937800"
)


plotPseudotimeGenes <- function(genes,
                                objs,
                                pseudotime_cols,
                                ncol = 4,
                                title = NULL,
                                use_rank = TRUE) {
    
  plots <- list()
  
  for (gene in genes) {
    for (obj_name in names(objs)) {
      obj <- objs[[obj_name]]
      pt_col <- pseudotime_cols[[obj_name]]
      
      df <- FetchData(obj, vars = c(pt_col, gene))
      colnames(df) <- c("Pseudotime", "Expression")
      df <- df[is.finite(df$Pseudotime), , drop = FALSE]
      
      # Replace pseudotime values by cell order (1 = earliest cell)
      if (use_rank) {
        df$Pseudotime <- rank(df$Pseudotime, ties.method = "first")
      }
      
      p <- ggplot(df, aes(x = Pseudotime, y = Expression)) +
        geom_point(size = 0.3, alpha = 0.15, stroke = 0) +
        geom_smooth(method = "loess", span = 0.3, color = "#f00f90", se = FALSE) +
        labs(title = paste0(gene, " (", obj_name, ")"),
             x = if (use_rank) "Cell order (by pseudotime)" else "Pseudotime",
             y = "Norm. expr.") +
        theme_minimal(base_size = 10)
      
      plots[[paste(gene, obj_name, sep = "_")]] <- p
    }
  }
  
  out <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(title = title)
  
  attr(out, "n_panels") <- length(plots)   
  out
}

# Pseudotime column for each object 
pt_cols <- c(`All cells` = "idc_pseudotime",
             Asexual     = "idc_pseudotime",    
             `AP2-G-Up`  = "idc_pseudotime",
             Gametocyte = "gam_pseudotime")


jobs <- list(
  full = list(genes = genes,
              objs  = list(`All cells` = so.pk),
              title = "All cells"),
  asexual  = list(genes = genes,
                  objs  = list(Asexual = asexual),
                  title = "Asexual"),
  ap2_up  = list(genes = genes,
                        objs  = list(`AP2-G-Up` = ap2_up),
                        title = "AP2-G-Up"),
  gametocyte  = list(genes = genes,
                 objs  = list(Gametocyte = so.gam),
                 title = "Gametocyte")
)

plots <- lapply(jobs, function(j) {
  plotPseudotimeGenes(j$genes, j$objs,
                      pseudotime_cols = pt_cols,
                      ncol = 4, title = j$title)
})


save_panels <- function(p, filename, ncol = 4,
                        panel_w = 3.5, panel_h = 2.8,   
                        title_h = 0.6, dpi = 300) {
  n_panels <- attr(p, "n_panels")
  nrows <- ceiling(n_panels / ncol)
  ncols <- min(n_panels, ncol)
  
  ggsave(filename, p,
         width  = ncols * panel_w,
         height = nrows * panel_h + title_h,
         dpi = dpi, limitsize = FALSE)
}

for (n in names(plots)) {
  save_panels(plots[[n]], paste0(n, "_pseudotime.png"), ncol = 4)
}

FeaturePlot(so.pk, genes, order = T, dims = c(2,3))
