
library(Seurat)
library(ggplot2)
library(patchwork)

# =============================================================================
# LOAD DATA
# =============================================================================

so.pk <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk.integrated.rds")
so.h1a1 <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk_H1A1.rds")
ap2_up <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk_AP2G_Up.rds")
ap2_h1a1 <- readRDS("./Analysis/Paper/Pk_HIVE/S.O.Pk_AP2G_Up_H1A1.rds")

# =============================================================================
# COMPARE DATASETS
# =============================================================================

compare_expression <- function(obj1, obj2, fc_threshold = 10, top_n = 10,
                               pseudo = 1e-6, assay = NULL, layer = "data") {
  
  genes <- unique(intersect(rownames(obj1), rownames(obj2)))
  
  get_mat <- function(obj) {
    a <- if (is.null(assay)) DefaultAssay(obj) else assay
    GetAssayData(obj, assay = a, layer = layer)[genes, , drop = FALSE]  
  }
  
  avg1 <- Matrix::rowMeans(get_mat(obj1))
  avg2 <- Matrix::rowMeans(get_mat(obj2))
  
  fold_change <- (avg1 + pseudo) / (avg2 + pseudo)
  
  res <- data.frame(
    gene        = genes,
    avg_obj1    = as.numeric(avg1),
    avg_obj2    = as.numeric(avg2),
    fold_change = as.numeric(fold_change),
    log2FC      = log2(as.numeric(fold_change))
  )
  
  top1 <- res[res$fold_change >= fc_threshold, ]
  top1 <- head(top1[order(top1$avg_obj1, decreasing = TRUE), ], top_n)
  
  top2 <- res[res$fold_change <= 1 / fc_threshold, ]
  top2 <- head(top2[order(top2$avg_obj2, decreasing = TRUE), ], top_n)
  
  rownames(top1) <- rownames(top2) <- NULL
  list(higher_in_obj1 = top1, higher_in_obj2 = top2)
}

full_obj <- compare_expression(so.pk, so.h1a1)
ap2_obj <- compare_expression(ap2_up, ap2_h1a1)

combined <- rbind(
  cbind(direction = "higher_in_obj1", full_obj$higher_in_obj1),
  cbind(direction = "higher_in_obj2", full_obj$higher_in_obj2)
)
write.csv(combined, "expression_comparison_hvh1.csv", row.names = FALSE)

combined <- rbind(
  cbind(direction = "higher_in_obj1", ap2_obj$higher_in_obj1),
  cbind(direction = "higher_in_obj2", ap2_obj$higher_in_obj2)
)
write.csv(combined, "expression_comparison_ap2_hvh1.csv", row.names = FALSE)


plotCompareGenes <- function(genes,
                             obj1,
                             obj2,
                             layer = "data",
                             group_labels = c("Gameto", "Belgica"),
                             ncol = 4,
                             log2_transform = FALSE,
                             pseudo = 1e-6) {
  
  m1 <- GetAssayData(obj1, layer = layer)
  m2 <- GetAssayData(obj2, layer = layer)
  
  # Named vector: name = gene ID, value = display label. Unnamed: ID is the label.
  nm  <- names(genes)
  if (is.null(nm)) nm <- rep("", length(genes))
  val <- unname(genes)
  ids   <- ifelse(nm == "", val, nm)
  disps <- val
  
  find_id <- function(g, mat) {
    if (g %in% rownames(mat)) return(g)
    hit <- grep(g, rownames(mat), value = TRUE, fixed = TRUE)
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  build_gene_plot <- function(g, disp) {
    id1 <- find_id(g, m1)
    id2 <- find_id(g, m2)
    if (is.na(id1) || is.na(id2)) {
      warning("Gene not found, skipping: ", g)
      return(NULL)
    }
    
    v1 <- as.numeric(m1[id1, ])
    v2 <- as.numeric(m2[id2, ])
    
    # FC on the untransformed layer values (same logic as compare_expression)
    fc <- (mean(v1) + pseudo) / (mean(v2) + pseudo)
    fc_lab <- if (fc >= 1) sprintf("FC = %.1f", fc) else sprintf("FC = 1/%.1f", 1 / fc)
    
    if (log2_transform) { v1 <- log2(v1 + 1); v2 <- log2(v2 + 1) }
    
    df <- data.frame(
      expr  = c(v1, v2),
      group = factor(rep(group_labels, c(length(v1), length(v2))),
                     levels = group_labels)
    )
    
    pval <- tryCatch(wilcox.test(v1, v2)$p.value, error = function(e) NA_real_)
    lab <- paste0(
      fc_lab, "\n",
      if (!is.na(pval)) {
        if (pval < 1e-4) "p < 0.0001" else sprintf("p = %.4f", pval)
      } else "p = NA"
    )
    ytop <- max(df$expr, na.rm = TRUE)
    if (ytop == 0) ytop <- 1
    
    ggplot(df, aes(group, expr, color = group)) +
      geom_jitter(width = 0.18, height = 0, size = 1.4, alpha = 0.8) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.5,
                   linewidth = 0.4, color = "grey25") +
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                   geom = "errorbar", width = 0.25, color = "grey25") +
      annotate("text", x = 1.5, y = ytop * 1.12, label = lab,
               size = 3, lineheight = 0.9) +
      scale_color_manual(values = c("#f00f90", "#f777c0")) +
      #scale_color_manual(values = c("#b415ea", "#e4abf8")) +
      labs(title = disp, x = NULL,
           y = if (log2_transform) expression(Expression~(log[2])) else "Expression") +
      expand_limits(y = ytop * 1.2) +
      theme_classic(base_size = 11) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
            axis.text.x = element_text(size = 10))
  }
  
  plots <- Map(build_gene_plot, ids, disps)
  plots <- Filter(Negate(is.null), plots)
  wrap_plots(plots, ncol = ncol)
}

plotCompareGenes(c(full_obj$higher_in_obj1$gene, full_obj$higher_in_obj2$gene, "PKNH-1324400"), 
                 so.pk, so.h1a1, group_labels = c("H", "H1A1"))
plotCompareGenes(c(ap2_obj$higher_in_obj1$gene, ap2_obj$higher_in_obj2$gene, "PKNH-1324400"), 
                 ap2_up, ap2_h1a1, group_labels = c("H", "H1A1"))










