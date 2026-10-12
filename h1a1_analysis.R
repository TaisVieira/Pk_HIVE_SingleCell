library(Seurat)
library(ggplot2)
library(patchwork)

so.h1a1 <- readRDS("C:/Users/taisb/OneDrive/Desktop/Analysis/Paper/pk/S.O.Pk_erin.rds")

# Check AP2 expression
ap2 <- c(
  'PKNH-1441700' = 'AP2-G',
  'PKNH-1350200' = 'AP2-G2'
)

genes <- names(ap2)[names(ap2) %in% rownames(so.h1a1)]
plots <- FeaturePlot(so.h1a1, features = genes, combine = FALSE, order = T,
                     min.cutoff = "q97")
plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] +
    ggtitle(paste0(gsub("-", "_", genes[i]), " (", ap2[[genes[i]]], ")"))
})
wrap_plots(plots, ncol = 2)

# Create AP2-G Up dataset
genes <- c("PKNH-1441700", "PKNH-1350200")
expr <- FetchData(so.h1a1, vars = genes, layer = "data") 
cutoffs <- sapply(genes, function(g) quantile(expr[[g]], probs = 0.97, names = FALSE))
names(cutoffs)   
keep <- Reduce(`|`, lapply(genes, function(g) expr[[g]] >= cutoffs[[g]]))
sum(is.na(keep)); sum(keep)   
cells.keep <- rownames(expr)[keep]
ap2_up <- subset(so.h1a1, cells = cells.keep)
asexual <- subset(so.h1a1, cells = rownames(expr)[!keep])


#saveRDS(ap2_up, "C:/Users/taisb/OneDrive/Desktop/Analysis/Paper/pk/S.O.Pk_AP2G_Up_Erin.rds")
#saveRDS(asexual, "C:/Users/taisb/OneDrive/Desktop/Analysis/Paper/pk/S.O.Pk_asexual_Erin.rds")

# =============================================================================
# SUPPLEMENTARY FIGURE 8
# =============================================================================

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

FeaturePlot(so.h1a1, genes)






