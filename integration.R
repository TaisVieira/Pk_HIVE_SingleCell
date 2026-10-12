library(Seurat)
library(harmony)

setwd('./Analysis/R analysis/Scripts')
source('./util_funcs.R')
setwd('./Analysis/Paper/pk')

# =============================================================================
# LOAD DATA
# =============================================================================

pk1 <- readRDS("./Analysis/R analysis/RDS files/S.O.Pk1.rds")
pk1$dataset <- 'pk1'

pk2 <- readRDS("./Analysis/R analysis/RDS files/S.O.Pk2.rds")
pk2$dataset <- 'pk2'

# ============================================================================
# HARMONY INTEGRATION 
# ============================================================================

merged_obj <- merge(
  x = pk1, 
  y = pk2,
  add.cell.ids = c("pk1", "pk2")
)
merged_obj <- JoinLayers(merged_obj)
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 1330)  # 25% of total genes
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, verbose = FALSE)

so.integrated <- RunHarmony(
  merged_obj,
  group.by.vars = 'dataset',
  dims.use = 1:20, 
  verbose = FALSE
)

so.integrated <- RunUMAP(
  so.integrated,
  reduction = "harmony",
  dims = 1:20,
  n.components = 3, 
  min.dist = 0.4, 
  spread = 2.5, 
  repulsion.strength = 2,
  metric = "correlation",
  verbose = FALSE,
  return.model = TRUE
)

so.integrated <- FindNeighbors(
  so.integrated, 
  reduction = 'harmony', 
  dims = 1:20
)

so.integrated <- FindClusters(
  so.integrated, 
  resolution = 0.4, 
  algorithm = 4, 
  random.seed = 42
)

DimPlot(so.integrated, dims=c(2,3))

create_3DPlot(so.integrated)

# ============================================================================
# SAVE RESULTS
# ============================================================================

saveRDS(
  so.integrated,
  file = "C:/Users/taisb/OneDrive/Desktop/Analysis/Paper/pk/S.O.Pk_final.rds"
)

