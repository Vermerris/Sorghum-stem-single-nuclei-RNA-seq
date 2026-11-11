library(Seurat)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(Matrix)
library(BiocParallel)

cds <- readRDS("assay_RNA_50_Integrated_CellRanger_50PC_0.5res.rds")

UMAP <- Seurat::Embeddings(cds[["umap"]])

start_clust <- "8"

sling_res <- slingshot(
  UMAP,
  clusterLabels = cds@meta.data$seurat_clusters,
  start.clus = start_clust,
  omega = FALSE
)

pseudotime <- slingPseudotime(sling_res)
cellWeights <- slingCurveWeights(sling_res)

pseudotime_l13 <- pseudotime[, 1:3, drop = FALSE]
cellWeights_l13 <- cellWeights[, 1:3, drop = FALSE]

valid_cells <- rowSums(!is.na(pseudotime_l13)) > 0 & rowSums(cellWeights_l13 > 0) > 0

pseudotime_sub <- pseudotime_l13[valid_cells, , drop = FALSE]
cellWeights_sub <- cellWeights_l13[valid_cells, , drop = FALSE]

counts <- cds@assays$RNA@counts
counts_sub <- counts[, valid_cells]

counts_sub <- counts_sub[rowSums(counts_sub) > 0, ]

BPPARAM <- BiocParallel::MulticoreParam(workers = 23)

sce_fitted <- fitGAM(
  counts = counts_sub,
  pseudotime = pseudotime_sub,
  cellWeights = cellWeights_sub,
  nknots = 6,
  verbose = TRUE,
  parallel = TRUE,
  BPPARAM = BPPARAM
)

saveRDS(sce_fitted, file = "Sorghum_Tradeseq_lineages1to3.rds")
