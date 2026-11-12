library(Seurat)
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(BiocParallel)
library(ggplot2)
library(Matrix)

PARENT <- "assay_RNA_50_Integrated_CellRanger_50PC_0.5res.rds"
OUT_SUBSET <- "SC_Sorghum_Xy_Phl_TA_Subset.rds"
OUT_FIT <- "SCE_fitted_Sorghum_Xy_Phl_TA_SubsetTradeseq.rds"
CLUSTERS_KEEP <- c("8", "13", "4", "7", "14", "5", "10")
NPCS <- 50
DIMS <- 1:42
RESOLUTION <- 0.3
START_CLUSTER <- "3"
N_KNOTS <- 8
SEED <- 42
EXPECTED_LINEAGES <- list(
  Lineage1 = c("3", "7", "2", "0", "1"),
  Lineage2 = c("3", "7", "2", "0", "5"),
  Lineage3 = c("3", "7", "2", "6"),
  Lineage4 = c("3", "7", "2", "4")
)

if (any(file.exists(c(OUT_SUBSET, OUT_FIT)))) {
  stop("Output RDS files already exist. Run in a fresh working directory to preserve the previous analysis.")
}

cds <- readRDS(PARENT)
stopifnot(all(CLUSTERS_KEEP %in% as.character(cds$seurat_clusters)))
cds$parent_cluster <- as.character(cds$seurat_clusters)
cds_sub <- subset(cds, subset = seurat_clusters %in% CLUSTERS_KEEP)
if (ncol(cds_sub) != 4001) {
  stop("Expected 4,001 vascular-subset cells; check the parent object and its cluster identities.")
}
rm(cds)
invisible(gc())

DefaultAssay(cds_sub) <- "RNA"
if (inherits(cds_sub[["RNA"]], "Assay5")) {
  cds_sub <- JoinLayers(cds_sub, assay = "RNA")
}
cds_sub <- NormalizeData(cds_sub, normalization.method = "LogNormalize", scale.factor = 10000)
cds_sub <- FindVariableFeatures(cds_sub, selection.method = "vst", nfeatures = 2000)
cds_sub <- ScaleData(cds_sub)
cds_sub <- RunPCA(cds_sub, npcs = NPCS, seed.use = SEED)
cds_sub <- FindNeighbors(cds_sub, dims = DIMS)
cds_sub <- FindClusters(cds_sub, resolution = RESOLUTION, random.seed = 0)
cds_sub <- RunUMAP(cds_sub, dims = DIMS, seed.use = SEED)

if (!setequal(as.character(cds_sub$seurat_clusters), as.character(0:7))) {
  stop("Expected subset clusters 0–7. Validate cluster identities before using the manuscript lineage labels.")
}

print(DimPlot(cds_sub, reduction = "umap", label = TRUE,
              repel = TRUE, pt.size = 0.4) + theme_classic())
print(FeaturePlot(cds_sub, features = "Sobic.010G017400.v3.2"))

UMAP <- Embeddings(cds_sub, reduction = "umap")
sling_res <- slingshot(
  UMAP,
  clusterLabels = cds_sub$seurat_clusters,
  start.clus = START_CLUSTER,
  omega = FALSE
)

lineages <- slingLineages(sling_res)
endpoints <- vapply(lineages, function(x) tail(x, 1), character(1))
expected_endpoints <- vapply(EXPECTED_LINEAGES, function(x) tail(x, 1), character(1))
if (length(lineages) != 4L || anyDuplicated(endpoints) ||
    !setequal(endpoints, expected_endpoints)) {
  stop("Expected four lineages ending in clusters 1, 5, 6, and 4. Check the inferred trajectories before downstream testing.")
}
lineage_order <- match(expected_endpoints, endpoints)
lineages <- lineages[lineage_order]
names(lineages) <- names(EXPECTED_LINEAGES)
if (!identical(unname(lapply(lineages, as.character)), unname(EXPECTED_LINEAGES))) {
  stop("Inferred cluster paths differ from the manuscript lineage mapping. Check cluster identities and trajectories.")
}

pseudotime <- slingPseudotime(sling_res, na = FALSE)[, lineage_order, drop = FALSE]
cellWeights <- slingCurveWeights(sling_res)[, lineage_order, drop = FALSE]
colnames(pseudotime) <- colnames(cellWeights) <- names(EXPECTED_LINEAGES)
counts_mt <- GetAssayData(cds_sub, assay = "RNA", layer = "counts")
counts_mt <- counts_mt[Matrix::rowSums(counts_mt) > 0, , drop = FALSE]
stopifnot(setequal(colnames(counts_mt), rownames(pseudotime)),
          setequal(colnames(counts_mt), rownames(cellWeights)))
pseudotime <- pseudotime[colnames(counts_mt), , drop = FALSE]
cellWeights <- cellWeights[colnames(counts_mt), , drop = FALSE]
stopifnot(identical(colnames(counts_mt), rownames(pseudotime)),
          identical(colnames(counts_mt), rownames(cellWeights)),
          all(is.finite(pseudotime)), all(is.finite(cellWeights)))

analysis_info <- list(
  parent_file = PARENT,
  parent_clusters = CLUSTERS_KEEP,
  n_cells = ncol(cds_sub),
  assay = "RNA", npcs = NPCS, dims = DIMS, resolution = RESOLUTION,
  start_cluster = START_CLUSTER, lineages = lineages,
  nknots = N_KNOTS, seed = SEED,
  session_info = sessionInfo()
)
cds_sub@misc$sorghum_analysis <- analysis_info
saveRDS(cds_sub, OUT_SUBSET)
write.csv(data.frame(
  lineage = names(lineages),
  clusters = vapply(lineages, paste, character(1), collapse = " -> "),
  endpoint = unname(expected_endpoints)
), "Sorghum_vascular_lineages.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "Sorghum_trajectory_sessionInfo.txt")

BPPARAM <- if (.Platform$OS.type == "windows") {
  BiocParallel::SnowParam(workers = 4, RNGseed = SEED)
} else {
  BiocParallel::MulticoreParam(workers = 4, RNGseed = SEED)
}
set.seed(SEED)
png("SC_Sorghum_Xy_Phl_TA_Subset_Tradeseq_knots.png", width = 2000, height = 3000, res = 300)
icMat <- evaluateK(
  counts = counts_mt, pseudotime = pseudotime, cellWeights = cellWeights,
  k = 4:12, nGenes = 200, verbose = TRUE
)
dev.off()

set.seed(SEED)
sce_fitted <- fitGAM(
  counts = counts_mt, pseudotime = pseudotime, cellWeights = cellWeights,
  nknots = N_KNOTS, verbose = TRUE, parallel = TRUE, BPPARAM = BPPARAM
)
metadata(sce_fitted)$sorghum_analysis <- analysis_info
stopifnot(identical(colnames(sce_fitted), colnames(cds_sub)))
saveRDS(sce_fitted, OUT_FIT)
