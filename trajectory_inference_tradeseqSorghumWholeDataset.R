library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(Matrix)
library(Seurat)

cds <- assay_RNA_50_Integrated_CellRanger_50PC_0.5res

cds <- readRDS("SC_Sorghum_Xy_Phl_TA_Subset.rds")

UMAP <- Seurat::Embeddings( cds[["umap"]] )

start_clust <- "8"

sling_res <- slingshot::slingshot( UMAP,
                                   clusterLabels = cds@meta.data$seurat_clusters,
                                   start.clus = start_clust,
                                   omega = F,

)

slingshot::slingLineages(sling_res)

DimPlot(
  object = cds,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
)

clusters_remove <- c("5", "0", "3")

cds_sub <- subset(
  cds,
  subset = !seurat_clusters %in% clusters_remove
)

DimPlot(
  cds_sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 0.4,
  label = TRUE,
  repel = TRUE
)

UMAP <- Seurat::Embeddings( cds_sub[["umap"]] )

start_clust <- "8"

sling_res <- slingshot::slingshot( UMAP,
                                   clusterLabels = cds_sub@meta.data$seurat_clusters,
                                   start.clus = start_clust,
                                   omega = F,

)

slingshot::slingLineages(sling_res)

clusters_keep <- c("8", "13", "4", "7", "14", "10")

cds_sub <- subset(cds, subset = seurat_clusters %in% clusters_keep)

cds_sub <- NormalizeData(cds_sub)
cds_sub <- FindVariableFeatures(cds_sub)
cds_sub <- ScaleData(cds_sub)
cds_sub <- RunPCA(cds_sub, npcs = 50)

cds_sub <- FindNeighbors(
  cds_sub,
  dims = 1:30,
  graph.name = "RNA_snn"
)

cds_sub <- FindClusters(
  cds_sub,
  resolution = 0.3,
  graph.name = "RNA_snn"
)

cds_sub <- RunUMAP(cds_sub, dims = 1:30)

DimPlot(
  cds_sub,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  theme_classic()

FeaturePlot(cds_sub,features = "Sobic.010G017400.v3.2")

UMAP <- Seurat::Embeddings( cds_sub[["umap"]] )

start_clust <- "C3"

sling_res <- slingshot::slingshot( UMAP,
                                   clusterLabels = cds_sub@meta.data$seurat_clusters,
                                   start.clus = start_clust,
                                   omega = F,

)

slingshot::slingLineages(sling_res)

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 4

counts_mt <- cds@assays$RNA$counts

counts_mt <- counts_mt[rowSums(counts_mt) > 0, ]

png("SC_Sorghum_Xy_Phl_TA_Subset_Tradeseq_knots.png", width = 2000, height = 3000, res = 300)

icMat <- evaluateK(counts = counts_mt, sds = sling_res, k = 4:12,
                   nGenes = 200, verbose = T)
dev.off()

sce_fitted <- tradeSeq::fitGAM(
  counts = counts_mt,
  sds = sling_res,
  nknots = 8,
  verbose = TRUE,
  parallel = TRUE,
  BPPARAM = BPPARAM
)

saveRDS(sce_fitted, file = "SCE_fitted_Sorghum_Xy_Phl_TA_SubsetTradeseq.rds")
