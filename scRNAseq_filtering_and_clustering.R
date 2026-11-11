library(Seurat)

N609_D29 <- Read10X(data.dir = "filtered_feature_bc_matrix_N609_D29")
N609_D29 <- CreateSeuratObject(counts = N609_D29, project = "N609_D29", min.cells = 3, min.features = 200)
N609_D29 <- subset(N609_D29, subset = nFeature_RNA > 700 & nFeature_RNA < 3500)

N609_D29 <- NormalizeData(N609_D29, normalization.method = "LogNormalize", scale.factor = 10000)
N609_D29 <- FindVariableFeatures(N609_D29, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(N609_D29)
N609_D29 <- ScaleData(N609_D29, features = all.genes)
N609_D29 <- RunPCA(N609_D29, features = VariableFeatures(object = N609_D29))

N609_D31 <- Read10X(data.dir = "filtered_feature_bc_matrix_N609_D31")
N609_D31 <- CreateSeuratObject(counts = N609_D31, project = "N609_D31", min.cells = 3, min.features = 200)
N609_D31 <- subset(N609_D31, subset = nFeature_RNA > 700 & nFeature_RNA < 3500)

N609_D31 <- NormalizeData(N609_D31, normalization.method = "LogNormalize", scale.factor = 10000)
N609_D31 <- FindVariableFeatures(N609_D31, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(N609_D31)
N609_D31 <- ScaleData(N609_D31, features = all.genes)
N609_D31 <- RunPCA(N609_D31, features = VariableFeatures(object = N609_D31))

Tx430_D29 <- Read10X(data.dir = "filtered_feature_bc_matrix_Tx430_D29")
Tx430_D29 <- CreateSeuratObject(counts = Tx430_D29, project = "Tx430_D29", min.cells = 3, min.features = 200)
Tx430_D29 <- subset(Tx430_D29, subset = nFeature_RNA > 700 & nFeature_RNA < 3000)

Tx430_D29 <- NormalizeData(Tx430_D29, normalization.method = "LogNormalize", scale.factor = 10000)
Tx430_D29 <- FindVariableFeatures(Tx430_D29, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tx430_D29)
Tx430_D29 <- ScaleData(Tx430_D29, features = all.genes)
Tx430_D29 <- RunPCA(Tx430_D29, features = VariableFeatures(object = Tx430_D29))

Tx430_D31 <- Read10X(data.dir = "filtered_feature_bc_matrix_Tx430_D31")
Tx430_D31 <- CreateSeuratObject(counts = Tx430_D31, project = "Tx430_D31", min.cells = 3, min.features = 200)
Tx430_D31 <- subset(Tx430_D31, subset = nFeature_RNA > 700 & nFeature_RNA < 3500)

Tx430_D31 <- NormalizeData(Tx430_D31, normalization.method = "LogNormalize", scale.factor = 10000)
Tx430_D31 <- FindVariableFeatures(Tx430_D31, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tx430_D31)
Tx430_D31 <- ScaleData(Tx430_D31, features = all.genes)
Tx430_D31 <- RunPCA(Tx430_D31, features = VariableFeatures(object = Tx430_D31))

sing_cell_list <- list()

N609_D29_filtered <- N609_D29

N609_D29_filtered <- AddMetaData(N609_D29_filtered,
                                 "N609_D29_filtered",
                                 col.name = 'Sample')

N609_D29_filtered <- SCTransform(N609_D29_filtered, conserve.memory = T, verbose = T)

N609_D29_filtered <- FindVariableFeatures(N609_D29_filtered,
                                          selection.method = "vst",
                                          nfeatures = 2000,
                                          verbose = T)

sing_cell_list[[1]] <- N609_D29_filtered

N609_D31_filtered <- N609_D31

N609_D31_filtered <- AddMetaData(N609_D31_filtered,
                                 "N609_D31_filtered",
                                 col.name = 'Sample')

N609_D31_filtered <- SCTransform(N609_D31_filtered, conserve.memory = T, verbose = T)

N609_D31_filtered <- FindVariableFeatures(N609_D31_filtered,
                                          selection.method = "vst",
                                          nfeatures = 2000,
                                          verbose = T)

sing_cell_list[[2]] <- N609_D31_filtered

Tx430_D29_filtered <- Tx430_D29
Tx430_D29_filtered <- AddMetaData(Tx430_D29_filtered,
                                  "Tx430_D29_filtered",
                                  col.name = 'Sample')

Tx430_D29_filtered <- SCTransform(Tx430_D29_filtered, conserve.memory = T, verbose = T)

Tx430_D29_filtered <- FindVariableFeatures(Tx430_D29_filtered,
                                           selection.method = "vst",
                                           nfeatures = 2000,
                                           verbose = T)

sing_cell_list[[3]] <- Tx430_D29_filtered

Tx430_D31_filtered <- Tx430_D31
Tx430_D31_filtered <- AddMetaData(Tx430_D31_filtered,
                                  "Tx430_D31_filtered",
                                  col.name = 'Sample')

Tx430_D31_filtered <- SCTransform(Tx430_D31_filtered, conserve.memory = T, verbose = T)

Tx430_D31_filtered <- FindVariableFeatures(Tx430_D31_filtered,
                                           selection.method = "vst",
                                           nfeatures = 2000,
                                           verbose = T)

sing_cell_list[[4]] <- Tx430_D31_filtered

features <- SelectIntegrationFeatures(object.list = sing_cell_list)

sc_data.anchors <- FindIntegrationAnchors( object.list = sing_cell_list,
                                           anchor.features = features
)

sc_data.combined <- IntegrateData( anchorset = sc_data.anchors )

saveRDS(sc_data.combined, "Sorghum_CellRanger_filtered_sc_data.combined.rds")

Sorghum_CellRanger_filtered_sc_data.combined <- sc_data.combined
Sorghum_CellRanger_filtered_sc_data.combined <- FindVariableFeatures(Sorghum_CellRanger_filtered_sc_data.combined, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Sorghum_CellRanger_filtered_sc_data.combined)
Sorghum_CellRanger_filtered_sc_data.combined <- ScaleData(Sorghum_CellRanger_filtered_sc_data.combined, features = all.genes)

Sorghum_CellRanger_filtered_sc_data.combined <- RunPCA(Sorghum_CellRanger_filtered_sc_data.combined, features = VariableFeatures(object = Sorghum_CellRanger_filtered_sc_data.combined))
ElbowPlot(Sorghum_CellRanger_filtered_sc_data.combined, ndims = 50)

Sorghum_CellRanger_filtered_sc_data.combined <- FindNeighbors(Sorghum_CellRanger_filtered_sc_data.combined, dims = 1:50)

Sorghum_CellRanger_filtered_sc_data.combined <- FindClusters(Sorghum_CellRanger_filtered_sc_data.combined, resolution = 0.5)
Sorghum_CellRanger_filtered_sc_data.combined <- RunUMAP(Sorghum_CellRanger_filtered_sc_data.combined, dims = 1:50)

DimPlot(Sorghum_CellRanger_filtered_sc_data.combined, reduction = "umap", label = TRUE, label.size = 8)

saveRDS(Sorghum_CellRanger_filtered_sc_data.combined, file = 'assay_RNA_50_Integrated_CellRanger_50PC_0.5res.rds')
