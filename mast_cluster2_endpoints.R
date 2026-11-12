if (!requireNamespace("MAST", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("MAST", ask = FALSE)
}
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")

library(Seurat)
library(SingleCellExperiment)
library(MAST)
library(ggplot2)
library(dplyr)
library(svglite)

OUTPUT_DIR <- "."
dir.create(OUTPUT_DIR, showWarnings = FALSE)

DPI       <- 600
FDR_CUT   <- 0.05
LFC_CUT   <- 0.5
MIN_PCT   <- 0.1
TEST_LFC  <- 0.25

lin_cb <- c(Lineage1 = "#0072B2", Lineage2 = "#D55E00",
            Lineage3 = "#009E73", Lineage4 = "#CC79A7")

export_plot <- function(p, base, w, h) {
  f <- file.path(OUTPUT_DIR, base)
  ggsave(paste0(f, ".svg"),  p, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(f, ".pdf"),  p, width = w, height = h)
  ggsave(paste0(f, ".png"),  p, width = w, height = h, dpi = DPI, bg = "white")
  ggsave(paste0(f, ".tiff"), p, width = w, height = h, dpi = DPI,
         bg = "white", compression = "lzw")
  cat(sprintf("  saved: %s (svg/pdf/png@%d/tiff@%d)\n", base, DPI, DPI))
}

cat("Loading fitted tradeSeq SCE...\n")
sce <- readRDS("SCE_fitted_Sorghum_Xy_Phl_TA_SubsetTradeseq.rds")
cat("Loading Seurat subset...\n")
seu <- readRDS("SC_Sorghum_Xy_Phl_TA_Subset.rds")
stopifnot(setequal(colnames(sce), colnames(seu)),
          identical(metadata(sce)$sorghum_analysis, seu@misc$sorghum_analysis))

DefaultAssay(seu) <- "RNA"

if (inherits(seu[["RNA"]], "Assay5") && length(Layers(seu[["RNA"]], search = "data")) > 1) {
  cat("Joining Seurat v5 layers...\n")
  seu <- JoinLayers(seu)
}

cd  <- as.data.frame(colData(sce))
wts <- as.matrix(cd[, grep("^crv\\.cellWeights\\.", names(cd)), drop = FALSE])
stopifnot(identical(sub("^crv\\.cellWeights\\.", "", colnames(wts)),
                    paste0("Lineage", 1:4)))
colnames(wts) <- sub("^crv\\.cellWeights\\.", "", colnames(wts))

shared <- intersect(colnames(seu), rownames(wts))
cat(sprintf("Shared barcodes (Seurat n Slingshot): %d\n", length(shared)))

cl <- as.character(seu$seurat_clusters); names(cl) <- colnames(seu)

on_lineage <- function(li) {
  w <- setNames(rep(0, ncol(seu)), colnames(seu))
  w[shared] <- wts[shared, li]
  names(w)[w > 0.5]
}
L <- setNames(lapply(seq_len(ncol(wts)), on_lineage),
              paste0("Lineage", seq_len(ncol(wts))))

C2 <- names(cl)[cl == "2"]
C6 <- names(cl)[cl == "6"]
C4 <- names(cl)[cl == "4"]

cat("\n=========== CLUSTER 2 COMPOSITION ===========\n")
cat(sprintf("Cluster 2 total cells : %d\n", length(C2)))
cat(sprintf("Cluster 6 total cells : %d\n", length(C6)))
cat(sprintf("Cluster 4 total cells : %d\n\n", length(C4)))

c2_by_lin <- sapply(L, function(x) length(intersect(C2, x)))
print(c2_by_lin)
cat(sprintf("\nSum across lineages: %d  (vs %d actual cells)\n",
            sum(c2_by_lin), length(C2)))

memb <- sapply(L, function(x) C2 %in% x)
rownames(memb) <- C2
n_assign <- rowSums(memb)
cat("\nLineages assigned per cluster-2 cell:\n")
print(table(n_assign))
if (mean(n_assign == 1) < 0.5)
  cat("\n  NOTE: most cluster-2 cells map to >1 lineage, so the lineage-restricted\n",
      "        progenitor sets overlap substantially. Analyses 3 and 4 will largely\n",
      "        mirror 1 and 2 rather than being independent comparisons.\n")

cat(sprintf("\nL3 n L4 overlap within cluster 2: %d cells\n",
            length(intersect(intersect(C2, L$Lineage3), intersect(C2, L$Lineage4)))))

sets <- list(
  "1_C2all_vs_C6" = list(ref = C2,                      target = C6),
  "2_C2all_vs_C4" = list(ref = C2,                      target = C4),
  "3_C2L3_vs_C6"  = list(ref = intersect(C2, L$Lineage3),
                         target = intersect(C6, L$Lineage3)),
  "4_C2L4_vs_C4"  = list(ref = intersect(C2, L$Lineage4),
                         target = intersect(C4, L$Lineage4))
)

LABELS <- c("1_C2all_vs_C6" = "1. All cluster 2\nvs cluster 6",
            "2_C2all_vs_C4" = "2. All cluster 2\nvs cluster 4",
            "3_C2L3_vs_C6"  = "3. L3 cluster 2\nvs L3 cluster 6",
            "4_C2L4_vs_C4"  = "4. L4 cluster 2\nvs L4 cluster 4")

cat("\n=========== COMPARISON GROUP SIZES ===========\n")
for (nm in names(sets))
  cat(sprintf("  %-16s  progenitor n=%4d   endpoint n=%4d\n",
              nm, length(sets[[nm]]$ref), length(sets[[nm]]$target)))

run_mast <- function(nm, s) {
  cat(sprintf("\n=== %s ===\n", nm))
  if (length(s$ref) < 3 || length(s$target) < 3) {
    cat("  too few cells, skipping\n"); return(NULL)
  }

  grp <- setNames(rep(NA_character_, ncol(seu)), colnames(seu))
  grp[s$ref]     <- "ref"
  grp[s$target]  <- "target"
  seu$mast_group <- grp

  d <- FindMarkers(seu,
                   ident.1 = "target", ident.2 = "ref",
                   group.by = "mast_group",
                   test.use = "MAST",
                   latent.vars = "nFeature_RNA",
                   min.pct = MIN_PCT,
                   logfc.threshold = TEST_LFC)

  d$FDR <- p.adjust(d$p_val, method = "BH")
  d <- d[order(-abs(d$avg_log2FC)), ]
  write.csv(d, file.path(OUTPUT_DIR, paste0("MAST_", nm, ".csv")))

  sig <- subset(d, FDR < FDR_CUT & abs(avg_log2FC) > LFC_CUT)
  cat(sprintf("  %d genes tested | %d significant (FDR<%.2f, |log2FC|>%.1f)\n",
              nrow(d), nrow(sig), FDR_CUT, LFC_CUT))
  cat(sprintf("  up in endpoint: %d    down: %d\n",
              sum(sig$avg_log2FC > 0), sum(sig$avg_log2FC < 0)))
  print(head(sig[, c("avg_log2FC","pct.1","pct.2","p_val","FDR")], 10))
  invisible(d)
}

mast <- setNames(lapply(names(sets), function(nm) run_mast(nm, sets[[nm]])),
                 names(sets))

sig_genes <- function(nm) {
  d <- mast[[nm]]
  if (is.null(d)) return(character(0))
  rownames(subset(d, FDR < FDR_CUT & abs(avg_log2FC) > LFC_CUT))
}

overlap_summary <- do.call(rbind, lapply(
  list(c("1_C2all_vs_C6", "3_C2L3_vs_C6", "L3_cluster6"),
       c("2_C2all_vs_C4", "4_C2L4_vs_C4", "L4_cluster4")),
  function(x) {
    a <- sig_genes(x[1]); b <- sig_genes(x[2])
    cat(sprintf("\n%s\n  pooled cluster 2      : %d genes\n  lineage-restricted    : %d genes\n  shared                : %d\n  unique to restricted  : %d\n  unique to pooled      : %d\n",
                x[3], length(a), length(b), length(intersect(a,b)),
                length(setdiff(b,a)), length(setdiff(a,b))))
    if (length(setdiff(b, a)))
      write.csv(data.frame(gene = setdiff(b, a)),
                file.path(OUTPUT_DIR, paste0("MAST_unique_to_restricted_", x[3], ".csv")),
                row.names = FALSE)
    data.frame(comparison = x[3], pooled = length(a), restricted = length(b),
               shared = length(intersect(a,b)),
               only_restricted = length(setdiff(b,a)),
               only_pooled = length(setdiff(a,b)))
  }))
write.csv(overlap_summary,
          file.path(OUTPUT_DIR, "MAST_pooled_vs_restricted_summary.csv"),
          row.names = FALSE)

umap_mat <- Embeddings(seu, names(seu@reductions)[grep("umap|UMAP",
                       names(seu@reductions))[1]])
plot_df  <- data.frame(UMAP1 = umap_mat[, 1], UMAP2 = umap_mat[, 2],
                       cluster = cl[rownames(umap_mat)],
                       row.names = rownames(umap_mat))
centroids <- plot_df %>% group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

hl <- do.call(rbind, lapply(names(L), function(ln) {
  cells <- intersect(C2, L[[ln]])
  if (!length(cells)) return(NULL)
  data.frame(plot_df[cells, c("UMAP1","UMAP2")], lineage = ln, row.names = NULL)
}))
lab <- hl %>% count(lineage) %>% mutate(txt = sprintf("n = %d", n))

pA <- ggplot() +
  geom_point(data = plot_df, aes(UMAP1, UMAP2), colour = "grey88", size = 0.3) +
  geom_point(data = hl, aes(UMAP1, UMAP2, colour = lineage),
             size = 0.8, alpha = 0.85) +
  geom_text(data = centroids, aes(UMAP1, UMAP2, label = cluster),
            fontface = "bold", size = 3.2, colour = "grey25") +
  geom_text(data = lab, aes(x = Inf, y = -Inf, label = txt),
            hjust = 1.15, vjust = -0.8, size = 3.6, fontface = "bold") +
  facet_wrap(~lineage, ncol = 2) +
  scale_colour_manual(values = lin_cb, guide = "none") +
  labs(title = "Cluster 2 progenitor cells by lineage assignment",
       subtitle = "Coloured = cluster 2 cells on that lineage (Slingshot weight > 0.5)",
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_blank())

sizes <- do.call(rbind, lapply(names(sets), function(nm) data.frame(
  comparison = LABELS[nm],
  group = factor(c("Progenitor (cluster 2)", "Endpoint"),
                 levels = c("Progenitor (cluster 2)", "Endpoint")),
  n = c(length(sets[[nm]]$ref), length(sets[[nm]]$target))
)))

pB <- ggplot(sizes, aes(x = comparison, y = n, fill = group)) +
  geom_col(position = position_dodge(0.75), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.75),
            vjust = -0.4, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = c("Progenitor (cluster 2)" = "#0072B2",
                               "Endpoint" = "#D55E00"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Cells per group in each MAST comparison",
       x = NULL, y = "Number of cells") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey92"))

deg_counts <- do.call(rbind, lapply(names(sets), function(nm) {
  d <- mast[[nm]]
  if (is.null(d)) return(NULL)
  s <- subset(d, FDR < FDR_CUT & abs(avg_log2FC) > LFC_CUT)
  data.frame(comparison = LABELS[nm],
             direction = factor(c("Up in endpoint", "Down in endpoint"),
                                levels = c("Up in endpoint", "Down in endpoint")),
             n = c(sum(s$avg_log2FC > 0), sum(s$avg_log2FC < 0)))
}))

pC <- ggplot(deg_counts, aes(x = comparison, y = n, fill = direction)) +
  geom_col(position = position_dodge(0.75), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.75),
            vjust = -0.4, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = c("Up in endpoint"   = "#009E73",
                               "Down in endpoint" = "#CC79A7"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Differentially expressed genes per comparison",
       subtitle = sprintf("MAST, BH-adjusted FDR < %.2f and |log2FC| > %.1f",
                          FDR_CUT, LFC_CUT),
       x = NULL, y = "Number of genes") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey92"))

cat("\nExporting figures...\n")
export_plot(pA, "MAST_cluster2_lineage_split",   10, 9)
export_plot(pB, "MAST_comparison_cell_counts",    9, 6)
export_plot(pC, "MAST_DEG_counts",                9, 6)

print(pA); print(pB); print(pC)

cat("\n=========== DONE ===========\n")
cat("Result tables:\n")
for (nm in names(sets)) cat(sprintf("  MAST_%s.csv\n", nm))
cat("  MAST_pooled_vs_restricted_summary.csv\n")
cat("Figures: MAST_cluster2_lineage_split, MAST_comparison_cell_counts, MAST_DEG_counts\n")
cat("         (each as .svg / .pdf / .png@600 / .tiff@600)\n")
cat("\nUse the FDR column, not p_val_adj (that one is Bonferroni).\n")
