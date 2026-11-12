library(tradeSeq)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(gridExtra)
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
library(svglite)

OUTPUT_DIR <- "."
dir.create(OUTPUT_DIR, showWarnings = FALSE)

DPI         <- 600
VIRIDIS_OPT <- "plasma"

lin_cb <- c(Lineage1 = "#0072B2",
            Lineage2 = "#D55E00",
            Lineage3 = "#009E73",
            Lineage4 = "#CC79A7")

lin_lty <- c(Lineage1 = "solid", Lineage2 = "22",
             Lineage3 = "44",    Lineage4 = "1343")

export_figure <- function(draw_fn, basename, w, h) {
  p <- file.path(OUTPUT_DIR, paste0(basename, c(".svg", ".pdf", ".png", ".tiff")))
  svglite::svglite(p[1], width = w, height = h);                         draw_fn(); dev.off()
  pdf(p[2], width = w, height = h);                                      draw_fn(); dev.off()
  png(p[3], width = w, height = h, units = "in", res = DPI, bg = "white"); draw_fn(); dev.off()
  tiff(p[4], width = w, height = h, units = "in", res = DPI,
       compression = "lzw", bg = "white");                               draw_fn(); dev.off()
  cat(sprintf("  saved: %s  (svg / pdf / png@%d / tiff@%d)\n", basename, DPI, DPI))
}

cat("Loading fitted tradeSeq SCE...\n")
sce <- readRDS("SCE_fitted_Sorghum_Xy_Phl_TA_SubsetTradeseq.rds")

cat("Loading original Seurat object...\n")
seu <- readRDS("SC_Sorghum_Xy_Phl_TA_Subset.rds")
stopifnot(setequal(colnames(sce), colnames(seu)),
          identical(metadata(sce)$sorghum_analysis, seu@misc$sorghum_analysis))

cd      <- as.data.frame(colData(sce))
pt_cols <- grep("^crv\\.pseudotime\\.",  names(cd), value = TRUE)
wt_cols <- grep("^crv\\.cellWeights\\.", names(cd), value = TRUE)

pt  <- as.matrix(cd[, pt_cols, drop = FALSE])
wts <- as.matrix(cd[, wt_cols, drop = FALSE])
colnames(pt)  <- sub("^crv\\.pseudotime\\.",  "", colnames(pt))
colnames(wts) <- sub("^crv\\.cellWeights\\.", "", colnames(wts))
stopifnot(identical(colnames(pt), paste0("Lineage", 1:4)),
          identical(colnames(wts), paste0("Lineage", 1:4)))

cat("Pseudotime columns:", colnames(pt), "\n")
cat("Cells per lineage (weight > 0.5):", colSums(wts > 0.5), "\n")

n_k      <- nknots(sce)
knot_raw <- metadata(sce)$tradeSeq$knots
knot_vec <- sort(unique(as.numeric(unlist(knot_raw))))
cat(sprintf("\n%d knots found:\n", length(knot_vec)))
for (i in seq_along(knot_vec)) cat(sprintf("  K%d : %.3f\n", i, knot_vec[i]))

possible_cols <- c("seurat_clusters","cluster","Cluster","cluster_id",
                   "louvain","leiden","label","Label","celltype","RNA_snn_res.0.5")
CLUSTER_COL <- possible_cols[possible_cols %in% names(seu@meta.data)][1]
cat(sprintf("\nUsing cluster column: '%s'\n", CLUSTER_COL))
cat("Cluster levels:", sort(unique(as.character(seu@meta.data[[CLUSTER_COL]]))), "\n")

umap_name <- names(seu@reductions)[grep("umap|UMAP", names(seu@reductions))[1]]
umap_mat  <- Embeddings(seu, umap_name)

shared <- intersect(colnames(sce), rownames(umap_mat))
cat(sprintf("Shared barcodes: %d / %d SCE cells\n", length(shared), ncol(sce)))
if (length(shared) == 0) stop("No shared barcodes!")

umap_df <- data.frame(
  UMAP1   = umap_mat[shared, 1],
  UMAP2   = umap_mat[shared, 2],
  cluster = as.character(seu@meta.data[shared, CLUSTER_COL]),
  row.names = shared
)

pt_s  <- pt[shared,  , drop = FALSE]
wts_s <- wts[shared, , drop = FALSE]

knot_coords <- do.call(rbind, lapply(seq_len(ncol(pt_s)), function(li) {
  on_lin <- which(wts_s[, li] > 0.5)
  if (!length(on_lin)) return(NULL)
  do.call(rbind, lapply(seq_along(knot_vec), function(ki) {
    cidx <- on_lin[which.min(abs(pt_s[on_lin, li] - knot_vec[ki]))]
    data.frame(knot=ki, label=paste0("K",ki),
               lineage=colnames(pt_s)[li], lin_idx=li,
               pseudotime=knot_vec[ki],
               UMAP1=umap_df$UMAP1[cidx],
               UMAP2=umap_df$UMAP2[cidx],
               stringsAsFactors=FALSE)
  }))
}))

for (li in seq_len(ncol(pt_s))) {
  lin      <- colnames(pt_s)[li]
  on       <- which(wts_s[, li] > 0.5)
  bins     <- cut(pt_s[on, li], breaks = 10)
  dominant <- tapply(umap_df$cluster[on], bins,
                     function(x) names(sort(table(x), decreasing=TRUE)[1]))
  transitions <- rle(as.character(dominant[!is.na(dominant)]))$values
  cat(sprintf("\n%s: %s\n", lin, paste(transitions, collapse = " -> ")))
}

lineage_counts <- do.call(rbind, lapply(seq_len(ncol(pt_s)), function(li) {
  lin   <- colnames(pt_s)[li]
  cells <- rownames(wts_s)[wts_s[, li] > 0.5]
  ct    <- table(umap_df[cells, "cluster"])
  data.frame(lineage=lin, cluster=names(ct), n_cells=as.integer(ct),
             stringsAsFactors=FALSE)
}))
cat("\nCell counts per cluster per lineage:\n")
print(lineage_counts)

p_counts <- ggplot(lineage_counts, aes(x = cluster, y = n_cells, fill = lineage)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n_cells), hjust = -0.15, size = 3.5, fontface = "bold") +
  coord_flip(clip = "off") +
  facet_wrap(~lineage, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = lin_cb) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Cells per Cluster per Lineage", x = "Cluster", y = "Number of Cells") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.major.x = element_line(colour = "grey92"))

centroids <- umap_df %>%
  group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

p_cb <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(colour = "grey85", size = 0.3, alpha = 0.5) +
  geom_text(data = centroids, aes(UMAP1, UMAP2, label = cluster),
            fontface = "bold", size = 4.5, colour = "grey15", inherit.aes = FALSE) +
  geom_path(data = knot_coords %>% arrange(lin_idx, knot),
            aes(UMAP1, UMAP2, group = lineage, colour = lineage, linetype = lineage),
            linewidth = 1.2, inherit.aes = FALSE) +
  geom_point(data = knot_coords, aes(UMAP1, UMAP2, fill = lineage),
             shape = 21, size = 3.5, colour = "white", stroke = 0.8,
             inherit.aes = FALSE) +
  geom_label_repel(data = knot_coords,
                   aes(UMAP1, UMAP2, label = label, fill = lineage),
                   colour = "white", fontface = "bold", size = 3,
                   box.padding = 0.35, point.padding = 0.4,
                   max.overlaps = Inf, show.legend = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(values   = lin_cb,  name = "Lineage") +
  scale_fill_manual(values     = lin_cb,  name = "Lineage") +
  scale_linetype_manual(values = lin_lty, name = "Lineage") +
  labs(title = "TradeSeq Knot Positions per Lineage", x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.key.width = unit(1.4, "cm"))

pt_plots_cb <- lapply(seq_len(ncol(pt_s)), function(li) {
  lin <- colnames(pt_s)[li]
  df  <- umap_df
  df$pt <- pt_s[rownames(umap_df), li]
  df$on <- wts_s[rownames(umap_df), li] > 0.5
  kdf <- subset(knot_coords, lin_idx == li) %>% arrange(knot)

  ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(data = subset(df, !on), colour = "grey88", size = 0.25, alpha = 0.4) +
    geom_point(data = subset(df,  on), aes(colour = pt),  size = 0.5,  alpha = 0.85) +
    scale_colour_viridis_c(option = VIRIDIS_OPT, name = "Pseudotime") +
    geom_text(data = centroids, aes(UMAP1, UMAP2, label = cluster),
              fontface = "bold", size = 3.5, colour = "grey15", inherit.aes = FALSE) +
    geom_path(data = kdf, aes(UMAP1, UMAP2), colour = "white",
              linewidth = 2.4, lineend = "round", inherit.aes = FALSE) +
    geom_path(data = kdf, aes(UMAP1, UMAP2), colour = "black",
              linewidth = 1.1, lineend = "round", inherit.aes = FALSE) +
    geom_point(data = kdf, aes(UMAP1, UMAP2), shape = 21, size = 3.4,
               fill = "white", colour = "black", stroke = 0.9, inherit.aes = FALSE) +
    geom_label_repel(data = kdf, aes(UMAP1, UMAP2, label = label),
                     fill = lin_cb[lin], colour = "white", fontface = "bold",
                     size = 3.4, box.padding = 0.4, max.overlaps = Inf,
                     label.size = 0.3, segment.colour = "black",
                     inherit.aes = FALSE, show.legend = FALSE) +
    labs(title = lin, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 13, colour = lin_cb[lin]))
})

cat("\nExporting figures...\n")
export_figure(function() print(p_cb),
              "UMAP_knots_overview_colorblind", 14, 9)
export_figure(function() grid.arrange(grobs = pt_plots_cb, ncol = 2),
              "UMAP_knots_per_lineage_colorblind", 14, 12)
export_figure(function() print(p_counts),
              "cells_per_cluster_per_lineage_colorblind", 10, 8)

cat("\nRunning earlyDETest K3->K4 (exit cluster 2)...\n")
res_exit    <- earlyDETest(sce, knots = c(3, 4), pairwise = FALSE)
res_exit_pw <- earlyDETest(sce, knots = c(3, 4), pairwise = TRUE)

cat("Running earlyDETest K2->K3 (entry into cluster 2)...\n")
res_entry    <- earlyDETest(sce, knots = c(2, 3), pairwise = FALSE)
res_entry_pw <- earlyDETest(sce, knots = c(2, 3), pairwise = TRUE)

add_fdr <- function(r) {
  orig <- names(r)
  for (pc in grep("^pvalue", orig, value = TRUE))
    r[[sub("^pvalue", "FDR", pc)]] <- p.adjust(r[[pc]], method = "BH")
  ord <- unlist(lapply(orig, function(nm)
    if (grepl("^pvalue", nm)) c(nm, sub("^pvalue", "FDR", nm)) else nm))
  r[, ord]
}

res_exit     <- add_fdr(res_exit)
res_exit_pw  <- add_fdr(res_exit_pw)
res_entry    <- add_fdr(res_entry)
res_entry_pw <- add_fdr(res_entry_pw)

cat("\n=========== SIGNIFICANCE SUMMARY ===========\n")
for (nm in c("res_entry", "res_exit")) {
  r <- get(nm)
  cat(sprintf("\n%s  (%d genes tested, %d NA)\n", nm, nrow(r), sum(is.na(r$pvalue))))
  cat(sprintf("   raw p < 0.05 : %5d\n", sum(r$pvalue < 0.05, na.rm = TRUE)))
  cat(sprintf("   FDR   < 0.05 : %5d   <-- use this\n", sum(r$FDR < 0.05, na.rm = TRUE)))
  cat(sprintf("   FDR   < 0.01 : %5d\n", sum(r$FDR < 0.01, na.rm = TRUE)))
}
cat("\n")

PAIRS <- list(L1 = c("1vs2","1vs3","1vs4"),
              L2 = c("1vs2","2vs3","2vs4"),
              L3 = c("1vs3","2vs3","3vs4"),
              L4 = c("1vs4","2vs4","3vs4"))
ENDPOINT <- c(L1 = "cluster 1", L2 = "cluster 5",
              L3 = "cluster 6", L4 = "cluster 4")

score_lineages <- function(pw, tag) {

  sc  <- sapply(PAIRS, function(p)
           apply(pw[, paste0("waldStat_", p), drop = FALSE], 1, min))

  fdr <- sapply(PAIRS, function(p)
           apply(pw[, paste0("FDR_", p), drop = FALSE], 1, max))

  ok  <- complete.cases(sc) & complete.cases(fdr) & !is.na(pw$FDR)
  sc  <- sc[ok, , drop = FALSE]; fdr <- fdr[ok, , drop = FALSE]
  gl  <- pw$FDR[ok]; gn <- rownames(pw)[ok]

  best <- max.col(sc, ties.method = "first")
  idx  <- cbind(seq_len(nrow(sc)), best)

  out <- data.frame(
    outlier_lineage = names(PAIRS)[best],
    endpoint        = unname(ENDPOINT[best]),
    min_waldStat    = sc[idx],
    worst_pair_FDR  = fdr[idx],
    global_FDR      = gl,
    row.names = gn, stringsAsFactors = FALSE
  )
  out <- out[order(-out$min_waldStat), ]

  hits <- subset(out, worst_pair_FDR < 0.05 & global_FDR < 0.05)
  write.csv(hits, file.path(OUTPUT_DIR, paste0("lineage_specific_", tag, ".csv")))

  cat(sprintf("\n=== %s ===\n%d genes lineage-specific (all 3 pairs FDR<0.05)\n",
              tag, nrow(hits)))
  print(table(hits$outlier_lineage, hits$endpoint))
  cat("\nTop 15:\n"); print(head(hits, 15))
  invisible(out)
}

score_exit  <- score_lineages(res_exit_pw,  "K3_K4_exit")
score_entry <- score_lineages(res_entry_pw, "K2_K3_entry")

annotate_pw <- function(pw, sc_out) {
  for (col in c("outlier_lineage","endpoint","min_waldStat","worst_pair_FDR"))
    pw[[col]] <- sc_out[rownames(pw), col]
  pw
}
res_exit_pw  <- annotate_pw(res_exit_pw,  score_exit)
res_entry_pw <- annotate_pw(res_entry_pw, score_entry)

write.csv(res_entry,    file.path(OUTPUT_DIR, "earlyDE_K2_K3_entry_overall.csv"))
write.csv(res_entry_pw, file.path(OUTPUT_DIR, "earlyDE_K2_K3_entry_pairwise.csv"))
write.csv(res_exit,     file.path(OUTPUT_DIR, "earlyDE_K3_K4_exit_overall.csv"))
write.csv(res_exit_pw,  file.path(OUTPUT_DIR, "earlyDE_K3_K4_exit_pairwise.csv"))
cat("\nResult tables written.\n")

top_exit  <- rownames(res_exit)[order(-res_exit$waldStat)][1:9]
top_entry <- rownames(res_entry)[order(-res_entry$waldStat)][1:9]

STAT_COLS <- c("waldStat","df","pvalue","FDR","fcMedian")
cat("\nTop 20 — EXIT cluster 2 (K3->K4), ranked by waldStat:\n")
print(head(res_exit[order(-res_exit$waldStat), STAT_COLS], 20))
cat("\nTop 20 — ENTRY into cluster 2 (K2->K3), ranked by waldStat:\n")
print(head(res_entry[order(-res_entry$waldStat), STAT_COLS], 20))

pdf(file.path(OUTPUT_DIR, "earlyDE_K3_K4_exit_smoothers.pdf"), width = 14, height = 10)
for (g in top_exit) tryCatch({
  p <- plotSmoothers(sce, counts(sce), gene = g, border = TRUE)
  print(p + ggtitle(g))
}, error = function(e) message("skip ", g, ": ", e$message))
dev.off()

pdf(file.path(OUTPUT_DIR, "earlyDE_K2_K3_entry_smoothers.pdf"), width = 14, height = 10)
for (g in top_entry) tryCatch({
  p <- plotSmoothers(sce, counts(sce), gene = g, border = TRUE)
  print(p + ggtitle(g))
}, error = function(e) message("skip ", g, ": ", e$message))
dev.off()
cat("Saved smoother PDFs.\n")

GENE <- "Sobic.005G029700.v3.2"

p_gene <- plotSmoothers(sce, counts(sce), gene = GENE, border = TRUE) + ggtitle(GENE)
print(p_gene)
ggsave(file.path(OUTPUT_DIR, paste0("smoother_", GENE, ".png")),
       p_gene, width = 8, height = 5, dpi = DPI, bg = "white")
ggsave(file.path(OUTPUT_DIR, paste0("smoother_", GENE, ".svg")),
       p_gene, width = 8, height = 5, device = svglite::svglite)

gene_tables <- list(
  "K2->K3 entry (overall)"  = res_entry,
  "K2->K3 entry (pairwise)" = res_entry_pw,
  "K3->K4 exit  (overall)"  = res_exit,
  "K3->K4 exit  (pairwise)" = res_exit_pw
)
for (nm in names(gene_tables)) {
  r <- gene_tables[[nm]]
  cat("\n──", nm, "──\n")
  if (GENE %in% rownames(r)) print(t(r[GENE, , drop = FALSE]))
  else cat("  gene not found\n")
}

cat("\nDone. All output in:", OUTPUT_DIR, "\n")
