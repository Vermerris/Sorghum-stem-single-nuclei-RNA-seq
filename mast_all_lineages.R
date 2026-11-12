if (!requireNamespace("MAST", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("MAST", ask = FALSE)
}
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")

library(Seurat); library(SingleCellExperiment); library(MAST)
library(ggplot2); library(dplyr); library(svglite)

OUTPUT_DIR <- "."
dir.create(OUTPUT_DIR, showWarnings = FALSE)

DPI           <- 600
FDR_CUT       <- 0.05
LFC_CUT       <- 0.5
MIN_PCT       <- 0.1
TEST_LFC      <- 0.25
SEEDS         <- 1:3
SKIP_EXISTING <- TRUE
MIN_GROUP_N   <- 50
WARN_GROUP_N  <- 200

GROUPS <- list(
  L1 = list(label = "L1",    target = "1", restrict = 1,      colour = "#0072B2"),
  L2 = list(label = "L2",    target = "5", restrict = 2,      colour = "#D55E00"),
  L3 = list(label = "L3",    target = "6", restrict = 3,      colour = "#009E73"),
  L4 = list(label = "L4",    target = "4", restrict = 4,      colour = "#CC79A7"),
  C0 = list(label = "L1+L2", target = "0", restrict = c(1, 2), colour = "#56B4E9")
)

sig <- function(d) rownames(subset(d, FDR < FDR_CUT & abs(avg_log2FC) > LFC_CUT))
disp <- function(g) sprintf("%s → C%s", GROUPS[[g]]$label, GROUPS[[g]]$target)

export_plot <- function(p, base, w, h) {
  f <- file.path(OUTPUT_DIR, base)
  ggsave(paste0(f, ".svg"),  p, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(f, ".pdf"),  p, width = w, height = h)
  ggsave(paste0(f, ".png"),  p, width = w, height = h, dpi = DPI, bg = "white")
  ggsave(paste0(f, ".tiff"), p, width = w, height = h, dpi = DPI,
         bg = "white", compression = "lzw")
  cat(sprintf("  saved: %s\n", base))
}

cat("Loading fitted tradeSeq SCE...\n")
sce <- readRDS("SCE_fitted_Sorghum_Xy_Phl_TA_SubsetTradeseq.rds")
cat("Loading Seurat subset...\n")
seu <- readRDS("SC_Sorghum_Xy_Phl_TA_Subset.rds")
stopifnot(setequal(colnames(sce), colnames(seu)),
          identical(metadata(sce)$sorghum_analysis, seu@misc$sorghum_analysis))
DefaultAssay(seu) <- "RNA"
if (inherits(seu[["RNA"]], "Assay5") && length(Layers(seu[["RNA"]], search = "data")) > 1) {
  cat("Joining Seurat v5 layers...\n"); seu <- JoinLayers(seu)
}

cd  <- as.data.frame(colData(sce))
wts <- as.matrix(cd[, grep("^crv\\.cellWeights\\.", names(cd)), drop = FALSE])
stopifnot(identical(sub("^crv\\.cellWeights\\.", "", colnames(wts)),
                    paste0("Lineage", 1:4)))
shared <- intersect(colnames(seu), rownames(wts))
cl <- as.character(seu$seurat_clusters); names(cl) <- colnames(seu)

on_lineage <- function(li) {
  w <- setNames(rep(0, ncol(seu)), colnames(seu))
  w[shared] <- wts[shared, li]
  names(w)[w > 0.5]
}
L  <- setNames(lapply(1:ncol(wts), on_lineage), paste0("Lineage", 1:ncol(wts)))
C2 <- names(cl)[cl == "2"]

lineage_cells <- function(idx) Reduce(union, lapply(idx, function(i) L[[i]]))

in_cluster <- function(k) names(cl)[cl == k]

cat("\n=========== CELL COUNTS ===========\n")
cat(sprintf("Cluster 2 total cells: %d\n\n", length(C2)))

comp <- do.call(rbind, lapply(names(GROUPS), function(g) {
  cfg <- GROUPS[[g]]
  lc  <- lineage_cells(cfg$restrict)
  tgt <- in_cluster(cfg$target)
  data.frame(group = g, restriction = cfg$label, target_cluster = cfg$target,
             c2_all = length(C2), c2_restricted = length(intersect(C2, lc)),
             target_all = length(tgt), target_restricted = length(intersect(tgt, lc)),
             stringsAsFactors = FALSE)
}))
print(comp, row.names = FALSE)
write.csv(comp, file.path(OUTPUT_DIR, "cell_counts_per_comparison.csv"), row.names = FALSE)

small <- comp$group[comp$c2_restricted < WARN_GROUP_N]
if (length(small))
  cat(sprintf("\n  WARNING: %s %s fewer than %d restricted cluster 2 cells.\n  %s\n",
              paste(small, collapse = ", "), ifelse(length(small) == 1, "has", "have"),
              WARN_GROUP_N,
              "Restricted comparisons will be underpowered and their controls noisy."))

cat("\nLineages assigned per cluster-2 cell:\n")
print(table(rowSums(sapply(L, function(x) C2 %in% x))))

sets <- list()
for (g in names(GROUPS)) {
  cfg <- GROUPS[[g]]
  lc  <- lineage_cells(cfg$restrict)
  tgt <- in_cluster(cfg$target)

  sets[[paste0("pooled_",     g)]] <- list(
    ref = C2, target = tgt, group = g, kind = "pooled")
  sets[[paste0("restricted_", g)]] <- list(
    ref = intersect(C2, lc), target = intersect(tgt, lc), group = g, kind = "restricted")
}

cat("\n=========== COMPARISONS ===========\n")
for (nm in names(sets))
  cat(sprintf("  %-18s %-14s  ref n=%4d   target n=%4d\n",
              nm, disp(sets[[nm]]$group),
              length(sets[[nm]]$ref), length(sets[[nm]]$target)))

run_mast <- function(nm, s) {
  f <- file.path(OUTPUT_DIR, paste0("MAST_", nm, ".csv"))
  if (SKIP_EXISTING && file.exists(f)) {
    cat(sprintf("\n=== %s (%s) === cached\n", nm, disp(s$group)))
    return(read.csv(f, row.names = 1))
  }
  cat(sprintf("\n=== %s (%s) ===\n", nm, disp(s$group)))
  if (length(s$ref) < MIN_GROUP_N || length(s$target) < MIN_GROUP_N) {
    cat(sprintf("  group too small (ref %d, target %d), skipping\n",
                length(s$ref), length(s$target)))
    return(NULL)
  }

  grp <- setNames(rep(NA_character_, ncol(seu)), colnames(seu))
  grp[s$ref]     <- "ref"
  grp[s$target]  <- "target"
  seu$mast_group <- grp

  d <- FindMarkers(seu, ident.1 = "target", ident.2 = "ref",
                   group.by = "mast_group", test.use = "MAST",
                   latent.vars = "nFeature_RNA",
                   min.pct = MIN_PCT, logfc.threshold = TEST_LFC)
  d$FDR <- p.adjust(d$p_val, method = "BH")
  d <- d[order(-abs(d$avg_log2FC)), ]
  write.csv(d, f)

  s2 <- subset(d, FDR < FDR_CUT & abs(avg_log2FC) > LFC_CUT)
  cat(sprintf("  %d tested | %d significant | %d up, %d down\n",
              nrow(d), nrow(s2), sum(s2$avg_log2FC > 0), sum(s2$avg_log2FC < 0)))
  d
}

mast <- setNames(lapply(names(sets), function(nm) run_mast(nm, sets[[nm]])),
                 names(sets))

overlap <- do.call(rbind, lapply(names(GROUPS), function(g) {
  mp <- mast[[paste0("pooled_", g)]]
  mr <- mast[[paste0("restricted_", g)]]
  if (is.null(mp) || is.null(mr)) return(NULL)
  cfg <- GROUPS[[g]]

  a <- sig(mp); b <- sig(mr)
  only <- setdiff(b, a)

  if (length(only)) {
    reason <- vapply(only, function(gene) {
      if (!(gene %in% rownames(mp)))                   "not tested in pooled (below prefilter)"
      else if (mp[gene, "FDR"] >= FDR_CUT)             "FDR not significant in pooled"
      else if (abs(mp[gene, "avg_log2FC"]) <= LFC_CUT) "log2FC below threshold in pooled"
      else                                             "other"
    }, character(1))

    out <- data.frame(
      restricted_log2FC = mr[only, "avg_log2FC"],
      restricted_FDR    = mr[only, "FDR"],
      restricted_pct.1  = mr[only, "pct.1"],
      restricted_pct.2  = mr[only, "pct.2"],
      pooled_log2FC     = mp[only, "avg_log2FC"],
      pooled_FDR        = mp[only, "FDR"],
      reason_unique     = reason,
      row.names = only
    )
    out$lfc_gain <- abs(out$restricted_log2FC) - abs(out$pooled_log2FC)
    out <- out[order(out$restricted_FDR, -abs(out$restricted_log2FC)), ]
    write.csv(out, file.path(OUTPUT_DIR,
              sprintf("unique_to_restricted_%s_cluster%s.csv", g, cfg$target)))
    cat(sprintf("\n%s : %d unique-to-restricted genes\n", disp(g), nrow(out)))
    print(table(out$reason_unique))
  }

  data.frame(group = g, restriction = cfg$label, target_cluster = cfg$target,
             pooled = length(a), restricted = length(b),
             shared = length(intersect(a, b)),
             only_restricted = length(only), only_pooled = length(setdiff(a, b)),
             stringsAsFactors = FALSE)
}))
cat("\n=========== POOLED vs RESTRICTED ===========\n")
print(overlap, row.names = FALSE)
write.csv(overlap, file.path(OUTPUT_DIR, "pooled_vs_restricted_summary.csv"), row.names = FALSE)

DS_FILE <- file.path(OUTPUT_DIR, "downsample_control.csv")

if (SKIP_EXISTING && file.exists(DS_FILE)) {
  downsample <- read.csv(DS_FILE, stringsAsFactors = FALSE)
  cat("\nDownsample control: loaded from cache\n")
} else {
  downsample <- do.call(rbind, lapply(names(GROUPS), function(g) {
    cfg <- GROUPS[[g]]
    lc  <- lineage_cells(cfg$restrict)
    tgt <- intersect(in_cluster(cfg$target), lc)
    n_ref <- length(intersect(C2, lc))
    if (n_ref < MIN_GROUP_N || length(tgt) < MIN_GROUP_N) return(NULL)

    cat(sprintf("\nDownsample %s: %d random cluster 2 cells x %d seeds\n",
                disp(g), n_ref, length(SEEDS)))
    counts <- sapply(SEEDS, function(sd) {
      set.seed(sd)
      grp <- setNames(rep(NA_character_, ncol(seu)), colnames(seu))
      grp[sample(C2, n_ref)] <- "ref"
      grp[tgt]               <- "target"
      seu$mast_group <- grp
      d <- FindMarkers(seu, ident.1 = "target", ident.2 = "ref",
                       group.by = "mast_group", test.use = "MAST",
                       latent.vars = "nFeature_RNA",
                       min.pct = MIN_PCT, logfc.threshold = TEST_LFC)
      d$FDR <- p.adjust(d$p_val, method = "BH")
      length(sig(d))
    })
    cat(sprintf("  -> %s\n", paste(counts, collapse = ", ")))

    data.frame(group = g, n_ref = n_ref,
               random_mean = mean(counts), random_sd = sd(counts),
               random_values = paste(counts, collapse = "/"),
               stringsAsFactors = FALSE)
  }))
  write.csv(downsample, DS_FILE, row.names = FALSE)
}

decomp <- merge(overlap[, c("group","restriction","target_cluster","pooled","restricted")],
                downsample, by = "group")
decomp$loss_total       <- decomp$pooled      - decomp$restricted
decomp$loss_sample_size <- decomp$pooled      - decomp$random_mean
decomp$loss_composition <- decomp$random_mean - decomp$restricted
decomp$pct_composition  <- round(100 * decomp$loss_composition / decomp$loss_total, 1)

cat("\n=========== DECOMPOSITION ===========\n")
print(decomp, row.names = FALSE)
write.csv(decomp, file.path(OUTPUT_DIR, "decomposition_summary.csv"), row.names = FALSE)

umap_mat <- Embeddings(seu, names(seu@reductions)[grep("umap|UMAP",
                       names(seu@reductions))[1]])
plot_df <- data.frame(UMAP1 = umap_mat[, 1], UMAP2 = umap_mat[, 2],
                      cluster = cl[rownames(umap_mat)],
                      row.names = rownames(umap_mat))
centroids <- plot_df %>% group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

lin_cb <- c(Lineage1 = "#0072B2", Lineage2 = "#D55E00",
            Lineage3 = "#009E73", Lineage4 = "#CC79A7")

hl <- do.call(rbind, lapply(names(L), function(lin) {
  cells <- intersect(C2, L[[lin]])
  if (!length(cells)) return(NULL)
  data.frame(plot_df[cells, c("UMAP1","UMAP2")], lineage = lin, row.names = NULL)
}))
lab <- hl %>% count(lineage) %>% mutate(txt = sprintf("n = %d", n))

pA <- ggplot() +
  geom_point(data = plot_df, aes(UMAP1, UMAP2), colour = "grey88", size = 0.3) +
  geom_point(data = hl, aes(UMAP1, UMAP2, colour = lineage), size = 0.8, alpha = 0.85) +
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

sizes <- do.call(rbind, lapply(names(sets), function(nm) {
  s <- sets[[nm]]
  data.frame(label = sprintf("%s\n(%s)", disp(s$group), s$kind),
             group = factor(c("Progenitor (cluster 2)", "Target cluster"),
                            levels = c("Progenitor (cluster 2)", "Target cluster")),
             n = c(length(s$ref), length(s$target)))
}))

pB <- ggplot(sizes, aes(x = label, y = n, fill = group)) +
  geom_col(position = position_dodge(0.75), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.75),
            vjust = -0.4, size = 2.7, fontface = "bold") +
  scale_fill_manual(values = c("Progenitor (cluster 2)" = "#0072B2",
                               "Target cluster" = "#D55E00"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Cells per group in each comparison", x = NULL, y = "Number of cells") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        axis.text.x = element_text(size = 7.5),
        panel.grid.major.y = element_line(colour = "grey92"))

deg_counts <- do.call(rbind, lapply(names(sets), function(nm) {
  d <- mast[[nm]]; if (is.null(d)) return(NULL)
  s  <- sets[[nm]]
  s2 <- subset(d, FDR < FDR_CUT & abs(avg_log2FC) > LFC_CUT)
  data.frame(label = sprintf("%s\n(%s)", disp(s$group), s$kind),
             direction = factor(c("Up in target", "Down in target"),
                                levels = c("Up in target", "Down in target")),
             n = c(sum(s2$avg_log2FC > 0), sum(s2$avg_log2FC < 0)))
}))

pC <- ggplot(deg_counts, aes(x = label, y = n, fill = direction)) +
  geom_col(position = position_dodge(0.75), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.75),
            vjust = -0.4, size = 2.7, fontface = "bold") +
  scale_fill_manual(values = c("Up in target"   = "#009E73",
                               "Down in target" = "#CC79A7"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Differentially expressed genes per comparison",
       subtitle = sprintf("MAST, BH FDR < %.2f and |log2FC| > %.1f", FDR_CUT, LFC_CUT),
       x = NULL, y = "Number of genes") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        axis.text.x = element_text(size = 7.5),
        panel.grid.major.y = element_line(colour = "grey92"))

dd <- do.call(rbind, lapply(seq_len(nrow(decomp)), function(i) {
  r <- decomp[i, ]
  data.frame(
    panel = sprintf("%s → C%s   (n = %d)", r$restriction, r$target_cluster, r$n_ref),
    type = factor(c("All cluster 2", "Random subsample", "Lineage-restricted"),
                  levels = c("All cluster 2", "Random subsample", "Lineage-restricted")),
    n = c(r$pooled, r$random_mean, r$restricted))
}))

pD <- ggplot(dd, aes(x = type, y = n, fill = type)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = round(n)), vjust = -0.45, size = 3.1, fontface = "bold") +
  facet_wrap(~panel, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("All cluster 2"      = "#9CA3AF",
                               "Random subsample"   = "#0072B2",
                               "Lineage-restricted" = "#009E73"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Size-matched downsample control",
       subtitle = "Random vs lineage-restricted progenitor reference at identical n — the gap between them is compositional",
       x = NULL, y = "Significant genes") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"))

cat("\nExporting figures...\n")
export_plot(pA, "cluster2_lineage_split",    10, 9)
export_plot(pB, "comparison_cell_counts",    14, 6)
export_plot(pC, "DEG_counts",                14, 6)
export_plot(pD, "downsample_decomposition",  14, 8)

print(pA); print(pB); print(pC); print(pD)

cat("\n=========== DONE ===========\n")
cat("Output in:", OUTPUT_DIR, "\n")
cat("Use the FDR column (BH), not p_val_adj (Bonferroni).\n")
