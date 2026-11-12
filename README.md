# Sorghum stem single-nuclei RNA-seq

Scripts accompanying **Identification of developmental programs governing cell identity in stems of sorghum (*Sorghum bicolor* (L.) Moench) using single-nuclei transcriptomics.**

Use Cell Ranger 7.0.0, gffread, and R with `Seurat` (version 5), `tradeSeq`, `slingshot`, `SingleCellExperiment`, `BiocParallel`, `Matrix`, `RColorBrewer`, `ggplot2`, `dplyr`, `viridis`, `ggrepel`, `gridExtra`, `svglite`, and `MAST`. The shell scripts use SLURM and environment modules. The trajectory script records the installed package versions in `Sorghum_trajectory_sessionInfo.txt`.

Keep the scripts and input files in one working directory. Start your terminal or R/RStudio session there; all paths are relative. Set the SLURM email and resource options for your cluster.

Required inputs are not included:

- FASTQ files for `N609_D29`, `N609_D31`, `Tx430_D29`, and `Tx430_D31`, together in the working directory. The `--sample` values must match their FASTQ sample prefixes.
- Reference FASTA `Sbicolor_454_v3.0.1.softmasked.fa` and annotation `Sbicolor_454_v3.1.1.gene_exons.gff3`.

Follow this order:

| Step | Script | Purpose |
| --- | --- | --- |
| 1 | `cellranger_to_build_sorghum_reference_genome.sh` | Convert the annotation to GTF and build the reference. |
| 2 | `cellranger_to_generating_matrix_counts.sh` | Generate count matrices for the four samples. |
| 3 | `scRNAseq_filtering_and_clustering.R` | Filter, integrate, and cluster nuclei; save `assay_RNA_50_Integrated_CellRanger_50PC_0.5res.rds`. |
| 4 | `trajectory_inference_tradeseqSorghumWholeDataset.R` | Subset and recluster vascular cells, save the Seurat subset, infer Slingshot trajectories, and fit tradeSeq models. |
| 5 | `tradeseq_knot_analysis.R` | Generate knot plots, early differential-expression tests, and lineage-specific gene tables. |
| 6 (optional) | `visualize_knots_earlyDE_FDR_corrected.R` | Generate additional plot variants and tables from the same trajectory outputs. |
| 7 | `mast_cluster2_endpoints.R` | Compare cluster 2 with endpoint clusters 6 and 4 using MAST. |
| 8 | `mast_all_lineages.R` | Run MAST comparisons for all four lineages and cluster 0, with subsampling controls. |

Submit steps 1–2 with `sbatch SCRIPT.sh`, waiting for step 1 to finish first. The reference is created in `Sbicolor_v3.1.1/`. Skip these steps if you already have filtered matrices. For step 3, place each sample's Cell Ranger `outs/filtered_feature_bc_matrix` directory in the working directory and name it `filtered_feature_bc_matrix_SAMPLE`, such as `filtered_feature_bc_matrix_N609_D29`.

Run the R scripts with `Rscript SCRIPT.R` or in R/RStudio, following the order above. For example:

```sh
Rscript scRNAseq_filtering_and_clustering.R
Rscript trajectory_inference_tradeseqSorghumWholeDataset.R
```

Step 4 performs subsetting, reclustering, and trajectory analysis together. It fits eight-knot models for four vascular lineages using four workers and writes `SC_Sorghum_Xy_Phl_TA_Subset.rds` and `SCE_fitted_Sorghum_Xy_Phl_TA_SubsetTradeseq.rds`. Steps 5–8 read this matching pair. The lineage mapping is saved in `Sorghum_vascular_lineages.csv`. Checks stop the workflow if the expected cell count, cluster labels, or lineage paths differ; verify those identities before changing the checks. Package versions can affect the results despite fixed random seeds.

R results are written directly to the working directory. Use a fresh working directory for each analysis: step 4 preserves existing output RDS files by stopping, and step 8 reuses existing result tables. Step 6 replaces some files from step 5; run it when those additional plot variants are needed. Steps 7 and 8 can run independently after step 4.

`Subset_lineages_tradeSeq.R` is a separate analysis starting from `assay_RNA_50_Integrated_CellRanger_50PC_0.5res.rds`, the output of step 3. It fits six-knot models for lineages 1–3 and saves `Sorghum_Tradeseq_lineages1to3.rds` in the working directory. This three-lineage fit is separate from the four-lineage analysis in steps 5–8.
