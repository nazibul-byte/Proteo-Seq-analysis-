#!/usr/bin/env Rscript

# PDO manuscript-specific GeneScore analysis

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =============================
# User configuration
# =============================
archr_project_dir <- "/path/to/PDO/ArchRProject_QCfiltered"
cluster_col <- "ClustersKmeans_K5"

proj <- loadArchRProject(path = archr_project_dir)

if (!cluster_col %in% colnames(getCellColData(proj))) {
  stop("Cluster column not found: ", cluster_col)
}

message("PDO cells in project: ", length(getCellNames(proj)))
print(table(getCellColData(proj, cluster_col)[, 1], useNA = "ifany"))

# =============================
# 1. Define manuscript gene panel
# =============================
genes_plot <- c(
  "ESR1", "PGR", "GATA3", "FOXA1",
  "EPCAM", "KRT8", "HSP90AA1", "KRT14", "KRT5", "TP63", "KRT17",
  "ITGA6", "KRT15", "EGFR", "LGALS7", "ACTA2", "TAGLN", "MYLK", "CNN1"
)
genes_plot <- unique(toupper(genes_plot))

# =============================
# 2. Extract per-cell GeneScores
# =============================
gene_score_se <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix"
)

archr_gene_names <- as.character(rowData(gene_score_se)$name)
gene_index <- match(toupper(genes_plot), toupper(archr_gene_names))

genes_found <- genes_plot[!is.na(gene_index)]
genes_missing <- genes_plot[is.na(gene_index)]
gene_index_found <- gene_index[!is.na(gene_index)]

message("Requested genes: ", length(genes_plot))
message("Found genes: ", length(genes_found))
if (length(genes_missing) > 0) {
  message("Missing genes: ", paste(genes_missing, collapse = ", "))
}
if (length(genes_found) == 0) stop("None of the selected genes were found.")

gene_score_mat <- assay(
  gene_score_se,
  "GeneScoreMatrix"
)[gene_index_found, , drop = FALSE]
rownames(gene_score_mat) <- genes_found

gene_score_df <- data.frame(
  Cell = colnames(gene_score_mat),
  t(as.matrix(gene_score_mat)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cell_meta <- as.data.frame(
  getCellColData(ArchRProj = proj, select = cluster_col)
)
cluster_df <- data.frame(
  Cell = rownames(cell_meta),
  Cluster = as.character(cell_meta[[cluster_col]]),
  stringsAsFactors = FALSE
)

gene_score_df <- gene_score_df %>%
  left_join(cluster_df, by = "Cell") %>%
  relocate(Cell, Cluster)

if (any(is.na(gene_score_df$Cluster))) {
  warning(sum(is.na(gene_score_df$Cluster)), " cells did not receive a cluster label.")
}

plots_dir <- file.path(getOutputDirectory(proj), "Plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

write.xlsx(
  gene_score_df,
  file = file.path(plots_dir, "PDO_gene_scores_19_genes.xlsx"),
  overwrite = TRUE
)

# =============================
# 3. K=5 dot plot
# =============================
cluster_levels <- c("K1", "K2", "K3", "K4", "K5")
cluster_cols <- c(
  "K1" = "#E41A1C",
  "K2" = "#252A6B",
  "K3" = "#1B8A3A",
  "K4" = "#8E2492",
  "K5" = "#E69F00"
)

gene_score_df <- gene_score_df %>%
  filter(Cluster %in% cluster_levels) %>%
  mutate(Cluster = factor(Cluster, levels = cluster_levels))

genes_plot <- genes_plot[genes_plot %in% colnames(gene_score_df)]

gene_long <- gene_score_df %>%
  select(Cell, Cluster, all_of(genes_plot)) %>%
  pivot_longer(
    cols = all_of(genes_plot),
    names_to = "Gene",
    values_to = "GeneScore"
  )

gene_long$Gene <- factor(gene_long$Gene, levels = rev(genes_plot))

gene_long <- gene_long %>%
  group_by(Gene) %>%
  mutate(GeneScore_Z = as.numeric(scale(GeneScore))) %>%
  ungroup()

gene_long$GeneScore_Z[!is.finite(gene_long$GeneScore_Z)] <- 0

dot_df <- gene_long %>%
  group_by(Cluster, Gene) %>%
  summarise(
    mean_score = mean(GeneScore, na.rm = TRUE),
    mean_z = mean(GeneScore_Z, na.rm = TRUE),
    pct_nonzero = mean(GeneScore > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

p_dot <- ggplot(dot_df, aes(x = Cluster, y = Gene)) +
  geom_point(aes(size = pct_nonzero, color = mean_z)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Mean\nz-score"
  ) +
  scale_size_continuous(
    range = c(1, 8),
    limits = c(0, 100),
    name = "% cells\nscore > 0"
  ) +
  theme_classic(base_size = 14) +
  labs(x = "K-means cluster", y = NULL) +
  theme(
    axis.text.x = element_text(color = cluster_cols[cluster_levels]),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(plots_dir, "PDO_19gene_score_dotplot_K5.pdf"),
  plot = p_dot,
  width = 6,
  height = 4,
  units = "in"
)

write.xlsx(
  dot_df,
  file = file.path(plots_dir, "PDO_19gene_dotplot_summary_K5.xlsx"),
  overwrite = TRUE
)

message("PDO manuscript analysis complete. Outputs: ", plots_dir)
