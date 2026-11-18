#!/usr/bin/env Rscript

# MCF7 manuscript-specific GeneScore analysis

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

# =============================
# User configuration
# =============================
archr_project_dir <- "/path/to/MCF7/ArchRProject_QCfiltered"
cluster_col <- "ClustersKmeans_K4"

# Selected genes used for MCF7 GeneScore export and visualization.
selected_genes <- c(
  "ESR1", "PGR", "ERBB2", "GATA3", "FOXA1",
  "EPCAM", "KRT8", "KRT18", "HSP90AA1"
)

dotplot_genes <- c(
  "ESR1", "PGR", "GATA3", "FOXA1",
  "EPCAM", "KRT8", "HSP90AA1"
)

er_genes_requested <- c("ESR1", "PGR", "GATA3", "FOXA1")

cluster_levels <- c("K1", "K2", "K3", "K4")
cluster_cols <- c(
  "K1" = "#E41A1C",
  "K2" = "#252A6B",
  "K3" = "#1B8A3A",
  "K4" = "#8E2492"
)

# =============================
# Load project and GeneScoreMatrix
# =============================
proj <- loadArchRProject(path = archr_project_dir)

if (!cluster_col %in% colnames(getCellColData(proj))) {
  stop("Cluster column not found: ", cluster_col)
}

gse <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gene_names <- as.character(rowData(gse)$name)
gene_names_upper <- toupper(gene_names)

selected_genes <- unique(toupper(selected_genes))
idx <- match(selected_genes, gene_names_upper)
found <- !is.na(idx)

if (!any(found)) {
  stop("None of the selected MCF7 genes were found in GeneScoreMatrix.")
}

if (any(!found)) {
  message("Missing selected genes: ", paste(selected_genes[!found], collapse = ", "))
}

selected_genes <- selected_genes[found]
idx <- idx[found]

mat <- assay(gse, "GeneScoreMatrix")[idx, , drop = FALSE]
rownames(mat) <- selected_genes

cells <- getCellNames(proj)
mat <- mat[, cells, drop = FALSE]

cl_df <- getCellColData(proj, select = cluster_col)
cl <- as.character(cl_df[, 1])
names(cl) <- rownames(cl_df)
cl <- cl[colnames(mat)]

# Restrict to the K1-K4 manuscript clusters.
keep <- !is.na(cl) & cl %in% cluster_levels
mat <- mat[, keep, drop = FALSE]
cl <- factor(cl[keep], levels = cluster_levels)

# =============================
# Export per-cell GeneScores
# =============================
gene_score_df <- data.frame(
  Cell = colnames(mat),
  Cluster = as.character(cl),
  t(as.matrix(mat)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_score_xlsx <- file.path(
  getOutputDirectory(proj),
  "Plots",
  "MCF7_gene_scores_selected_genes_all_cells.xlsx"
)
dir.create(dirname(gene_score_xlsx), recursive = TRUE, showWarnings = FALSE)

write.xlsx(
  gene_score_df,
  file = gene_score_xlsx,
  rowNames = FALSE,
  overwrite = TRUE
)

# =============================
# Cluster-level dot plot
# =============================
dotplot_genes <- toupper(dotplot_genes)
dotplot_genes <- dotplot_genes[dotplot_genes %in% colnames(gene_score_df)]

if (length(dotplot_genes) == 0) {
  stop("No requested dot-plot genes are present in the exported GeneScore table.")
}

gene_long <- gene_score_df %>%
  select(Cell, Cluster, all_of(dotplot_genes)) %>%
  pivot_longer(
    cols = all_of(dotplot_genes),
    names_to = "Gene",
    values_to = "GeneScore"
  ) %>%
  mutate(
    Cluster = factor(Cluster, levels = cluster_levels),
    Gene = factor(Gene, levels = rev(dotplot_genes))
  ) %>%
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

plots_dir <- file.path(getOutputDirectory(proj), "Plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(plots_dir, "MCF7_gene_score_dotplot_by_cluster.pdf"),
  plot = p_dot,
  width = 5.5,
  height = 4,
  units = "in"
)

write.xlsx(
  dot_df,
  file = file.path(plots_dir, "MCF7_gene_score_dotplot_source_values.xlsx"),
  overwrite = TRUE
)

# =============================
# Composite ER/luminal score
# =============================
er_genes <- er_genes_requested[er_genes_requested %in% colnames(gene_score_df)]
if (length(er_genes) == 0) {
  stop("None of the requested ER/luminal genes were found in gene_score_df.")
}

er_z_df <- gene_score_df %>%
  select(Cell, Cluster, all_of(er_genes)) %>%
  mutate(across(all_of(er_genes), as.numeric))

er_z_mat <- scale(as.matrix(er_z_df[, er_genes, drop = FALSE]))
rownames(er_z_mat) <- er_z_df$Cell
er_z_mat[!is.finite(er_z_mat)] <- 0

er_score_df <- er_z_df %>%
  select(Cell, Cluster) %>%
  mutate(ER_luminal_score = rowMeans(er_z_mat))

er_gene_z_df <- as.data.frame(er_z_mat, check.names = FALSE)
colnames(er_gene_z_df) <- paste0(colnames(er_gene_z_df), "_Z")

er_per_cell_export <- bind_cols(
  er_z_df %>% select(Cell, Cluster),
  er_gene_z_df
) %>%
  mutate(ER_luminal_score = er_score_df$ER_luminal_score) %>%
  arrange(Cluster, Cell)

er_cluster_summary <- er_score_df %>%
  group_by(Cluster) %>%
  summarise(
    n_cells = sum(!is.na(ER_luminal_score)),
    mean = mean(ER_luminal_score, na.rm = TRUE),
    SD = sd(ER_luminal_score, na.rm = TRUE),
    SEM = SD / sqrt(n_cells),
    median = median(ER_luminal_score, na.rm = TRUE),
    Q1 = unname(quantile(ER_luminal_score, 0.25, na.rm = TRUE)),
    Q3 = unname(quantile(ER_luminal_score, 0.75, na.rm = TRUE)),
    IQR = IQR(ER_luminal_score, na.rm = TRUE),
    minimum = min(ER_luminal_score, na.rm = TRUE),
    maximum = max(ER_luminal_score, na.rm = TRUE),
    .groups = "drop"
  )

score_definition <- data.frame(
  Parameter = c(
    "Genes requested",
    "Genes used",
    "Standardization",
    "Composite score",
    "Included clusters"
  ),
  Value = c(
    paste(er_genes_requested, collapse = ", "),
    paste(er_genes, collapse = ", "),
    "Each gene z-scored across all included cells",
    "Row mean of the available gene-level z-scores",
    paste(cluster_levels, collapse = ", ")
  )
)

write.xlsx(
  list(
    Per_cell_scores = er_per_cell_export,
    Cluster_summary = er_cluster_summary,
    Raw_gene_scores = er_z_df %>% arrange(Cluster, Cell),
    Score_definition = score_definition
  ),
  file = file.path(plots_dir, "MCF7_ER_luminal_score_quantitative_values.xlsx"),
  overwrite = TRUE
)

p_er <- ggplot(
  er_score_df,
  aes(
    x = factor(Cluster, levels = cluster_levels),
    y = ER_luminal_score,
    fill = Cluster
  )
) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = cluster_cols) +
  theme_classic(base_size = 14) +
  labs(
    x = "K-means cluster",
    y = "ER/luminal accessibility score",
    title = "ER-associated chromatin accessibility by cluster"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(plots_dir, "MCF7_ER_luminal_score_by_cluster.pdf"),
  plot = p_er,
  width = 4.5,
  height = 4,
  units = "in"
)

message("MCF7 GeneScore analysis complete. Outputs: ", plots_dir)
