#!/usr/bin/env Rscript

# Proteo-Seq core ArchR workflow
#
# Edit the configuration block before running.

suppressPackageStartupMessages({
  library(ArchR)
  library(openxlsx)
})

# =============================
# User configuration
# =============================
fragments_file <- "/path/to/fragments.tsv.gz"
sample_name <- "Device7"                    # e.g. Device7 or Device11
output_dir <- "/path/to/ArchR_analysis"

threads <- 8
genome <- "hg38"

# Arrow creation uses permissive thresholds first; manuscript QC is applied below.
arrow_min_tss <- 0
arrow_min_frags <- 10

# Manuscript QC thresholds.
qc_tss_min <- 1
qc_frags_min <- 1000

# Final K-means setting for the dataset being analyzed.
# Manuscript analyses used K=4 for MCF7 and K=5 for PDO.
k_value <- 4
kmeans_seed <- 123
kmeans_nstart <- 50

# LSI/UMAP parameters.
lsi_var_features <- 2000
lsi_dims <- 1:10
umap_neighbors <- 30
umap_min_dist <- 0.5
umap_metric <- "cosine"

# =============================
# Setup
# =============================
stopifnot(file.exists(fragments_file))
stopifnot(file.exists(paste0(fragments_file, ".tbi")))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

addArchRGenome(genome)
addArchRThreads(threads = threads)

# Remove only ArchR temporary folders from interrupted runs.
unlink(c("tmp", "ArchRTemp"), recursive = TRUE, force = TRUE)

# =============================
# 1. Create Arrow file
# =============================
input_files <- setNames(fragments_file, sample_name)

ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = sample_name,
  minTSS = arrow_min_tss,
  minFrags = arrow_min_frags,
  addTileMat = FALSE,
  addGeneScoreMat = FALSE
)

print(ArrowFiles)

# =============================
# 2. Create ArchR project
# =============================
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = file.path(output_dir, "ArchRProject_all_cells"),
  copyArrows = TRUE
)

message("Initial project:")
print(proj)

# Export all-cell metadata for traceability.
all_meta <- as.data.frame(getCellColData(proj))
all_meta$Cell <- rownames(all_meta)
write.xlsx(
  all_meta,
  file = file.path(output_dir, "ArchR_all_cells_metadata.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

# =============================
# 3. Apply manuscript QC
# =============================
cell_meta <- as.data.frame(getCellColData(proj))

cells_keep <- rownames(cell_meta)[
  cell_meta$TSSEnrichment > qc_tss_min &
    cell_meta$nFrags > qc_frags_min
]

message(
  "Cells retained after QC: ", length(cells_keep), " / ", nrow(cell_meta),
  " (TSSEnrichment > ", qc_tss_min,
  "; nFrags > ", qc_frags_min, ")"
)

filtered_dir <- file.path(output_dir, "ArchRProject_QCfiltered")
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells_keep,
  outputDirectory = filtered_dir,
  force = TRUE
)

# QC plots.
p_frag <- plotFragmentSizes(proj)
p_tss <- plotTSSEnrichment(proj)
plotPDF(
  p_frag,
  p_tss,
  name = "QC_fragment_size_and_TSS.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)

filtered_meta <- as.data.frame(getCellColData(proj))
filtered_meta$Cell <- rownames(filtered_meta)
write.xlsx(
  filtered_meta,
  file = file.path(filtered_dir, "ArchR_QCfiltered_metadata.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

# =============================
# 4. TileMatrix + Iterative LSI
# =============================
proj <- addTileMatrix(proj)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 1,
  varFeatures = lsi_var_features,
  dimsToUse = lsi_dims,
  clusterParams = list(
    resolution = NULL,
    sampleCells = NULL,
    n.start = 10
  )
)

# =============================
# 5. K-means clustering
# =============================
lsi_mat <- getReducedDims(proj, "IterativeLSI")

set.seed(kmeans_seed)
km_res <- kmeans(lsi_mat, centers = k_value, nstart = kmeans_nstart)
clust_labels <- paste0("K", km_res$cluster)
cluster_col <- paste0("ClustersKmeans_K", k_value)

proj <- addCellColData(
  ArchRProj = proj,
  data = clust_labels,
  name = cluster_col,
  cells = rownames(lsi_mat),
  force = TRUE
)

message("Cluster sizes for K=", k_value)
print(table(getCellColData(proj, cluster_col)[, 1]))

# =============================
# 6. UMAP
# =============================
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = umap_neighbors,
  minDist = umap_min_dist,
  metric = umap_metric,
  force = TRUE
)

p_umap <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = cluster_col,
  embedding = "UMAP"
)

plotPDF(
  p_umap,
  name = paste0("UMAP_", cluster_col, ".pdf"),
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)

# =============================
# 7. GeneScoreMatrix
# =============================
proj <- addGeneScoreMatrix(proj)
message("Available matrices: ", paste(getAvailableMatrices(proj), collapse = ", "))

# =============================
# 8. Save final project
# =============================
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = getOutputDirectory(proj),
  load = FALSE
)

message("Saved filtered ArchR project: ", getOutputDirectory(proj))
