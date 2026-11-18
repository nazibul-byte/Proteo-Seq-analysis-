#!/usr/bin/env Rscript

# ESR1 chromatin-accessibility coverage analysis for the MCF7 manuscript dataset.

suppressPackageStartupMessages({
  library(ArchR)
})

# =============================
# User configuration
# =============================
archr_project_dir <- "/path/to/MCF7/ArchRProject_QCfiltered"

proj <- loadArchRProject(path = archr_project_dir)
n_total <- length(getCellNames(proj))
proj$AllCells <- rep(paste0("All_", n_total, "_cells"), n_total)

# Aggregated ESR1 accessibility track with 20-kb flanking sequence and gene annotation.
p_ESR1 <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "AllCells",
  geneSymbol = "ESR1",
  upstream = 20000,
  downstream = 20000,
  plotSummary = c("bulkTrack", "geneTrack"),
  tileSize = 100,
  maxCells = 500,
  normMethod = "ReadsInTSS",
  title = "ESR1 locus"
)

plotPDF(
  p_ESR1,
  name = "ESR1_AllCells_full_locus_coverage.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 9,
  height = 5
)

message("Saved ESR1 coverage plot to: ", file.path(getOutputDirectory(proj), "Plots"))
