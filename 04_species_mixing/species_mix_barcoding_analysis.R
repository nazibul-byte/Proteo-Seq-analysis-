#!/usr/bin/env Rscript

# Proteo-Seq species-mixing / barcoding validation.
#
# Inputs are per-well mapped-read tables with columns Sample and Mapped.
# Sample names encode the designed species and well, for example:
#   d4_human_r11_c10_REP1
#   d4_mouse_r8_c7_REP1

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
})

# Locate this script so the bundled input files work regardless of the
# directory from which Rscript is called.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  script_dir <- getwd()
}

# =============================
# User configuration
# =============================
human_map_xlsx <- file.path(script_dir, "input_data", "device4_human_map.xlsx")
mouse_map_xlsx <- file.path(script_dir, "input_data", "device4_mouse_map.xlsx")
output_dir <- file.path(script_dir, "output")

expected_rows <- 16
expected_cols <- 32
min_reads <- 10000
ratio_threshold <- 5

# =============================
# 1. Read and validate mapping-count tables
# =============================
read_mapping_table <- function(path, mapped_name) {
  if (!file.exists(path)) stop("Input file not found: ", path)

  df <- openxlsx::read.xlsx(path)
  required <- c("Sample", "Mapped")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s) in ", path, ": ", paste(missing_cols, collapse = ", "))
  }

  out <- data.frame(
    Sample = as.character(df$Sample),
    value = as.numeric(df$Mapped),
    stringsAsFactors = FALSE
  )
  names(out)[2] <- mapped_name
  out
}

human_df <- read_mapping_table(human_map_xlsx, "human_map")
mouse_df <- read_mapping_table(mouse_map_xlsx, "mouse_map")

if (anyDuplicated(human_df$Sample)) stop("Duplicate Sample values in human mapping table.")
if (anyDuplicated(mouse_df$Sample)) stop("Duplicate Sample values in mouse mapping table.")
if (!setequal(human_df$Sample, mouse_df$Sample)) {
  stop("Human- and mouse-reference mapping tables do not contain the same Sample IDs.")
}

analysis_df <- full_join(human_df, mouse_df, by = "Sample")

# =============================
# 2. Parse experimental identity and well coordinates
# =============================
analysis_df$sample_name <- sub("_REP[0-9]+$", "", analysis_df$Sample)
valid_name <- grepl(
  "^d[0-9]+_(human|mouse)_r[0-9]+_c[0-9]+$",
  analysis_df$sample_name
)
if (!all(valid_name)) {
  stop("One or more Sample IDs do not match the expected d#_{human|mouse}_r#_c#_REP# naming pattern.")
}

analysis_df$expected_species <- sub(
  "^d[0-9]+_(human|mouse)_.*$", "\\1", analysis_df$sample_name
)
analysis_df$row <- as.integer(sub(
  "^.*_r([0-9]+)_c[0-9]+$", "\\1", analysis_df$sample_name
))
analysis_df$column <- as.integer(sub(
  "^.*_c([0-9]+)$", "\\1", analysis_df$sample_name
))

expected_n <- expected_rows * expected_cols
if (nrow(analysis_df) != expected_n) {
  stop("Expected ", expected_n, " wells, found ", nrow(analysis_df), ".")
}
if (any(analysis_df$row < 1 | analysis_df$row > expected_rows) ||
    any(analysis_df$column < 1 | analysis_df$column > expected_cols)) {
  stop("Parsed row/column values fall outside the expected device dimensions.")
}
if (anyDuplicated(analysis_df[, c("row", "column")])) {
  stop("Duplicate physical well coordinates detected.")
}

expected_from_row <- ifelse(analysis_df$row %% 2 == 1, "human", "mouse")
if (!all(analysis_df$expected_species == expected_from_row)) {
  stop("Sample-name species labels are inconsistent with the alternating-row experimental design.")
}

# =============================
# 3. Apply low-signal filter and calculate species specificity
# =============================
analysis_df <- analysis_df %>%
  mutate(
    both_low = !is.na(human_map) & !is.na(mouse_map) &
      human_map < min_reads & mouse_map < min_reads,
    human_map_filtered = ifelse(both_low, NA_real_, human_map),
    mouse_map_filtered = ifelse(both_low, NA_real_, mouse_map),
    pass_10k_filter = !both_low,
    human_reads_per_1e5 = human_map_filtered / 1e5,
    mouse_reads_per_1e5 = mouse_map_filtered / 1e5,
    expected_to_other_ratio = case_when(
      expected_species == "human" ~ human_map_filtered / mouse_map_filtered,
      expected_species == "mouse" ~ mouse_map_filtered / human_map_filtered,
      TRUE ~ NA_real_
    ),
    ratio_gt_threshold = expected_to_other_ratio > ratio_threshold
  ) %>%
  select(
    sample_name, row, column, expected_species,
    human_map, mouse_map,
    pass_10k_filter,
    human_reads_per_1e5, mouse_reads_per_1e5,
    expected_to_other_ratio, ratio_gt_threshold
  ) %>%
  arrange(row, column)

# When the off-species count is zero, R returns Inf; this is retained.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.xlsx(
  analysis_df,
  file = file.path(output_dir, "SpeciesMix_device4_mapping_analysis.xlsx"),
  overwrite = TRUE
)
write.table(
  analysis_df,
  file = file.path(output_dir, "SpeciesMix_device4_mapping_analysis.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =============================
# 4. Human-versus-mouse mapped-read scatter plot
# =============================
p_scatter <- ggplot(
  analysis_df,
  aes(
    x = human_reads_per_1e5,
    y = mouse_reads_per_1e5,
    colour = expected_species
  )
) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(
    values = c(human = "thistle", mouse = "lightpink1"),
    breaks = c("human", "mouse"),
    labels = c("Human", "Mouse")
  ) +
  labs(
    x = "Human mapped reads (×10^5)",
    y = "Mouse mapped reads (×10^5)",
    colour = "Expected species"
  ) +
  theme_minimal()

ggsave(
  file.path(output_dir, "SpeciesMix_device4_human_vs_mouse_mapped_reads.png"),
  plot = p_scatter,
  width = 5,
  height = 4,
  units = "in"
)

# =============================
# 5. 16 x 32 expected/off-species ratio heatmap
# =============================
color_mapping <- function(ratio) {
  if (is.na(ratio)) {
    "white"
  } else if (ratio <= 2) {
    "#FFB3B3"
  } else if (ratio <= 3) {
    "#FFCC99"
  } else if (ratio <= 4) {
    "#FFFF99"
  } else if (ratio <= ratio_threshold) {
    "#B3E5B3"
  } else {
    "#66FF66"
  }
}

heatmap_data <- analysis_df %>%
  mutate(fill_color = vapply(expected_to_other_ratio, color_mapping, character(1)))

valid_n <- sum(!is.na(heatmap_data$expected_to_other_ratio))
pass_n <- sum(heatmap_data$expected_to_other_ratio > ratio_threshold, na.rm = TRUE)
pass_pct <- if (valid_n == 0) 0 else 100 * pass_n / valid_n
pct_label <- sprintf("Ratio > %g: %.1f%%", ratio_threshold, pass_pct)

p_heatmap <- ggplot(
  heatmap_data,
  aes(x = column, y = row, fill = fill_color)
) +
  geom_tile(color = "black") +
  geom_text(
    aes(label = ifelse(is.na(expected_to_other_ratio), "NA", sprintf("%.2f", expected_to_other_ratio))),
    size = 1.5,
    colour = "black"
  ) +
  scale_fill_identity() +
  scale_y_reverse() +
  annotate(
    "text",
    x = expected_cols + 0.4,
    y = expected_rows + 0.4,
    label = pct_label,
    hjust = 1,
    vjust = 1,
    size = 3.5
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_dir, "SpeciesMix_device4_ratio_heatmap.png"),
  plot = p_heatmap,
  width = 6,
  height = 4,
  units = "in"
)

message(
  "Species-mixing analysis complete: ",
  pass_n, "/", valid_n,
  " valid wells have expected/off-species ratio > ",
  ratio_threshold,
  " (", sprintf("%.1f", pass_pct), "%)."
)
