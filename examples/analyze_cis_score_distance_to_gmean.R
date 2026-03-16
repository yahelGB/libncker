Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# Paths and parameters
# -----------------------------
gff_path <- "/Users/yahelgb/Desktop/projects/libncker/examples/data_test/GCF_003789085.1_ASM378908v1_genomic.gff"
output_dir <- "/Users/yahelgb/Desktop/projects/libncker/examples/output_test"

cis_score_target <- 3
selected_tolerance <- 0.10
report_tolerances <- c(0.05, 0.10, 0.25, 0.50)

out_path_full_log <- file.path(output_dir, "intergenic_distance_full_log10_with_cis_score3.png")
out_path_zoom <- file.path(output_dir, "intergenic_distance_zoom_99_with_cis_score3.png")
out_path_cis_hist <- file.path(output_dir, "cis_score3_distance_histogram.png")
out_path_rows <- file.path(output_dir, "cis_score3_distance_to_gmean_rows.tsv")
out_path_lnc <- file.path(output_dir, "cis_score3_distance_to_gmean_lnc.tsv")
out_path_summary <- file.path(output_dir, "cis_score3_distance_to_gmean_summary.tsv")

# -----------------------------
# Read GFF
# -----------------------------
lines <- readLines(gff_path)
lines <- lines[!grepl("^#", lines)]

gff <- fread(
  text = lines,
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  data.table = FALSE,
  col.names = c(
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes"
  )
)

# -----------------------------
# Extract genes
# -----------------------------
genes <- gff %>%
  filter(type == "gene") %>%
  mutate(
    ID = ifelse(
      grepl("ID=", attributes),
      sub(".*ID=([^;]+).*", "\\1", attributes),
      NA_character_
    )
  ) %>%
  filter(!is.na(ID)) %>%
  select(seqid, start, end, strand, ID)

# -----------------------------
# Compute intergenic distances
# -----------------------------
out <- list()

for (chr in unique(genes$seqid)) {
  chr_genes <- genes %>%
    filter(seqid == chr) %>%
    arrange(start)

  if (nrow(chr_genes) < 2) next

  for (i in seq_len(nrow(chr_genes) - 1)) {
    a <- chr_genes[i, ]
    b <- chr_genes[i + 1, ]

    dist <- b$start - a$end - 1

    out[[length(out) + 1]] <- data.frame(
      seqid = chr,
      gene_A = a$ID,
      gene_B = b$ID,
      distance = dist
    )
  }
}

intergenic <- bind_rows(out) %>%
  filter(distance > 0)

# -----------------------------
# Summary statistics
# -----------------------------
gmean_bp <- exp(mean(log(intergenic$distance)))
gmean_kb <- gmean_bp / 1000

summary_stats <- quantile(
  intergenic$distance,
  probs = c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
)

xmax_99 <- unname(summary_stats["99%"])

scaffold_stats <- intergenic %>%
  group_by(seqid) %>%
  summarise(
    n_intervals = n(),
    median_distance = median(distance),
    max_distance = max(distance),
    geometric_mean_distance = exp(mean(log(distance))),
    .groups = "drop"
  ) %>%
  arrange(desc(max_distance))

# -----------------------------
# Read all cis tables
# -----------------------------
cis_files <- list.files(
  output_dir,
  pattern = "_cis_regulation_module_output\\.txt$",
  full.names = TRUE
)

if (length(cis_files) == 0) {
  stop("No *cis_regulation_module_output.txt files were found in output_dir.")
}

cis_all <- lapply(cis_files, function(path) {
  tissue <- sub("_cis_regulation_module_output\\.txt$", "", basename(path))
  tbl <- fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  tbl$tissue <- tissue
  tbl
}) %>%
  bind_rows() %>%
  mutate(
    cis_score = as.numeric(cis_score),
    distance_Kb = as.numeric(distance_Kb),
    distance_bp = distance_Kb * 1000
  )

cis_score_hits <- cis_all %>%
  filter(cis_score == cis_score_target) %>%
  mutate(
    delta_bp = abs(distance_bp - gmean_bp),
    delta_Kb = abs(distance_Kb - gmean_kb),
    rel_delta = delta_bp / gmean_bp
  )

if (nrow(cis_score_hits) == 0) {
  stop(sprintf("No rows were found with cis_score = %s.", cis_score_target))
}

for (tol in report_tolerances) {
  flag_name <- paste0("within_", gsub("\\.", "_", format(tol, nsmall = 2)))
  cis_score_hits[[flag_name]] <- cis_score_hits$rel_delta <= tol
}

row_summary <- data.frame(
  tolerance = report_tolerances,
  lower_kb = gmean_kb * (1 - report_tolerances),
  upper_kb = gmean_kb * (1 + report_tolerances),
  row_count = sapply(report_tolerances, function(tol) sum(cis_score_hits$rel_delta <= tol)),
  unique_lnc_count = sapply(report_tolerances, function(tol) {
    cis_score_hits %>%
      filter(rel_delta <= tol) %>%
      summarise(n = n_distinct(lncRNA_ID)) %>%
      pull(n)
  })
)

lnc_summary <- cis_score_hits %>%
  group_by(tissue, lncRNA_ID, Contig, lncRNA_strand) %>%
  summarise(
    n_score3_matches = n(),
    min_distance_Kb = min(distance_Kb),
    max_distance_Kb = max(distance_Kb),
    mean_distance_Kb = mean(distance_Kb),
    closest_distance_Kb = distance_Kb[which.min(delta_Kb)][1],
    closest_mRNA_ID = mRNA_ID[which.min(delta_Kb)][1],
    closest_product = product[which.min(delta_Kb)][1],
    min_delta_Kb = min(delta_Kb),
    min_rel_delta = min(rel_delta),
    within_selected_tolerance = any(rel_delta <= selected_tolerance),
    .groups = "drop"
  ) %>%
  arrange(min_rel_delta, tissue, lncRNA_ID)

selected_band_rows <- cis_score_hits %>%
  filter(rel_delta <= selected_tolerance)

selected_band_lnc <- lnc_summary %>%
  filter(within_selected_tolerance)

tissue_summary <- cis_score_hits %>%
  group_by(tissue) %>%
  summarise(
    row_count = n(),
    unique_lnc_count = n_distinct(lncRNA_ID),
    rows_within_selected_tolerance = sum(rel_delta <= selected_tolerance),
    lncs_within_selected_tolerance = n_distinct(lncRNA_ID[rel_delta <= selected_tolerance]),
    .groups = "drop"
  )

summary_export <- bind_rows(
  data.frame(
    section = "global",
    metric = c(
      "geometric_mean_bp",
      "geometric_mean_kb",
      "cis_score_target",
      "selected_tolerance",
      "score3_row_count",
      "score3_unique_lnc_count",
      "rows_within_selected_tolerance",
      "lncs_within_selected_tolerance"
    ),
    value = c(
      gmean_bp,
      gmean_kb,
      cis_score_target,
      selected_tolerance,
      nrow(cis_score_hits),
      n_distinct(cis_score_hits$lncRNA_ID),
      nrow(selected_band_rows),
      nrow(selected_band_lnc)
    ),
    stringsAsFactors = FALSE
  ) %>%
    mutate(value = as.character(value)),
  row_summary %>%
    mutate(section = "tolerance_band", metric = paste0("within_", tolerance)) %>%
    transmute(
      section = section,
      metric = metric,
      value = paste0(
        "rows=", row_count,
        ";unique_lnc=", unique_lnc_count,
        ";lower_kb=", round(lower_kb, 4),
        ";upper_kb=", round(upper_kb, 4)
      )
    )
)

# -----------------------------
# Plots
# -----------------------------
p_full_log <- ggplot(intergenic, aes(x = distance)) +
  geom_histogram(
    bins = 100,
    fill = "gray",
    color = "black",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = gmean_bp,
    color = "black",
    linetype = "dotted",
    linewidth = 1.2
  ) +
  geom_rug(
    data = cis_score_hits,
    aes(x = distance_bp, color = tissue),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.75
  ) +
  scale_x_log10(limits = range(intergenic$distance, finite = TRUE)) +
  labs(
    title = "Intergenic distance distribution with cis_score = 3 hits",
    subtitle = paste0(
      "Geometric mean = ",
      round(gmean_bp, 2), " bp (",
      round(gmean_kb, 3), " Kb)"
    ),
    x = "Distance (bp, log10 scale)",
    y = "Count",
    color = "Tissue"
  ) +
  theme_minimal()

p_zoom_99 <- ggplot(intergenic, aes(x = distance)) +
  geom_histogram(
    binwidth = 100,
    fill = "#4C72B0",
    color = "black",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = gmean_bp,
    color = "darkgreen",
    linetype = "dotted",
    linewidth = 1.1
  ) +
  annotate(
    "rect",
    xmin = gmean_bp * (1 - selected_tolerance),
    xmax = gmean_bp * (1 + selected_tolerance),
    ymin = -Inf,
    ymax = Inf,
    fill = "gold",
    alpha = 0.12
  ) +
  geom_rug(
    data = cis_score_hits,
    aes(x = distance_bp, color = tissue),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.75
  ) +
  coord_cartesian(xlim = c(0, xmax_99)) +
  labs(
    title = "Intergenic distances zoomed to the 99th percentile",
    subtitle = paste0(
      "Gold band = +/- ",
      round(selected_tolerance * 100), "% around the geometric mean"
    ),
    x = "Distance (bp)",
    y = "Count",
    color = "Tissue"
  ) +
  theme_minimal()

p_cis_hist <- ggplot(cis_score_hits, aes(x = distance_bp, fill = tissue)) +
  geom_histogram(
    bins = 40,
    color = "black",
    alpha = 0.7,
    position = "identity"
  ) +
  geom_vline(
    xintercept = gmean_bp,
    color = "black",
    linetype = "dotted",
    linewidth = 1.2
  ) +
  annotate(
    "rect",
    xmin = gmean_bp * (1 - selected_tolerance),
    xmax = gmean_bp * (1 + selected_tolerance),
    ymin = -Inf,
    ymax = Inf,
    fill = "gold",
    alpha = 0.12
  ) +
  labs(
    title = paste0("cis_score = ", cis_score_target, " distances"),
    subtitle = paste0(
      "Rows within +/- ",
      round(selected_tolerance * 100),
      "% of geometric mean: ",
      nrow(selected_band_rows),
      " | Unique lncRNAs: ",
      nrow(selected_band_lnc)
    ),
    x = "Distance (bp)",
    y = "Count",
    fill = "Tissue"
  ) +
  theme_minimal()

# -----------------------------
# Save outputs
# -----------------------------
fwrite(cis_score_hits, out_path_rows, sep = "\t")
fwrite(lnc_summary, out_path_lnc, sep = "\t")
fwrite(summary_export, out_path_summary, sep = "\t")

ggsave(out_path_full_log, p_full_log, width = 10, height = 6, dpi = 300)
ggsave(out_path_zoom, p_zoom_99, width = 10, height = 6, dpi = 300)
ggsave(out_path_cis_hist, p_cis_hist, width = 10, height = 6, dpi = 300)

# -----------------------------
# Console output
# -----------------------------
cat("\nSummary of intergenic distances:\n")
print(summary_stats)

cat("\nGeometric mean (bp):", round(gmean_bp, 2), "\n")
cat("Geometric mean (Kb):", round(gmean_kb, 3), "\n")

cat("\nTissue summary for cis_score =", cis_score_target, ":\n")
print(tissue_summary)

cat("\nTolerance summary:\n")
print(row_summary)

cat(
  "\nSelected tolerance (+/-", round(selected_tolerance * 100),
  "%) rows:", nrow(selected_band_rows),
  "| unique lncRNAs:", nrow(selected_band_lnc), "\n"
)

cat("\nTop lncRNAs closest to the geometric mean:\n")
print(head(lnc_summary, 10))

cat("\nSaved row-level table to:", out_path_rows, "\n")
cat("Saved lnc-level table to:", out_path_lnc, "\n")
cat("Saved summary table to:", out_path_summary, "\n")
cat("Saved full-range log10 plot to:", out_path_full_log, "\n")
cat("Saved zoomed plot to:", out_path_zoom, "\n")
cat("Saved cis_score histogram to:", out_path_cis_hist, "\n")
