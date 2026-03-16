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
# Paths
# -----------------------------
gff_path <- "/Users/yahelgb/Desktop/projects/libncker/examples/data_test/GCF_003789085.1_ASM378908v1_genomic.gff"
output_dir <- "/Users/yahelgb/Desktop/projects/libncker/examples/output_test"

out_plot_lnc <- file.path(output_dir, "cis_score3_first_neighbor_lncRNA_distances_vs_gmean.png")
out_plot_rows <- file.path(output_dir, "cis_score3_first_neighbor_all_rows_vs_gmean.png")
out_tsv_lnc <- file.path(output_dir, "cis_score3_first_neighbor_lncRNA_distances_vs_gmean.tsv")
out_tsv_rows <- file.path(output_dir, "cis_score3_first_neighbor_all_rows_vs_gmean.tsv")

cis_score_target <- 3
first_neighbor_positions <- c("1_U", "1_D")

# -----------------------------
# Read GFF and compute intergenic geometric mean
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

gmean_bp <- exp(mean(log(intergenic$distance)))
gmean_kb <- gmean_bp / 1000

# -----------------------------
# Read cis output files
# -----------------------------
cis_files <- list.files(
  output_dir,
  pattern = "_cis_regulation_module_output\\.txt$",
  full.names = TRUE
)

if (length(cis_files) == 0) {
  stop("No *cis_regulation_module_output.txt files were found.")
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
    distance_Kb = as.numeric(distance_Kb)
  )

first_neighbor_rows <- cis_all %>%
  filter(
    cis_score == cis_score_target,
    mRNA_position %in% first_neighbor_positions
  ) %>%
  mutate(
    distance_class = case_when(
      distance_Kb <= 0 ~ "Overlap (<= 0 Kb)",
      distance_Kb < gmean_kb ~ "Minor (< geometric mean)",
      TRUE ~ "Not minor (>= geometric mean)"
    ),
    abs_distance_kb = abs(distance_Kb)
  ) %>%
  arrange(distance_Kb, tissue, lncRNA_ID)

if (nrow(first_neighbor_rows) == 0) {
  stop("No rows matched cis_score = 3 and positions 1_U / 1_D.")
}

# One representative point per lncRNA: the closest first neighbor to zero.
first_neighbor_lnc <- first_neighbor_rows %>%
  group_by(tissue, lncRNA_ID) %>%
  arrange(abs_distance_kb, distance_Kb, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    minor_vs_gmean = ifelse(distance_Kb < gmean_kb, "Below geometric mean", "At/above geometric mean"),
    lnc_label = paste(tissue, lncRNA_ID, sep = " | ")
  ) %>%
  arrange(distance_Kb) %>%
  mutate(lnc_label = factor(lnc_label, levels = lnc_label))

summary_table <- first_neighbor_lnc %>%
  count(minor_vs_gmean, name = "unique_lnc_count")

row_summary <- first_neighbor_rows %>%
  count(distance_class, name = "row_count")

# -----------------------------
# Plot 1: one point per lncRNA
# -----------------------------
p_lnc <- ggplot(first_neighbor_lnc, aes(x = distance_Kb, y = lnc_label, color = minor_vs_gmean)) +
  geom_vline(
    xintercept = gmean_kb,
    color = "black",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_point(size = 2.6, alpha = 0.9) +
  facet_wrap(~ tissue, scales = "free_y", ncol = 1) +
  scale_color_manual(
    values = c(
      "Below geometric mean" = "#1b9e77",
      "At/above geometric mean" = "#d95f02"
    )
  ) +
  labs(
    title = "Closest first-neighbor cis distances per lncRNA",
    subtitle = paste0(
      "Filtered to cis_score = ", cis_score_target,
      " and positions ", paste(first_neighbor_positions, collapse = ", "),
      " | Geometric mean = ", round(gmean_kb, 3), " Kb"
    ),
    x = "Distance to DE mRNA (Kb)",
    y = "lncRNA",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 7)
  )

# -----------------------------
# Plot 2: all first-neighbor rows
# -----------------------------
p_rows <- ggplot(first_neighbor_rows, aes(x = distance_Kb, fill = distance_class)) +
  geom_histogram(binwidth = 2, color = "black", alpha = 0.8) +
  geom_vline(
    xintercept = gmean_kb,
    color = "black",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "Overlap (<= 0 Kb)" = "#7570b3",
      "Minor (< geometric mean)" = "#1b9e77",
      "Not minor (>= geometric mean)" = "#d95f02"
    )
  ) +
  labs(
    title = "All first-neighbor cis distances",
    subtitle = paste0(
      "Rows below the geometric mean (including overlaps) visually separated at ",
      round(gmean_kb, 3), " Kb"
    ),
    x = "Distance to DE mRNA (Kb)",
    y = "Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 11)

# -----------------------------
# Save outputs
# -----------------------------
fwrite(first_neighbor_rows, out_tsv_rows, sep = "\t")
fwrite(first_neighbor_lnc, out_tsv_lnc, sep = "\t")
ggsave(out_plot_lnc, p_lnc, width = 11, height = 12, dpi = 300)
ggsave(out_plot_rows, p_rows, width = 10, height = 6, dpi = 300)

# -----------------------------
# Console output
# -----------------------------
cat("\nGeometric mean (Kb):", round(gmean_kb, 6), "\n")
cat("Unique lncRNAs with cis_score = 3 and 1_U/1_D:", nrow(first_neighbor_lnc), "\n")
cat("Unique lncRNAs below geometric mean:", sum(first_neighbor_lnc$distance_Kb < gmean_kb), "\n")
cat("Unique lncRNAs at or above geometric mean:", sum(first_neighbor_lnc$distance_Kb >= gmean_kb), "\n")
cat("\nUnique-lnc summary:\n")
print(summary_table)
cat("\nRow summary:\n")
print(row_summary)
cat("\nSaved lncRNA plot to:", out_plot_lnc, "\n")
cat("Saved row plot to:", out_plot_rows, "\n")
cat("Saved lncRNA table to:", out_tsv_lnc, "\n")
cat("Saved row table to:", out_tsv_rows, "\n")
