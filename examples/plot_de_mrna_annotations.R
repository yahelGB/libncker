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
# Paths — override via command-line args:
#   Rscript plot_de_mrna_annotations.R <annotations_dir> <output_dir>
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
annotations_dir <- if (length(args) >= 1) args[1] else "results/annotations"
output_dir      <- if (length(args) >= 2) args[2] else "results"

out_plot_categories <- file.path(output_dir, "de_mrna_functional_categories.png")
out_plot_proportion <- file.path(output_dir, "de_mrna_characterized_proportion.png")
out_plot_top        <- file.path(output_dir, "de_mrna_top_targeted.png")
out_tsv_classified  <- file.path(output_dir, "de_mrna_classified.tsv")

top_n_mrnas <- 10

# -----------------------------
# Read annotation files
# -----------------------------
ann_files <- list.files(
  annotations_dir,
  pattern = "_de_mrna_annotations\\.tsv$",
  full.names = TRUE
)

if (length(ann_files) == 0) {
  stop("No *_de_mrna_annotations.tsv files found in: ", annotations_dir)
}

ann_all <- lapply(ann_files, function(path) {
  tissue <- sub("_de_mrna_annotations\\.tsv$", "", basename(path))
  tbl <- fread(path, sep = "\t", header = TRUE, data.table = FALSE,
               na.strings = c("", "NA"))
  tbl$tissue <- tissue
  tbl
}) %>%
  bind_rows() %>%
  mutate(
    n_lncrnas            = lengths(strsplit(lncRNA_IDs, ";")),
    go_biological_process = ifelse(is.na(go_biological_process), "", go_biological_process),
    go_molecular_function = ifelse(is.na(go_molecular_function), "", go_molecular_function),
    go_cellular_component = ifelse(is.na(go_cellular_component), "", go_cellular_component)
  )

# -----------------------------
# Assign functional category
#
# Priority:
#   1. Use GO Biological Process terms (controlled vocabulary — most reliable)
#   2. Fall back to keyword matching on product name when GO BP is absent
# -----------------------------

# Map GO BP terms to broad categories (GO slim-like grouping)
go_bp_category <- function(go_bp) {
  g <- tolower(go_bp)
  case_when(
    grepl("immune|defense|innate|inflammatory|apoptot|phagocyt|autophagy", g)    ~ "Immune / defense",
    grepl("translation|ribosom", g)                                               ~ "Translation / ribosome",
    grepl("transcription|chromatin|histone|gene expression|epigenetic", g)        ~ "Transcription / chromatin",
    grepl("rna process|mrna|splicing|rna metabol|rna export", g)                  ~ "RNA processing",
    grepl("dna repair|dna replication|cell cycle|mitosis|meiosis|chromosome", g)  ~ "Cell cycle / DNA",
    grepl("proteolysis|ubiquitin|proteasom|protein degradation|protein catabol",g) ~ "Protein degradation",
    grepl("signal transduction|signaling|phosphorylation|kinase|gtpase", g)       ~ "Signaling",
    grepl("transport|transmembrane|ion channel|secretion|vesicl|trafficking", g)  ~ "Transport",
    grepl("cytoskeleton|actin|microtubule|myosin|motor", g)                       ~ "Cytoskeleton / motor",
    grepl("metabolic process|biosynthetic|catabolic|oxidation|lipid|fatty acid|carbohydrate|amino acid|glucos|glycolysis|citric acid|oxidative phosph", g) ~ "Metabolism",
    grepl("development|morphogenesis|differentiation|organogenesis", g)           ~ "Development",
    TRUE                                                                           ~ NA_character_
  )
}

# Fallback: keyword matching on product name when no GO BP is available
product_category <- function(product) {
  p <- tolower(product)
  case_when(
    grepl("uncharacterized|hypothetical", p)                                      ~ "Uncharacterized",
    grepl("actin|myosin|obscurin|tubulin|titin|filamin|tropomyosin|unc-89", p)   ~ "Cytoskeleton / motor",
    grepl("collagen|mucin|laminin|adhesin|fibronectin|integrin|cadherin", p)      ~ "Transport",
    grepl("ubiquitin|proteasome|ring finger", p)                                  ~ "Protein degradation",
    grepl("ribosom|argonaute", p)                                                 ~ "Translation / ribosome",
    grepl("rna|spliceosom|exonuclease|maelstrom", p)                             ~ "RNA processing",
    grepl("histone|methyltransferase|chromatin|bromodomain|homeobox", p)          ~ "Transcription / chromatin",
    grepl("lectin|chitinase|mannose|immune|antimicrobial|defensin", p)            ~ "Immune / defense",
    grepl("transporter|channel|abc sub|carrier|permease", p)                     ~ "Transport",
    grepl("dehydrogenase|synthase|reductase|oxidase|hydroxylase|lipase|transferase|glucuronidase|carboxylase|decarboxylase|phosphatase|glucosyl|glycosyl|kinase", p) ~ "Metabolism",
    grepl("receptor|gtpase|ras |erk|mapk", p)                                    ~ "Signaling",
    TRUE                                                                          ~ "No GO annotation"
  )
}

classify_mrna <- function(go_bp, product) {
  # Try GO BP first — take the first term that maps to a known category
  if (!is.na(go_bp) && nchar(go_bp) > 0) {
    terms <- strsplit(go_bp, ";")[[1]]
    for (term in terms) {
      cat <- go_bp_category(term)
      if (!is.na(cat)) return(cat)
    }
    return("No GO slim match")  # has GO BP but none mapped to our categories
  }
  # Fallback to product name
  product_category(product)
}

ann_classified <- ann_all %>%
  rowwise() %>%
  mutate(category = classify_mrna(go_biological_process, product)) %>%
  ungroup()

# -----------------------------
# Plot 1: functional category bar chart per tissue
# -----------------------------
cat_counts <- ann_classified %>%
  count(tissue, category) %>%
  group_by(tissue) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

category_order <- cat_counts %>%
  group_by(category) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(category)

cat_counts <- cat_counts %>%
  mutate(category = factor(category, levels = category_order))

category_colors <- c(
  "Uncharacterized"          = "#b0b0b0",
  "No GO annotation"         = "#d3d3d3",
  "No GO slim match"         = "#e8e8e8",
  "Metabolism"               = "#e6ab02",
  "Signaling"                = "#1b9e77",
  "Cytoskeleton / motor"     = "#d95f02",
  "Translation / ribosome"   = "#7570b3",
  "RNA processing"           = "#a6761d",
  "Transcription / chromatin" = "#66a61e",
  "Transport"                = "#56b4e9",
  "Protein degradation"      = "#e7298a",
  "Immune / defense"         = "#e41a1c",
  "Cell cycle / DNA"         = "#984ea3",
  "Development"              = "#ff7f00"
)

p_categories <- ggplot(cat_counts, aes(x = category, y = n, fill = category)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = n), vjust = -0.4, size = 3) +
  facet_wrap(~ tissue, scales = "free_y") +
  scale_fill_manual(values = category_colors, guide = "none") +
  labs(
    title = "Functional categories of DE mRNAs per tissue",
    subtitle = "Categories from GO Biological Process terms; product name used as fallback",
    x = NULL,
    y = "Number of DE mRNAs"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x    = element_text(angle = 40, hjust = 1, size = 9),
    strip.text     = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

# -----------------------------
# Plot 2: annotation source breakdown per tissue
# (GO-annotated vs product-name fallback vs uncharacterized)
# -----------------------------
prop_data <- ann_classified %>%
  mutate(
    source = case_when(
      category == "Uncharacterized"  ~ "Uncharacterized",
      go_biological_process != ""    ~ "GO-annotated",
      TRUE                           ~ "Named, no GO"
    )
  ) %>%
  count(tissue, source) %>%
  group_by(tissue) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(source = factor(source, levels = c("GO-annotated", "Named, no GO", "Uncharacterized")))

p_proportion <- ggplot(prop_data, aes(x = tissue, y = pct, fill = source)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(
    aes(label = paste0(round(pct, 1), "%\n(n=", n, ")")),
    position = position_stack(vjust = 0.5),
    size = 3.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "GO-annotated"   = "#1b9e77",
      "Named, no GO"   = "#d95f02",
      "Uncharacterized" = "#b0b0b0"
    )
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = "Annotation coverage of DE mRNAs per tissue",
    subtitle = "GO-annotated = has GO Biological Process terms; Named, no GO = known product but no GO",
    x = NULL,
    y = "Percentage of DE mRNAs",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.x = element_blank()
  )

# -----------------------------
# Plot 3: top N most-targeted characterized mRNAs per tissue (lollipop)
# -----------------------------
top_targeted <- ann_classified %>%
  filter(!category %in% c("Uncharacterized", "No GO annotation", "No GO slim match")) %>%
  group_by(tissue) %>%
  slice_max(order_by = n_lncrnas, n = top_n_mrnas, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    short_label = ifelse(
      nchar(product) > 45,
      paste0(substr(product, 1, 43), ".."),
      product
    )
  ) %>%
  arrange(tissue, n_lncrnas) %>%
  mutate(short_label = factor(short_label, levels = unique(short_label)))

p_top <- ggplot(top_targeted, aes(x = n_lncrnas, y = short_label, color = tissue)) +
  geom_segment(aes(x = 0, xend = n_lncrnas, yend = short_label), linewidth = 0.7) +
  geom_point(size = 3) +
  facet_wrap(~ tissue, scales = "free_y", ncol = 2) +
  scale_color_manual(
    values = c(
      "gills"     = "#1b9e77",
      "hemocytes" = "#d95f02",
      "hep"       = "#7570b3",
      "muscle"    = "#e7298a"
    ),
    guide = "none"
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  labs(
    title = paste0("Top ", top_n_mrnas, " most-targeted DE mRNAs per tissue"),
    subtitle = "Ranked by number of associated lncRNAs; uncharacterized excluded",
    x = "Number of associated lncRNAs",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text         = element_text(face = "bold"),
    axis.text.y        = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# -----------------------------
# Save outputs
# -----------------------------
fwrite(ann_classified, out_tsv_classified, sep = "\t")
ggsave(out_plot_categories, p_categories, width = 13, height = 8,  dpi = 300)
ggsave(out_plot_proportion, p_proportion, width = 8,  height = 6,  dpi = 300)
ggsave(out_plot_top,        p_top,        width = 13, height = 10, dpi = 300)

# -----------------------------
# Console output
# -----------------------------
cat("\nTissues loaded:", paste(unique(ann_all$tissue), collapse = ", "), "\n")
cat("Total unique mRNA rows:", nrow(ann_all), "\n")
cat("Rows with GO BP terms:", sum(ann_all$go_biological_process != ""), "\n")

cat("\nCategory counts (all tissues combined):\n")
print(
  ann_classified %>%
    count(category, sort = TRUE) %>%
    mutate(pct = round(n / sum(n) * 100, 1))
)

cat("\nSaved classified table to:", out_tsv_classified, "\n")
cat("Saved category bar chart to:", out_plot_categories, "\n")
cat("Saved annotation coverage plot to:", out_plot_proportion, "\n")
cat("Saved top-targeted lollipop to:", out_plot_top, "\n")
