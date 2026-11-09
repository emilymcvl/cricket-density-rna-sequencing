# summarys-for-your-convenience.R

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

cat("║   GO SYNTHESIS: density NO-SHRINK RESULTS (high vs low)  ║\n")

# ─── Directories ──────────────────────────────────────────────────────────────
sig_dir         <- "deseq2_results_density_ns/annotated/significant"
full_dir        <- "deseq2_results_density_ns/annotated/full"

go_up_dir       <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_UP_annotated"
go_down_dir     <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_DOWN_annotated"
go_up_red_dir   <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_UP_reduced"
go_down_red_dir <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_DOWN_reduced"

output_dir <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns"

dir.create(output_dir,                                       showWarnings = FALSE)
dir.create(paste0(output_dir, "/01_individual_comparisons"), showWarnings = FALSE)

# ─── Define comparisons ───────────────────────────────────────────────────────
# Single main-model contrast for the density experiment.
comparisons <- list(
  high_vs_low_ns = list(file = "high_vs_low_ns", source = "main")
)

# ─── Helper functions ─────────────────────────────────────────────────────────

# Load significant DE genes for a comparison.
load_sig_genes <- function(entry) {
  file_path <- paste0(sig_dir, "/annotated_sig_", entry$file, ".csv")
  if (!file.exists(file_path)) return(data.frame())
  read_csv(file_path, show_col_types = FALSE)
}

# Load raw GO enrichment terms for one comparison + direction.
load_go_terms <- function(entry, direction) {
  dir_path  <- if (direction == "UP") go_up_dir else go_down_dir
  file_path <- paste0(dir_path, "/go_annotated_", entry$file, ".csv")
  if (!file.exists(file_path)) return(data.frame())
  read_csv(file_path, show_col_types = FALSE)
}

# Load rrvgo-reduced GO terms for one comparison + direction.
load_go_reduced <- function(entry, direction) {
  dir_path  <- if (direction == "UP") go_up_red_dir else go_down_red_dir
  file_path <- paste0(dir_path, "/go_annotated_", entry$file, ".csv")
  if (!file.exists(file_path)) return(data.frame())
  df <- read_csv(file_path, show_col_types = FALSE)
  if (nrow(df) == 0) return(data.frame())
  df
}

# Format a GO term table block (raw OR reduced).
format_go_block <- function(go_df, reduced = FALSE) {
  if (nrow(go_df) == 0) return("  (No significant GO terms)")

  go_df <- go_df %>% arrange(p.adjust) %>% head(20)

  if (reduced) {
    header <- sprintf("  %-12s %-40s %10s %12s %12s %8s %8s",
                      "GO_ID", "GO Term", "P.adjust",
                      "GeneRatio", "BgRatio", "Enrich", "ClusterSz")
  } else {
    header <- sprintf("  %-12s %-40s %10s %12s %12s %8s",
                      "GO_ID", "GO Term", "P.adjust",
                      "GeneRatio", "BgRatio", "Enrichment")
  }
  divider <- paste0("  ", strrep("─", nchar(header) - 2))

  rows <- character(nrow(go_df))
  for (i in seq_len(nrow(go_df))) {
    go_id      <- go_df$ID[i]
    go_name    <- substr(go_df$GO_Name[i], 1, 38)
    gene_ratio <- go_df$GeneRatio[i]
    bg_ratio   <- go_df$BgRatio[i]
    gene_parts <- as.numeric(strsplit(gene_ratio, "/")[[1]])
    bg_parts   <- as.numeric(strsplit(bg_ratio,   "/")[[1]])
    fold_enr   <- (gene_parts[1] / gene_parts[2]) / (bg_parts[1] / bg_parts[2])

    if (reduced) {
      cluster_sz <- if ("cluster_size" %in% names(go_df)) go_df$cluster_size[i] else NA
      rows[i] <- sprintf("  %-12s %-40s %10.2e %12s %12s %8.2f %8s",
                         go_id, go_name, go_df$p.adjust[i],
                         gene_ratio, bg_ratio, fold_enr,
                         ifelse(is.na(cluster_sz), "?",
                                as.character(cluster_sz)))
    } else {
      rows[i] <- sprintf("  %-12s %-40s %10.2e %12s %12s %8.2f",
                         go_id, go_name, go_df$p.adjust[i],
                         gene_ratio, bg_ratio, fold_enr)
    }
  }

  c(header, divider, rows)
}

# ─── Load all data ────────────────────────────────────────────────────────────
cat("Loading comparison and annotations...\n\n")

all_genes       <- lapply(comparisons, load_sig_genes)
all_go_up       <- lapply(comparisons, function(e) load_go_terms(e,  "UP"))
all_go_down     <- lapply(comparisons, function(e) load_go_terms(e,  "DOWN"))
all_go_up_red   <- lapply(comparisons, function(e) load_go_reduced(e, "UP"))
all_go_down_red <- lapply(comparisons, function(e) load_go_reduced(e, "DOWN"))

cat("Data loaded successfully!\n\n")

#============================================================================
# INDIVIDUAL COMPARISON SUMMARY
#============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("ANALYZING THE HIGH vs LOW CONTRAST (NO-SHRINK)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (comp_name in names(comparisons)) {
  entry          <- comparisons[[comp_name]]
  comp_full_name <- entry$file
  genes_df       <- all_genes[[comp_name]]
  go_up          <- all_go_up[[comp_name]]
  go_down        <- all_go_down[[comp_name]]
  go_up_red      <- all_go_up_red[[comp_name]]
  go_down_red    <- all_go_down_red[[comp_name]]

  if (nrow(genes_df) == 0) {
    cat(sprintf("Skipping %s (no DE genes)\n", comp_name))
    next
  }

  cat(sprintf("Analyzing: %s\n", comp_name))

  genes_up   <- genes_df %>% filter(log2FoldChange > 0)
  genes_down <- genes_df %>% filter(log2FoldChange < 0)

  comparison_ref <- "HIGH vs LOW density (no shrinkage)"
  up_means       <- "UP = Higher in HIGH density | DOWN = Higher in LOW density"
  bio_question   <- "Which genes respond to rearing density (high vs low) in adult crickets?"
  interpretation <- "Genes here are differentially expressed between high-density and low-density treatments, using the raw (unshrunk) DESeq2 MLE log2 fold changes. UP genes = elevated under high density. DOWN genes = elevated under low density. This is the NO-SHRINK variant; compare against the apeglm-shrunk pipeline to assess how much the fold-change estimates and rankings shift."

  n_up_raw       <- nrow(go_up)
  n_down_raw     <- nrow(go_down)
  n_up_reduced   <- nrow(go_up_red)
  n_down_reduced <- nrow(go_down_red)

  summary_text <- c(
    "═══════════════════════════════════════════════════════════════",
    sprintf("COMPARISON: %s", comp_full_name),
    "═══════════════════════════════════════════════════════════════",
    "",
    "BIOLOGICAL QUESTION:",
    sprintf("  %s", bio_question),
    "",
    "INTERPRETATION GUIDE:",
    sprintf("  %s", interpretation),
    "",
    "CONTRAST DIRECTION:",
    sprintf("  %s", comparison_ref),
    sprintf("  %s", up_means),
    "",
    "OVERVIEW:",
    sprintf("  Total DE genes:  %d", nrow(genes_df)),
    sprintf("  UP-regulated:    %d genes", nrow(genes_up)),
    sprintf("  DOWN-regulated:  %d genes", nrow(genes_down)),
    sprintf("  GO terms (UP):   %d raw -> %d non-redundant (rrvgo)",
            n_up_raw, n_up_reduced),
    sprintf("  GO terms (DOWN): %d raw -> %d non-redundant (rrvgo)",
            n_down_raw, n_down_reduced),
    "",
    "═══════════════════════════════════════════════════════════════",
    "RAW GO TERMS",
    "═══════════════════════════════════════════════════════════════",
    "",
    "UP-REGULATED GO TERMS:",
    ""
  )

  summary_text <- c(summary_text, format_go_block(go_up, reduced = FALSE))
  summary_text <- c(summary_text, "", "DOWN-REGULATED GO TERMS:", "")
  summary_text <- c(summary_text, format_go_block(go_down, reduced = FALSE))

  summary_text <- c(
    summary_text,
    "",
    "═══════════════════════════════════════════════════════════════",
    "REDUCED GO TERMS (rrvgo, similarity threshold 0.8)",
    "═══════════════════════════════════════════════════════════════",
    "",
    "Each row below is a cluster representative — the most significant",
    "term of a group of semantically similar terms. The ClusterSz",
    "column shows how many raw redundant terms collapsed into each",
    "representative (ClusterSz = 1 means the term was already distinct).",
    "",
    "REDUCED UP-REGULATED GO TERMS:",
    ""
  )
  summary_text <- c(summary_text, format_go_block(go_up_red, reduced = TRUE))

  summary_text <- c(summary_text, "", "REDUCED DOWN-REGULATED GO TERMS:", "")
  summary_text <- c(summary_text, format_go_block(go_down_red, reduced = TRUE))

  output_file <- paste0(output_dir, "/01_individual_comparisons/",
                        comp_name, "_summary.txt")
  writeLines(summary_text, output_file)
  cat(sprintf("  Saved: %s\n", basename(output_file)))
}

cat("\n")

#============================================================================
# FINAL SUMMARY
#============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║          GO SYNTHESIS COMPLETE — density NO-SHRINK             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Output directory: deseq2_results_density_ns/comprehensive_go_analysis_density_ns/\n\n")
cat("  01_individual_comparisons/  Summary for the high-vs-low contrast (no shrink)\n")
cat("                              (raw + rrvgo-reduced GO term sections)\n\n")

cat("Comparison analysed:\n")
for (comp_name in names(comparisons)) {
  entry   <- comparisons[[comp_name]]
  genes_n <- nrow(all_genes[[comp_name]])
  cat(sprintf("  %-40s -> %d DE genes\n", entry$file, genes_n))
}
