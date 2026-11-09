# annotate_de_multifactorial_density.R
# Join eggNOG annotations with DESeq2 results from the NO-SHRINK density analysis


library(dplyr)
library(readr)
library(stringr)

cat("\n=== EggNOG Annotation Pipeline (density, NO-SHRINK) ===\n\n")

# ─── File paths ───────────────────────────────────────────────────────────────
eggnog_file <- "../eggnog/acheta_transdecoder_eggnog.emapper.annotations"
results_dir <- "deseq2_results_density_ns"
output_dir  <- "deseq2_results_density_ns/annotated"

dir.create(output_dir,                        showWarnings = FALSE)
dir.create(paste0(output_dir, "/full"),        showWarnings = FALSE)
dir.create(paste0(output_dir, "/significant"), showWarnings = FALSE)

# ─── Load eggNOG annotations ─────────────────────────────────────────────────
cat("Reading eggNOG annotations...\n")
eggnog <- read_tsv(eggnog_file,
                   comment    = "#",
                   col_names  = c("query", "seed_ortholog", "evalue", "score",
                                  "eggNOG_OGs", "max_annot_lvl", "COG_category",
                                  "Description", "Preferred_name", "GOs", "EC",
                                  "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
                                  "KEGG_Reaction", "KEGG_rclass", "BRITE",
                                  "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"),
                   skip       = 1,
                   show_col_types = FALSE)

# Remove TransDecoder .p suffix
eggnog <- eggnog %>%
  mutate(transcript_id = str_replace(query, "\\.p[0-9]+$", ""))

cat(sprintf("Loaded %d eggNOG annotations\n\n", nrow(eggnog)))

# ─── Find result files ────────────────────────────────────────────────────────
result_files <- list.files(results_dir, pattern = "^results_.*\\.csv$",
                           full.names = TRUE)

cat(sprintf("Found %d DESeq2 result files to annotate:\n", length(result_files)))
for (f in result_files) cat(sprintf("  - %s\n", basename(f)))
cat("\n")

# ─── Annotate function ────────────────────────────────────────────────────────
annotate_results <- function(result_file, eggnog_data, output_dir) {

  comparison <- basename(result_file) %>%
    str_replace("^results_", "") %>%
    str_replace("\\.csv$", "")

  cat(sprintf("Processing: %s\n", comparison))

  de <- read_csv(result_file, show_col_types = FALSE)

  de_annotated <- de %>%
    left_join(eggnog_data, by = c("name" = "transcript_id"))

  # Column ordering: DESeq2 stats → key annotation → everything else
  deseq_cols <- intersect(c("name", "baseMean", "log2FoldChange", "lfcSE",
                             "stat", "PValue", "PAdj", "FDR", "falsePos"),
                          names(de_annotated))
  annot_cols <- intersect(c("Preferred_name", "Description", "GOs", "KEGG_ko",
                             "KEGG_Pathway", "COG_category", "EC", "PFAMs"),
                          names(de_annotated))
  other_cols <- setdiff(names(de_annotated), c(deseq_cols, annot_cols))

  de_annotated <- de_annotated %>%
    select(all_of(c(deseq_cols, annot_cols, other_cols)))

  # Save full annotated results
  write_csv(de_annotated,
            paste0(output_dir, "/full/annotated_full_", comparison, ".csv"))

  # Filter significant: FDR < 0.05 and |log2FC| > 1
  de_sig <- de_annotated %>%
    filter(FDR < 0.05, abs(log2FoldChange) > 1)

  write_csv(de_sig,
            paste0(output_dir, "/significant/annotated_sig_", comparison, ".csv"))

  annotated_count  <- sum(!is.na(de_annotated$Preferred_name))
  annotation_rate  <- round(100 * annotated_count / nrow(de_annotated), 1)

  cat(sprintf("  Total genes:                          %d\n",  nrow(de)))
  cat(sprintf("  Annotated:                            %d (%.1f%%)\n",
              annotated_count, annotation_rate))
  cat(sprintf("  Significant (FDR<0.05, |log2FC|>1):  %d\n",  nrow(de_sig)))

  if (nrow(de_sig) > 0) {
    sig_annotated <- sum(!is.na(de_sig$Preferred_name))
    cat(sprintf("  Significant + annotated:              %d (%.1f%%)\n",
                sig_annotated,
                round(100 * sig_annotated / nrow(de_sig), 1)))
  }
  cat("\n")

  return(list(
    comparison  = comparison,
    total       = nrow(de),
    annotated   = annotated_count,
    significant = nrow(de_sig)
  ))
}

# ─── Run annotation ───────────────────────────────────────────────────────────
cat("=== Annotating all result files ===\n\n")
summary_list <- lapply(result_files, function(f) {
  annotate_results(f, eggnog, output_dir)
})

# ─── Summary table ────────────────────────────────────────────────────────────
cat("\n=== ANNOTATION SUMMARY ===\n\n")
summary_df <- do.call(rbind, lapply(summary_list, function(x) {
  data.frame(
    Comparison  = x$comparison,
    Total_Genes = x$total,
    Annotated   = x$annotated,
    Significant_DE = x$significant,
    stringsAsFactors = FALSE
  )
}))

print(summary_df, row.names = FALSE)
write_csv(summary_df, paste0(output_dir, "/annotation_summary.csv"))

cat("\n=== COMPLETE ===\n")
cat(sprintf("All annotated results saved to: %s/\n\n", output_dir))
cat("  - full/:       All genes with annotations\n")
cat("  - significant/ Only significant genes (FDR<0.05, |log2FC|>1)\n")
cat("  - annotation_summary.csv\n\n")
cat("Filtering criteria:\n")
cat("  FDR < 0.05    | log2FC > 1 → 2-fold UP\n")
cat("  |log2FC| > 1  | log2FC < -1 → 2-fold DOWN\n\n")
