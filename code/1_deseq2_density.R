# deseq2_multifactorial_density_ns.R
# Differential expression analysis for density experiment (high vs low density).
# NO-SHRINK variant: log2FoldChange is the raw DESeq2 MLE estimate (no apeglm
# lfcShrink). Use this alongside the shrunk pipeline to compare how much the
# LFC estimates and rankings change without shrinkage.
# Single factor: condition (low = reference, high vs low contrast).

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))

# ─── File paths ───────────────────────────────────────────────────────────────
counts_file  <- "../counts.csv"
design_file  <- "design.csv"
output_dir   <- "deseq2_results_density_ns"

dir.create(output_dir, showWarnings = FALSE)

# ─── Load data ────────────────────────────────────────────────────────────────
colData <- read.csv(design_file, stringsAsFactors = FALSE)

# Single experimental factor: condition (low vs high density).
# low is the reference level -> the named coefficient is condition_high_vs_low,
# so a positive log2FoldChange means UP in high density.
colData$condition <- relevel(factor(colData$condition), ref = "low")

# Read counts
df           <- read.csv(counts_file, header = TRUE, row.names = 1)
sample_names <- colData$sample
countData    <- round(df[, sample_names])
otherCols    <- df[!(names(df) %in% sample_names)]

cat("\n=== EXPERIMENTAL DESIGN SUMMARY (NO-SHRINK) ===\n")
cat("Total samples:", nrow(colData), "\n")
cat("Conditions:", paste(levels(colData$condition), collapse = ", "),
    "(reference: low)\n\n")
cat("Condition counts:\n")
print(table(colData$condition))
cat("\n")

# ─── Helper: format DESeq2 results ────────────────────────────────────────────
format_results <- function(res, normed, otherCols) {
  res_df    <- data.frame(res)
  ids       <- rownames(res_df)
  otherCols <- otherCols[ids, , drop = FALSE]
  data <- cbind(otherCols, res_df)
  data$foldChange <- 2 ^ data$log2FoldChange
  names(data)[names(data) == "pvalue"] <- "PValue"
  names(data)[names(data) == "padj"]   <- "FDR"
  data$PAdj     <- p.adjust(data$PValue, method = "hochberg")
  data          <- data[with(data, order(PValue, -foldChange)), ]
  data$falsePos <- 1:nrow(data) * data$FDR
  total         <- merge(data, normed, by = 0)
  colnames(total)[1] <- "name"
  total         <- total[with(total, order(PValue, -foldChange)), ]
  total$foldChange      <- round(total$foldChange, 3)
  total$log2FoldChange  <- round(total$log2FoldChange, 1)
  total$baseMean        <- round(total$baseMean, 1)
  total$lfcSE           <- round(total$lfcSE, 2)
  total$stat            <- round(total$stat, 2)
  total$FDR             <- round(total$FDR, 4)
  total$falsePos        <- round(total$falsePos, 0)
  total$PAdj  <- formatC(total$PAdj,  format = "e", digits = 1)
  total$PValue <- formatC(total$PValue, format = "e", digits = 1)
  return(total)
}

# ─── Low-count pre-filter ─────────────────────────────────────────────────────
filter_low_counts <- function(dds, min_count = 10, min_samples = NULL) {
  if (is.null(min_samples)) {
    min_samples <- min(table(colData(dds)$condition))
  }
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  cat(sprintf("  Low-count filter: kept %d of %d transcripts (>=%d counts in >=%d samples)\n",
              sum(keep), length(keep), min_count, min_samples))
  dds[keep, ]
}

summarize_results <- function(results_df, comparison_name, alpha = 0.05) {
  fdr_num  <- as.numeric(results_df$FDR)
  lfc      <- results_df$log2FoldChange
  sig      <- sum(fdr_num < alpha, na.rm = TRUE)
  up       <- sum(fdr_num < alpha & lfc > 0, na.rm = TRUE)
  down     <- sum(fdr_num < alpha & lfc < 0, na.rm = TRUE)
  cat(sprintf("  %-55s: %5d DE genes [%d up, %d down]\n",
              comparison_name, sig, up, down))
}

# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS: HIGH vs LOW DENSITY (NO SHRINKAGE)
# ─────────────────────────────────────────────────────────────────────────────
cat("=== ANALYSIS: High vs Low density (NO-SHRINK) ===\n")
cat("Design: ~ condition\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)
dds <- filter_low_counts(dds)

cat("Background after filter:", nrow(dds), "transcripts\n\n")

dds <- DESeq(dds)

normed <- round(counts(dds, normalized = TRUE), 1)

cat("Model coefficients available:\n")
print(resultsNames(dds))
cat("\n")

# === NO shrinkage: keep the raw DESeq2 MLE log2FoldChange ===
# This is the only logical difference from the shrunk pipeline. results() with
# the named coefficient returns the unshrunk Wald-test LFC and matching FDR; we
# write it as-is, with no lfcShrink/apeglm step.
cat("Extracting: High vs Low density effect (raw MLE, no shrinkage)...\n")
res <- results(dds, name = "condition_high_vs_low")
r   <- format_results(res, normed, otherCols)
write.csv(r, paste0(output_dir, "/results_high_vs_low_ns.csv"),
          row.names = FALSE, quote = FALSE)

cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
cat("=== SUMMARY OF ALL RESULTS (NO-SHRINK) ===\n\n")
summarize_results(read.csv(paste0(output_dir, "/results_high_vs_low_ns.csv")),
                  "High vs Low density (no shrink)")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(paste("All results saved to:", output_dir, "\n"))
cat("\nFiles generated:\n")
list.files(output_dir, pattern = "*.csv", full.names = FALSE) %>%
  paste0("  - ", .) %>%
  cat(sep = "\n")
