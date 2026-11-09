# enrich_multifactorial_density.R
# Directional GO and KEGG enrichment for the NO-SHRINK density DESeq2 results
# Run from ~/work/space/code

library(clusterProfiler)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

cat("\n=== Directional Functional Enrichment Analysis (density, NO-SHRINK) ===\n\n")

# ─── Directories ──────────────────────────────────────────────────────────────
sig_dir    <- "deseq2_results_density_ns/annotated/significant"
full_dir   <- "deseq2_results_density_ns/annotated/full"
enrich_dir <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns"

dir.create(enrich_dir,                           showWarnings = FALSE)
dir.create(paste0(enrich_dir, "/GO_UP"),         showWarnings = FALSE)
dir.create(paste0(enrich_dir, "/GO_DOWN"),       showWarnings = FALSE)
dir.create(paste0(enrich_dir, "/KEGG_UP"),       showWarnings = FALSE)
dir.create(paste0(enrich_dir, "/KEGG_DOWN"),     showWarnings = FALSE)
dir.create(paste0(enrich_dir, "/GO_UP_annotated"),   showWarnings = FALSE)
dir.create(paste0(enrich_dir, "/GO_DOWN_annotated"), showWarnings = FALSE)

# ─── Build background universe ────────────────────────────────────────────────
cat("Building background gene universe...\n")
full_file    <- list.files(full_dir, pattern = "^annotated_full_.*\\.csv$",
                           full.names = TRUE)[1]
de_annotated <- read_csv(full_file, show_col_types = FALSE)

background <- de_annotated$name[
  !is.na(de_annotated$Preferred_name) |
  !is.na(de_annotated$GOs) |
  !is.na(de_annotated$KEGG_ko)
]
cat(sprintf("Background universe: %d genes with annotations\n\n", length(background)))

# ─── Prepare GO and KEGG term mappings ────────────────────────────────────────
cat("Preparing GO term mappings...\n")
go_data <- de_annotated %>%
  filter(!is.na(GOs)) %>%
  select(name, GOs) %>%
  mutate(GOs = str_split(GOs, ",")) %>%
  unnest(GOs) %>%
  mutate(GOs = str_trim(GOs))

cat(sprintf("GO annotations: %d unique GO terms across %d genes\n",
            length(unique(go_data$GOs)), length(unique(go_data$name))))

cat("Preparing KEGG KO mappings...\n")
kegg_data <- de_annotated %>%
  filter(!is.na(KEGG_ko)) %>%
  select(name, KEGG_ko) %>%
  mutate(KEGG_ko = str_split(KEGG_ko, ",")) %>%
  unnest(KEGG_ko) %>%
  mutate(KEGG_ko = str_trim(KEGG_ko))

cat(sprintf("KEGG annotations: %d unique KO terms across %d genes\n\n",
            length(unique(kegg_data$KEGG_ko)), length(unique(kegg_data$name))))

# ─── Helper: annotate GO result with term names ───────────────────────────────
annotate_go_result <- function(go_result_df, go_data_full) {
  if ("Description" %in% names(go_result_df)) {
    go_result_df <- go_result_df %>%
      rename(GO_Name = Description)
  }
  return(go_result_df)
}

# ─── Directional enrichment function ─────────────────────────────────────────
directional_enrichment <- function(sig_file, background, go_data, kegg_data,
                                   enrich_dir) {

  comparison <- basename(sig_file) %>%
    str_replace("^annotated_sig_", "") %>%
    str_replace("\\.csv$", "")

  cat(sprintf("\n=== Processing: %s ===\n", comparison))

  sig_genes_df <- read_csv(sig_file, show_col_types = FALSE)
  up_genes     <- sig_genes_df %>% filter(log2FoldChange > 0) %>% pull(name)
  down_genes   <- sig_genes_df %>% filter(log2FoldChange < 0) %>% pull(name)

  cat(sprintf("Total DE genes: %d  (UP: %d | DOWN: %d)\n",
              nrow(sig_genes_df), length(up_genes), length(down_genes)))

  up_go   <- sum(up_genes   %in% go_data$name)
  down_go <- sum(down_genes %in% go_data$name)
  up_kegg <- sum(up_genes   %in% kegg_data$name)
  down_kegg <- sum(down_genes %in% kegg_data$name)

  cat(sprintf("  UP   with GO: %d (%.0f%%)  |  DOWN with GO: %d (%.0f%%)\n",
              up_go,   100 * up_go   / max(length(up_genes),   1),
              down_go, 100 * down_go / max(length(down_genes), 1)))
  cat(sprintf("  UP with KEGG: %d (%.0f%%)  |  DOWN with KEGG: %d (%.0f%%)\n",
              up_kegg,   100 * up_kegg   / max(length(up_genes),   1),
              down_kegg, 100 * down_kegg / max(length(down_genes), 1)))

  results <- list(comparison = comparison,
                  n_up = length(up_genes), n_down = length(down_genes),
                  n_go_up = 0, n_go_down = 0,
                  n_kegg_up = 0, n_kegg_down = 0)

# === GO UP ===
  if (up_go >= 5) {
    cat("  Running GO enrichment (UP)...\n")
    tryCatch({
      go_up <- enricher(gene = up_genes, 
                        universe = background,
                        TERM2GENE = go_data %>% select(GOs, name),
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.2,
                        minGSSize = 5,     # ← This is the key parameter
                        maxGSSize = 500)   # prevents huge overly general terms

      if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
        df <- as.data.frame(go_up)
        df <- annotate_go_result(df, go_data)
        
        # Extra safety filter (in case any slipped through)
        df <- df %>% filter(Count >= 5)
        
        write_csv(df, paste0(enrich_dir, "/GO_UP/go_UP_", comparison, ".csv"))
        write_csv(df, paste0(enrich_dir, "/GO_UP_annotated/go_annotated_", comparison, ".csv"))
        
        results$n_go_up <- nrow(df)
        cat(sprintf("    GO terms enriched (UP): %d\n", nrow(df)))
      } else { 
        cat("    No significant GO enrichment (UP)\n") 
      }
    }, error = function(e) cat(sprintf("    GO (UP) failed: %s\n", e$message)))
  } else { 
    cat("  Skipping GO (UP): too few annotated genes (< 5)\n") 
  }

# === GO DOWN ===
  if (down_go >= 5) {
    cat("  Running GO enrichment (DOWN)...\n")
    tryCatch({
      go_down <- enricher(gene = down_genes, 
                          universe = background,
                          TERM2GENE = go_data %>% select(GOs, name),
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.2,
                          minGSSize = 5,
                          maxGSSize = 500)

      if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
        df <- as.data.frame(go_down)
        df <- annotate_go_result(df, go_data)
        df <- df %>% filter(Count >= 5)
        
        write_csv(df, paste0(enrich_dir, "/GO_DOWN/go_DOWN_", comparison, ".csv"))
        write_csv(df, paste0(enrich_dir, "/GO_DOWN_annotated/go_annotated_", comparison, ".csv"))
        
        results$n_go_down <- nrow(df)
        cat(sprintf("    GO terms enriched (DOWN): %d\n", nrow(df)))
      } else { 
        cat("    No significant GO enrichment (DOWN)\n") 
      }
    }, error = function(e) cat(sprintf("    GO (DOWN) failed: %s\n", e$message)))
  } else { 
    cat("  Skipping GO (DOWN): too few annotated genes (< 5)\n") 
  }

# === KEGG ===
  if (up_kegg >= 3) {
    cat("  Running KEGG enrichment (UP)...\n")
    tryCatch({
      kegg_up <- enricher(gene = up_genes, universe = background,
                          TERM2GENE = kegg_data %>% select(KEGG_ko, name),
                          pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                          minGSSize = 3, maxGSSize = 300)
      if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
        write_csv(as.data.frame(kegg_up),
                  paste0(enrich_dir, "/KEGG_UP/kegg_UP_", comparison, ".csv"))
        results$n_kegg_up <- nrow(as.data.frame(kegg_up))
        cat(sprintf("    KEGG terms enriched (UP): %d\n", results$n_kegg_up))
      } else { cat("    No significant KEGG enrichment (UP)\n") }
    }, error = function(e) cat(sprintf("    KEGG (UP) failed: %s\n", e$message)))
  } else { cat("  Skipping KEGG (UP): too few genes\n") }

  if (down_kegg >= 3) {
    cat("  Running KEGG enrichment (DOWN)...\n")
    tryCatch({
      kegg_down <- enricher(gene = down_genes, universe = background,
                            TERM2GENE = kegg_data %>% select(KEGG_ko, name),
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                            minGSSize = 3, maxGSSize = 300)
      if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
        write_csv(as.data.frame(kegg_down),
                  paste0(enrich_dir, "/KEGG_DOWN/kegg_DOWN_", comparison, ".csv"))
        results$n_kegg_down <- nrow(as.data.frame(kegg_down))
        cat(sprintf("    KEGG terms enriched (DOWN): %d\n", results$n_kegg_down))
      } else { cat("    No significant KEGG enrichment (DOWN)\n") }
    }, error = function(e) cat(sprintf("    KEGG (DOWN) failed: %s\n", e$message)))
  } else { cat("  Skipping KEGG (DOWN): too few genes\n") }

  return(results)
}

# ─── Run enrichment ───────────────────────────────────────────────────────────
sig_files <- list.files(sig_dir, pattern = "^annotated_sig_.*\\.csv$",
                        full.names = TRUE)

cat(sprintf("Found %d significant gene files to process\n", length(sig_files)))

cat("\n=== Running Directional Enrichment Analysis ===\n")
results_list <- lapply(sig_files, function(f) {
  directional_enrichment(f, background, go_data, kegg_data, enrich_dir)
})

# ─── Summary table ────────────────────────────────────────────────────────────
cat("\n\n=== DIRECTIONAL ENRICHMENT SUMMARY ===\n\n")
summary_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(
    Comparison = x$comparison,
    UP_genes   = x$n_up,
    DOWN_genes = x$n_down,
    GO_UP      = x$n_go_up,
    GO_DOWN    = x$n_go_down,
    KEGG_UP    = x$n_kegg_up,
    KEGG_DOWN  = x$n_kegg_down,
    stringsAsFactors = FALSE
  )
}))

print(summary_df, row.names = FALSE)
write_csv(summary_df,
          paste0(enrich_dir, "/directional_enrichment_summary.csv"))

cat("\n=== COMPLETE ===\n\n")
cat(sprintf("Enrichment results saved to: %s/\n\n", enrich_dir))
cat("Directory structure:\n")
cat("  GO_UP/            GO terms enriched in UP-regulated genes\n")
cat("  GO_DOWN/          GO terms enriched in DOWN-regulated genes\n")
cat("  GO_UP_annotated/  Same as GO_UP but with term name column\n")
cat("  GO_DOWN_annotated/ Same as GO_DOWN with term name column\n")
cat("  KEGG_UP/          KEGG pathways enriched in UP-regulated genes\n")
cat("  KEGG_DOWN/        KEGG pathways enriched in DOWN-regulated genes\n")
cat("  directional_enrichment_summary.csv\n\n")
