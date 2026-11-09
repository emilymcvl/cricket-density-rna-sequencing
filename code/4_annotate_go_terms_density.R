# annotate_go_terms_density.R - Add GO term descriptions to NO-SHRINK enrichment results
# Run from ~/work/space/code

library(GO.db)
library(dplyr)
library(readr)
library(stringr)

cat("\n=== GO Term Annotation Pipeline (density, NO-SHRINK) ===\n\n")

# Directories
go_up_dir       <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_UP"
go_down_dir     <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_DOWN"
output_up_dir   <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_UP_annotated"
output_down_dir <- "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_DOWN_annotated"

dir.create(output_up_dir,   showWarnings = FALSE)
dir.create(output_down_dir, showWarnings = FALSE)

# Get all GO enrichment files
go_up_files   <- list.files(go_up_dir,   pattern = "^go_UP_.*\\.csv$",   full.names = TRUE)
go_down_files <- list.files(go_down_dir, pattern = "^go_DOWN_.*\\.csv$", full.names = TRUE)

cat(sprintf("Found %d GO UP files and %d GO DOWN files to annotate\n\n",
            length(go_up_files), length(go_down_files)))

# Function to get GO term information from GO.db
get_go_info <- function(go_id) {
  tryCatch({
    term_name  <- Term(GOTERM[[go_id]])
    definition <- Definition(GOTERM[[go_id]])
    ontology   <- Ontology(GOTERM[[go_id]])
    ontology_full <- switch(ontology,
                            "BP" = "Biological Process",
                            "MF" = "Molecular Function",
                            "CC" = "Cellular Component",
                            ontology)
    return(data.frame(GO_ID = go_id, GO_Name = term_name,
                      GO_Category = ontology_full, GO_Definition = definition,
                      stringsAsFactors = FALSE))
  }, error = function(e) {
    return(data.frame(GO_ID = go_id, GO_Name = NA,
                      GO_Category = NA, GO_Definition = NA,
                      stringsAsFactors = FALSE))
  })
}

# Function to annotate one GO enrichment file
annotate_go_file <- function(go_file, output_dir) {

  comparison <- basename(go_file) %>%
    str_replace("^go_UP_",   "") %>%
    str_replace("^go_DOWN_", "") %>%
    str_replace("\\.csv$",   "")

  cat(sprintf("Processing: %s\n", comparison))

  go_enrichment <- read_csv(go_file, show_col_types = FALSE)

  if (nrow(go_enrichment) == 0) {
    cat("  No enriched terms - skipping\n\n")
    return(NULL)
  }

  go_ids    <- go_enrichment$ID
  cat(sprintf("  Enriched GO terms: %d\n", length(go_ids)))

  go_info_list <- lapply(go_ids, get_go_info)
  go_info      <- do.call(rbind, go_info_list)

  go_annotated <- go_enrichment %>%
    select(-any_of("GO_Name")) %>%
    left_join(go_info, by = c("ID" = "GO_ID"))

  col_order    <- c("ID", "GO_Name", "GO_Category", "GO_Definition",
                    setdiff(names(go_annotated),
                            c("ID", "GO_Name", "GO_Category", "GO_Definition")))
  go_annotated <- go_annotated[, col_order]

  category_counts <- table(go_annotated$GO_Category)
  cat("  Categories:\n")
  for (cat_name in names(category_counts)) {
    cat(sprintf("    - %s: %d\n", cat_name, category_counts[cat_name]))
  }

  output_file <- paste0(output_dir, "/go_annotated_", comparison, ".csv")
  write_csv(go_annotated, output_file)
  cat(sprintf("  Saved: %s\n\n", basename(output_file)))

  return(data.frame(
    Comparison  = comparison,
    Total_Terms = nrow(go_annotated),
    BP = sum(go_annotated$GO_Category == "Biological Process", na.rm = TRUE),
    MF = sum(go_annotated$GO_Category == "Molecular Function", na.rm = TRUE),
    CC = sum(go_annotated$GO_Category == "Cellular Component", na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

# Process UP files
cat("=== Annotating GO UP Terms ===\n\n")
summary_up_list <- lapply(go_up_files, function(f) {
  result <- annotate_go_file(f, output_up_dir)
  if (!is.null(result)) result$Direction <- "UP"
  return(result)
})

# Process DOWN files
cat("\n=== Annotating GO DOWN Terms ===\n\n")
summary_down_list <- lapply(go_down_files, function(f) {
  result <- annotate_go_file(f, output_down_dir)
  if (!is.null(result)) result$Direction <- "DOWN"
  return(result)
})

# Summary
summary_list <- c(summary_up_list, summary_down_list)
summary_list <- summary_list[!sapply(summary_list, is.null)]

if (length(summary_list) > 0) {
  cat("\n=== GO ANNOTATION SUMMARY ===\n\n")
  summary_df <- do.call(rbind, summary_list)
  print(summary_df, row.names = FALSE)
  write_csv(summary_df,
            "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/go_annotation_summary.csv")
}

cat("\n=== COMPLETE ===\n")
cat(sprintf("Annotated GO files saved to:\n"))
cat(sprintf("  UP:   %s/\n", output_up_dir))
cat(sprintf("  DOWN: %s/\n", output_down_dir))
cat("\nEach file now includes:\n")
cat("  GO_Name:       Human-readable term name\n")
cat("  GO_Category:   Biological Process / Molecular Function / Cellular Component\n")
cat("  GO_Definition: Full description of what the term means\n\n")
