#!/usr/bin/env Rscript
# rrvgo_reduce_density_.R


suppressPackageStartupMessages({
  library(rrvgo)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(GO.db)
  library(org.Dm.eg.db)
})

setwd("~/work/space/code")

# ---- config ----
THRESHOLD <- 0.8
PROXY_ORG <- "org.Dm.eg.db"
VERBOSE   <- TRUE

DIRS <- list(
  list("deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_UP_annotated",
       "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_UP_reduced",
       "comprehensive_UP"),
  list("deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_DOWN_annotated",
       "deseq2_results_density_ns/comprehensive_go_analysis_density_ns/GO_DOWN_reduced",
       "comprehensive_DOWN")
)

SUMMARY_PATH <- "deseq2_results_density_ns/rrvgo_reduction_summary.csv"

# ---- helpers ----

msg <- function(...) if (VERBOSE) cat("    ", sprintf(...), "\n", sep = "")

# Look up the ontology (BP / MF / CC) for each GO ID.
get_ontology <- function(go_ids) {
  if (length(go_ids) == 0) return(character(0))

  res <- tryCatch(
    AnnotationDbi::select(GO.db, keys = go_ids,
                          columns = "ONTOLOGY", keytype = "GOID"),
    error = function(e) NULL
  )

  if (is.null(res) || nrow(res) == 0) {
    return(rep(NA_character_, length(go_ids)))
  }

  res <- res %>% dplyr::distinct(GOID, .keep_all = TRUE)
  res$ONTOLOGY[match(go_ids, res$GOID)]
}

# Every input term becomes its own "cluster". Safe pass-through.
pass_through <- function(df) {
  df %>% dplyr::mutate(parentTerm = GO_Name, cluster_size = 1L)
}

# Reduce one ontology slice. Guarantees every input term appears in the
# output as either a representative or a singleton.
reduce_one_ontology <- function(df, ontology) {
  msg("[%s] %d input term(s)", ontology, nrow(df))

  if (nrow(df) < 2) {
    msg("only %d term — pass-through", nrow(df))
    return(pass_through(df))
  }

  sim <- tryCatch(
    calculateSimMatrix(df$ID, orgdb = PROXY_ORG, ont = ontology, method = "Rel"),
    error = function(e) {
      msg("calculateSimMatrix error: %s", conditionMessage(e))
      NULL
    }
  )

  if (is.null(sim) || !is.matrix(sim) || nrow(sim) < 2) {
    msg("similarity matrix unusable (rows=%s) — pass-through",
        if (is.null(sim)) "NULL" else as.character(nrow(sim)))
    return(pass_through(df))
  }

  msg("similarity matrix: %d x %d", nrow(sim), ncol(sim))

  sim_ids     <- rownames(sim)
  dropped_ids <- setdiff(df$ID, sim_ids)
  if (length(dropped_ids) > 0) {
    msg("%d term(s) absent from Dm orthology — kept as singletons: %s",
        length(dropped_ids),
        paste(dropped_ids, collapse = ", "))
  }

  scores <- setNames(-log10(pmax(df$p.adjust, 1e-300)), df$ID)
  scores <- scores[sim_ids]

  reduced <- tryCatch(
    reduceSimMatrix(sim, scores, threshold = THRESHOLD, orgdb = PROXY_ORG),
    error = function(e) {
      msg("reduceSimMatrix error: %s", conditionMessage(e))
      NULL
    }
  )

  if (is.null(reduced) || !is.data.frame(reduced) || nrow(reduced) == 0) {
    msg("reduceSimMatrix returned empty — pass-through")
    return(pass_through(df))
  }

  msg("reduceSimMatrix returned %d row(s)", nrow(reduced))
  msg("returned columns: %s", paste(names(reduced), collapse = ", "))

  # Defensive column detection — rrvgo's API has changed names across versions
  go_col       <- intersect(c("go", "term", "ID", "GO"),            names(reduced))[1]
  par_col      <- intersect(c("parent", "cluster"),                 names(reduced))[1]
  par_name_col <- intersect(c("parentTerm", "term_name", "Term"),   names(reduced))[1]
  size_col     <- intersect(c("size", "cluster_size", "termSize"),  names(reduced))[1]

  if (is.na(go_col) || is.na(par_col)) {
    msg("required columns (go, parent) not found — pass-through")
    return(pass_through(df))
  }

  msg("using columns: go=%s, parent=%s, parentTerm=%s, size=%s",
      go_col, par_col,
      if (is.na(par_name_col)) "(fallback to GO_Name)" else par_name_col,
      if (is.na(size_col))     "(computed from table)" else size_col)

  reduced[[go_col]]  <- as.character(reduced[[go_col]])
  reduced[[par_col]] <- as.character(reduced[[par_col]])

  # Representatives: rows where the term IS its own parent.
  is_rep   <- reduced[[go_col]] == reduced[[par_col]]
  reps_raw <- reduced[is_rep, , drop = FALSE]
  msg("rows where go == parent: %d", nrow(reps_raw))

  if (nrow(reps_raw) == 0) {
    msg("no representatives found — pass-through")
    return(pass_through(df))
  }

  # Compute cluster sizes if not provided
  cluster_sizes <- if (!is.na(size_col)) {
    setNames(reps_raw[[size_col]], reps_raw[[go_col]])
  } else {
    as.integer(table(reduced[[par_col]]))[reps_raw[[go_col]]] %>%
      setNames(reps_raw[[go_col]])
  }

  # Use parentTerm column if present, otherwise we'll fall back to GO_Name from df
  parent_terms <- if (!is.na(par_name_col)) {
    as.character(reps_raw[[par_name_col]])
  } else {
    rep(NA_character_, nrow(reps_raw))
  }

  reps <- tibble(
    ID                = reps_raw[[go_col]],
    parentTerm_rrvgo  = parent_terms,
    cluster_size      = as.integer(cluster_sizes[reps_raw[[go_col]]])
  )

  msg("representatives identified: %d (compressed from %d clustered terms)",
      nrow(reps), length(sim_ids))

  # Join representatives back to df, fall back to GO_Name if rrvgo didn't supply parentTerm
  clustered <- df %>%
    dplyr::inner_join(reps, by = "ID") %>%
    dplyr::mutate(parentTerm = dplyr::coalesce(parentTerm_rrvgo, GO_Name)) %>%
    dplyr::select(-parentTerm_rrvgo)

  # Pass-through anything calculateSimMatrix didn't see
  passed_through <- df %>%
    dplyr::filter(ID %in% dropped_ids) %>%
    pass_through()

  result <- bind_rows(clustered, passed_through) %>% dplyr::arrange(p.adjust)

  msg("final: %d term(s) in output (%d cluster reps + %d singletons)",
      nrow(result), nrow(clustered), nrow(passed_through))

  result
}

# Process one input file end-to-end.
process_file <- function(in_path, out_path) {
  cat(sprintf("\n%s\n", basename(in_path)))

  df <- read_csv(in_path, show_col_types = FALSE)

  needed <- c("ID", "GO_Name", "p.adjust")
  if (!all(needed %in% names(df))) {
    warning(sprintf("Skipping %s: missing columns (%s)",
                    basename(in_path),
                    paste(setdiff(needed, names(df)), collapse = ", ")))
    return(tibble(file = basename(in_path),
                  n_in = nrow(df), n_out = NA_integer_))
  }

  if (nrow(df) == 0) {
    write_csv(df, out_path)
    msg("empty input — empty output")
    return(tibble(file = basename(in_path), n_in = 0L, n_out = 0L))
  }

  df$ID       <- as.character(df$ID)
  df$ontology <- get_ontology(df$ID)

  n_no_ont <- sum(is.na(df$ontology))
  if (n_no_ont > 0) {
    msg("%d term(s) have no ontology assignment — kept as singletons", n_no_ont)
  }

  ont_split <- df %>%
    dplyr::filter(!is.na(ontology)) %>%
    split(.$ontology)

  reduced_parts <- imap(ont_split, ~ reduce_one_ontology(.x, .y))

  no_ont <- df %>%
    dplyr::filter(is.na(ontology)) %>%
    pass_through()

  reduced_df <- bind_rows(c(reduced_parts, list(no_ontology = no_ont))) %>%
    dplyr::arrange(p.adjust)

  write_csv(reduced_df, out_path)

  cat(sprintf("    %d in -> %d out\n", nrow(df), nrow(reduced_df)))

  tibble(
    file  = basename(in_path),
    n_in  = nrow(df),
    n_out = nrow(reduced_df)
  )
}

# ---- run ----

summary_rows <- list()

for (cfg in DIRS) {
  in_dir  <- cfg[[1]]
  out_dir <- cfg[[2]]
  label   <- cfg[[3]]

  if (!dir.exists(in_dir)) {
    cat(sprintf("\n[skip] %s does not exist\n", in_dir))
    next
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
  cat(sprintf("\n=== [%s] %d file(s) in %s ===\n", label, length(files), in_dir))

  for (f in files) {
    out_path <- file.path(out_dir, basename(f))
    res <- tryCatch(
      process_file(f, out_path),
      error = function(e) {
        warning(sprintf("Error processing %s: %s",
                        basename(f), conditionMessage(e)))
        tibble(file = basename(f),
               n_in = NA_integer_, n_out = NA_integer_)
      }
    )
    res$source <- label
    summary_rows[[length(summary_rows) + 1]] <- res
  }
}

if (length(summary_rows) == 0) {
  cat("\n\n[done] No GO enrichment files were found to reduce.\n")
  cat("       This is expected when the enrichment step produced no GO_UP/\n")
  cat("       GO_DOWN terms (e.g. too few significant annotated genes at this\n")
  cat("       sample size). Nothing to reduce — writing an empty summary.\n\n")
  empty_summary <- tibble(source = character(0), file = character(0),
                          n_in = integer(0), n_out = integer(0),
                          reduction_pct = numeric(0))
  write_csv(empty_summary, SUMMARY_PATH)
  cat(sprintf("Empty summary written to %s\n\n", SUMMARY_PATH))
} else {
  summary_df <- bind_rows(summary_rows) %>%
    dplyr::select(source, file, n_in, n_out) %>%
    dplyr::mutate(reduction_pct = round(100 * (1 - n_out / pmax(n_in, 1)), 1))

  write_csv(summary_df, SUMMARY_PATH)
  cat(sprintf("\n\nSummary written to %s\n\n", SUMMARY_PATH))
  print(summary_df, n = Inf)
}
