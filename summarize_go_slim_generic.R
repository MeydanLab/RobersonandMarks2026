#!/usr/bin/env Rscript

# Summarize clusterProfiler::gseGO CSV output with a GO slim subset.
#
# How to run:
#   1. Run run_rnaseq_de_go_rif_analysis.R first to generate one or more
#      Tables/GSEA_*.csv files.
#   2. Pass those GSEA CSV files to this script as a comma-separated list.
#   3. Provide matching comma-separated labels in the same order.
#
# Example:
#   Rscript summarize_go_slim_generic.R \
#     --gsea-files "Tables/GSEA_comparison_a_BP.csv,Tables/GSEA_comparison_b_BP.csv" \
#     --labels "comparison A,comparison B" \
#     --outdir "GO_slim_generic" \
#     --prefix GO_BP \
#     --cache-dir "data/go_slim_cache" \
#     --slim-name goslim_generic \
#     --go-basic-source "https://current.geneontology.org/ontology/go-basic.obo" \
#     --goslim-source "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"
#
# Outputs:
#   - *_<slim-name>_expanded_*.csv: each original GO term mapped to its
#     nearest GO slim ancestor.
#   - *_<slim-name>_summary_*.csv: slim-level summary table.

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript summarize_go_slim_generic.R --gsea-files <file1.csv,file2.csv,...> --labels <label1,label2,...> --outdir <dir> --prefix <name>",
      "[--cache-dir data/go_slim_cache]",
      "[--slim-name goslim_generic]",
      "[--go-basic-source https://current.geneontology.org/ontology/go-basic.obo]",
      "[--goslim-source https://current.geneontology.org/ontology/subsets/goslim_generic.obo]",
      sep = "\n"
    ),
    "\n"
  )
  quit(status = 1)
}

parse_args <- function(args) {
  parsed <- list(
    cache_dir = "data/go_slim_cache",
    slim_name = "goslim_generic",
    go_basic_source = "https://current.geneontology.org/ontology/go-basic.obo",
    goslim_source = "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      usage()
    }
    parsed[[gsub("-", "_", substring(key, 3), fixed = TRUE)]] <- args[[i + 1]]
    i <- i + 2
  }
  parsed
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

required_cols <- c("ID", "Description", "setSize", "NES", "p.adjust", "core_enrichment")

source_basename <- function(source, default_name) {
  if (grepl("^https?://", source, ignore.case = TRUE)) {
    base <- basename(source)
  } else {
    base <- basename(normalizePath(source, mustWork = FALSE))
  }
  ifelse(nzchar(base), base, default_name)
}

resolve_source <- function(source, cache_dir, default_name) {
  if (file.exists(source)) {
    return(normalizePath(source))
  }

  if (!grepl("^https?://", source, ignore.case = TRUE)) {
    stop("Ontology source is not a file or URL: ", source, call. = FALSE)
  }

  dest <- file.path(cache_dir, source_basename(source, default_name))
  if (!file.exists(dest)) {
    message("Downloading ", source)
    download.file(source, destfile = dest, mode = "wb", quiet = FALSE)
  }
  normalizePath(dest)
}

parse_obo_terms <- function(path, include_part_of = TRUE) {
  lines <- readLines(path, warn = FALSE)
  term_starts <- which(lines == "[Term]")
  if (length(term_starts) == 0) {
    stop("No [Term] entries found in ", path, call. = FALSE)
  }
  term_ends <- c(term_starts[-1] - 1, length(lines))

  terms <- vector("list", length(term_starts))
  for (i in seq_along(term_starts)) {
    block <- lines[term_starts[[i]]:term_ends[[i]]]
    id_line <- grep("^id: GO:", block, value = TRUE)
    if (length(id_line) == 0) {
      next
    }

    id <- sub("^id: ", "", id_line[[1]])
    name_line <- grep("^name: ", block, value = TRUE)
    namespace_line <- grep("^namespace: ", block, value = TRUE)
    obsolete_line <- grep("^is_obsolete: true", block, value = TRUE)
    subset_lines <- grep("^subset: ", block, value = TRUE)
    is_a_lines <- grep("^is_a: GO:", block, value = TRUE)
    relationship_lines <- character()
    if (include_part_of) {
      relationship_lines <- grep("^relationship: part_of GO:", block, value = TRUE)
    }

    is_a_parents <- sub("^is_a: (GO:[0-9]+).*", "\\1", is_a_lines)
    part_of_parents <- sub("^relationship: part_of (GO:[0-9]+).*", "\\1", relationship_lines)

    terms[[i]] <- list(
      id = id,
      name = if (length(name_line) > 0) sub("^name: ", "", name_line[[1]]) else NA_character_,
      namespace = if (length(namespace_line) > 0) sub("^namespace: ", "", namespace_line[[1]]) else NA_character_,
      obsolete = length(obsolete_line) > 0,
      subsets = sub("^subset: ", "", subset_lines),
      parents = unique(c(is_a_parents, part_of_parents))
    )
  }

  terms <- terms[!vapply(terms, is.null, logical(1))]
  names(terms) <- vapply(terms, `[[`, character(1), "id")
  terms
}

nearest_slim_ancestors <- function(go_id, slim_ids, parent_map) {
  if (is.na(go_id) || !nzchar(go_id) || !(go_id %in% names(parent_map))) {
    return(data.frame(Slim_ID = NA_character_, Slim_Distance = NA_integer_))
  }

  queue <- data.frame(ID = go_id, Distance = 0L, stringsAsFactors = FALSE)
  seen <- character()

  while (nrow(queue) > 0) {
    current_distance <- min(queue$Distance)
    current <- queue[queue$Distance == current_distance, , drop = FALSE]
    queue <- queue[queue$Distance != current_distance, , drop = FALSE]
    current <- current[!(current$ID %in% seen), , drop = FALSE]
    if (nrow(current) == 0) {
      next
    }
    seen <- unique(c(seen, current$ID))

    hits <- current$ID[current$ID %in% slim_ids]
    if (length(hits) > 0) {
      return(data.frame(
        Slim_ID = unique(hits),
        Slim_Distance = current_distance,
        stringsAsFactors = FALSE
      ))
    }

    parents <- unique(unlist(parent_map[current$ID], use.names = FALSE))
    parents <- parents[nzchar(parents) & !(parents %in% seen)]
    if (length(parents) > 0) {
      queue <- rbind(
        queue,
        data.frame(ID = parents, Distance = current_distance + 1L, stringsAsFactors = FALSE)
      )
    }
  }

  data.frame(Slim_ID = NA_character_, Slim_Distance = NA_integer_)
}

count_core_genes <- function(x) {
  vapply(
    strsplit(as.character(x), "/", fixed = TRUE),
    function(parts) sum(nzchar(parts)),
    integer(1)
  )
}

core_gene_union_count <- function(x) {
  genes <- unique(unlist(strsplit(as.character(x), "/", fixed = TRUE), use.names = FALSE))
  sum(nzchar(genes))
}

safe_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("_+$", "", gsub("^_+", "", x))
}

safe_file_stub <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+$", "", gsub("^_+", "", x))
  ifelse(nzchar(x), x, "go_slim")
}

read_gsea_file <- function(path, label) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  df$Comparison <- label
  df$setSize <- as.numeric(df$setSize)
  df$NES <- as.numeric(df$NES)
  df$p.adjust <- as.numeric(df$`p.adjust`)
  df$core_count <- count_core_genes(df$core_enrichment)
  df
}

expand_to_slim <- function(df, slim_ids, slim_names, parent_map) {
  mapped <- lapply(df$ID, nearest_slim_ancestors, slim_ids = slim_ids, parent_map = parent_map)
  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    map_df <- mapped[[i]]
    repeated <- df[rep(i, nrow(map_df)), , drop = FALSE]
    repeated$Slim_ID <- map_df$Slim_ID
    repeated$Slim_Description <- unname(slim_names[map_df$Slim_ID])
    repeated$Slim_Distance <- map_df$Slim_Distance
    out[[i]] <- repeated
  }
  do.call(rbind, out)
}

summarize_slim <- function(df) {
  df <- df[!is.na(df$Slim_ID) & !is.na(df$p.adjust) & !is.na(df$NES), , drop = FALSE]
  if (nrow(df) == 0) {
    return(df)
  }

  split_df <- split(df, paste(df$Comparison, df$Slim_ID, sep = "\r"))
  rows <- lapply(split_df, function(group) {
    group <- group[order(group$p.adjust, -abs(group$NES)), , drop = FALSE]
    best <- group[1, , drop = FALSE]
    best$Slim_source_term_count <- length(unique(group$ID))
    best$Slim_core_gene_count <- core_gene_union_count(group$core_enrichment)
    best$Best_source_ID <- best$ID
    best$Best_source_Description <- best$Description
    best$neg_log10_padj <- -log10(pmax(best$p.adjust, .Machine$double.xmin))
    best
  })

  summary <- do.call(rbind, rows)
  summary <- summary[order(summary$Comparison, summary$p.adjust, -abs(summary$NES)), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("gsea_files", "labels", "outdir", "prefix")
missing <- required_args[!nzchar(vapply(required_args, function(x) args[[x]] %||% "", character(1)))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

gsea_files <- trimws(unlist(strsplit(args$gsea_files, ",", fixed = TRUE)))
labels <- trimws(unlist(strsplit(args$labels, ",", fixed = TRUE)))
slim_name <- safe_file_stub(args$slim_name)

if (length(gsea_files) != length(labels)) {
  stop("--gsea-files and --labels must have the same number of entries.", call. = FALSE)
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$cache_dir, recursive = TRUE, showWarnings = FALSE)

go_basic_path <- resolve_source(
  args$go_basic_source,
  args$cache_dir,
  "go-basic.obo"
)
goslim_path <- resolve_source(
  args$goslim_source,
  args$cache_dir,
  paste0(slim_name, ".obo")
)

terms <- parse_obo_terms(go_basic_path)
slim_terms <- parse_obo_terms(goslim_path, include_part_of = FALSE)
slim_ids <- names(slim_terms)
term_names <- vapply(terms, `[[`, character(1), "name")
slim_names <- term_names[slim_ids]
missing_slim_names <- is.na(slim_names)
if (any(missing_slim_names)) {
  slim_names[missing_slim_names] <- vapply(slim_terms[missing_slim_names], `[[`, character(1), "name")
}
parent_map <- lapply(terms, `[[`, "parents")

for (i in seq_along(gsea_files)) {
  label <- labels[[i]]
  df <- read_gsea_file(gsea_files[[i]], label)
  slim_expanded <- expand_to_slim(df, slim_ids = slim_ids, slim_names = slim_names, parent_map = parent_map)
  slim_summary <- summarize_slim(slim_expanded)

  expanded_path <- file.path(args$outdir, paste0(args$prefix, "_", slim_name, "_expanded_", safe_label(label), ".csv"))
  summary_path <- file.path(args$outdir, paste0(args$prefix, "_", slim_name, "_summary_", safe_label(label), ".csv"))
  write.csv(slim_expanded, expanded_path, row.names = FALSE)
  write.csv(slim_summary, summary_path, row.names = FALSE)
}

message("Wrote GO slim expanded and summary tables to ", normalizePath(args$outdir))
