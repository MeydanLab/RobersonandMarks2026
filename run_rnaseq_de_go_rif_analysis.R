#!/usr/bin/env Rscript

# Run an RNA-seq DESeq2, GO GSEA, and pathway gene-set workflow on a
# tab-delimited count matrix with user-supplied sample metadata.
#
# Example:
#   Rscript run_rnaseq_de_go_rif_analysis.R \
#     --counts path/to/counts.tsv \
#     --metadata path/to/metadata.tsv \
#     --outdir path/to/output_dir \
#     --gene-id-type SYMBOL \
#     --rif-files path/to/gene_set_a.csv,path/to/gene_set_b.csv \
#     --rif-names gene_set_a,gene_set_b
#
# Optional sample subset:
#   --samples "sample_1,sample_2,sample_3,sample_4"

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(DESeq2)
  library(clusterProfiler)
  library(AnnotationDbi)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript run_rnaseq_de_go_rif_analysis.R --counts <counts.tsv> --metadata <metadata.tsv> --outdir <output_dir>",
      "[--gene-id-col GeneID] [--length-col Length] [--gene-id-type SYMBOL|ENSEMBL]",
      "[--filter-reads 20] [--filter-min-samples 3] [--reference-inducer None]",
      "[--go-ontology MF] [--samples sample1,sample2,...]",
      "[--rif-files gene_set1.csv,gene_set2.csv] [--rif-names name1,name2]",
      "[--rif-up-padj 0.05] [--rif-up-log2fc 1]",
      "[--de-up-padj 0.05] [--de-up-log2fc 1]",
      sep = "\n"
    ),
    "\n"
  )
  quit(status = 1)
}

parse_args <- function(args) {
  parsed <- list(
    gene_id_col = "GeneID",
    length_col = "Length",
    gene_id_type = "SYMBOL",
    filter_reads = "20",
    filter_min_samples = "3",
    reference_inducer = "None",
    go_ontology = "MF",
    samples = "",
    rif_files = "",
    rif_names = "",
    rif_up_padj = "0.05",
    rif_up_log2fc = "1",
    de_up_padj = "0.05",
    de_up_log2fc = "1"
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

stop_if_missing <- function(opts, fields) {
  missing <- fields[!nzchar(vapply(fields, function(x) opts[[x]] %||% "", character(1)))]
  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

dir.create_safe <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

remove_all_zero_rows <- function(df, min_total_count = 0) {
  df[rowSums(df) > min_total_count, , drop = FALSE]
}

read_metadata_file <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

normalize_gene_id_type <- function(x) {
  toupper(trimws(x))
}

add_gene_name_column <- function(df, gene_ids, gene_id_type) {
  gene_id_type <- normalize_gene_id_type(gene_id_type)
  if (gene_id_type == "ENSEMBL") {
    symbols <- mapIds(
      org.Hs.eg.db,
      keys = gene_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    df$Gene_name <- unname(symbols[gene_ids])
  } else if (gene_id_type == "SYMBOL") {
    df$Gene_name <- gene_ids
  } else {
    df$Gene_name <- gene_ids
  }
  df
}

parse_csv_arg <- function(x) {
  values <- trimws(unlist(strsplit(x %||% "", ",", fixed = TRUE)))
  values[nzchar(values)]
}

safe_file_stub <- function(x) {
  stub <- gsub("[^A-Za-z0-9._-]+", "_", x)
  stub <- gsub("_+$", "", gsub("^_+", "", stub))
  ifelse(nzchar(stub), stub, "gene_set")
}

format_cutoff_label <- function(x) {
  label <- format(x, scientific = FALSE, trim = TRUE)
  gsub("[^A-Za-z0-9]+", "p", label)
}

read_gene_symbol_file <- function(path) {
  if (!file.exists(path)) {
    stop("RIF gene-set file does not exist: ", path, call. = FALSE)
  }
  symbols <- unique(trimws(as.character(read.csv(path, header = FALSE, stringsAsFactors = FALSE)[[1]])))
  symbols[nzchar(symbols)]
}

load_rif_modules <- function(rif_files_arg, rif_names_arg) {
  rif_files <- parse_csv_arg(rif_files_arg)
  rif_names <- parse_csv_arg(rif_names_arg)

  if (length(rif_files) == 0) {
    return(list())
  }
  if (length(rif_names) == 0) {
    rif_names <- tools::file_path_sans_ext(basename(rif_files))
  }
  if (length(rif_names) != length(rif_files)) {
    stop("--rif-files and --rif-names must have the same number of entries.", call. = FALSE)
  }

  modules <- vector("list", length(rif_files))
  for (i in seq_along(rif_files)) {
    module_name <- safe_file_stub(rif_names[[i]])
    modules[[i]] <- list(
      name = module_name,
      path = rif_files[[i]],
      symbols = read_gene_symbol_file(rif_files[[i]])
    )
  }
  names(modules) <- vapply(modules, `[[`, character(1), "name")
  modules
}

write_result_table <- function(res_obj, out_path, gene_info, gene_id_type) {
  res_df <- as.data.frame(res_obj)
  res_df$GeneID <- rownames(res_df)
  if (!is.null(gene_info) && "Length" %in% colnames(gene_info)) {
    res_df$Length <- gene_info[res_df$GeneID, "Length"]
  }
  res_df <- add_gene_name_column(res_df, res_df$GeneID, gene_id_type)
  res_df <- res_df[, c("GeneID", "Gene_name", intersect("Length", colnames(res_df)),
                       setdiff(colnames(res_df), c("GeneID", "Gene_name", "Length")))]
  res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
  write.table(res_df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  res_df
}

run_symbol_set_module <- function(
  res_df,
  contrast_label,
  gene_symbols,
  outdir,
  module_stub,
  status_label,
  padj_cutoff,
  log2fc_cutoff
) {
  cutoff_suffix <- paste0(
    "padj",
    format_cutoff_label(padj_cutoff),
    "_log2FC",
    format_cutoff_label(log2fc_cutoff)
  )
  sig_any_suffix <- paste0("padj", format_cutoff_label(padj_cutoff), "_anyLFC")

  gene_set_path <- file.path(
    outdir,
    "Tables",
    paste0(module_stub, "_symbol_set_", contrast_label, ".tsv")
  )
  filtered_path <- file.path(
    outdir,
    "Tables",
    paste0("filtered_", module_stub, "_", contrast_label, ".tsv")
  )
  filtered_up_path <- file.path(
    outdir,
    "Tables",
    paste0("filtered_", module_stub, "_up_", cutoff_suffix, "_", contrast_label, ".tsv")
  )
  filtered_sig_any_path <- file.path(
    outdir,
    "Tables",
    paste0("filtered_", module_stub, "_", sig_any_suffix, "_", contrast_label, ".tsv")
  )
  gene_set_df <- data.frame(
    Gene_name = unique(gene_symbols),
    stringsAsFactors = FALSE
  )
  write.table(
    gene_set_df,
    gene_set_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (!"Gene_name" %in% colnames(res_df)) {
    return(paste(status_label, "module skipped because Gene_name column is missing."))
  }

  res_df$Gene_name_upper <- toupper(trimws(as.character(res_df$Gene_name)))
  gene_set_upper <- toupper(trimws(unique(gene_symbols)))
  filtered_df <- res_df[res_df$Gene_name_upper %in% gene_set_upper, , drop = FALSE]
  filtered_df$Gene_name_upper <- NULL

  write.table(
    filtered_df,
    filtered_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  filtered_up_df <- subset(
    filtered_df,
    !is.na(padj) & padj <= padj_cutoff & !is.na(log2FoldChange) & log2FoldChange >= log2fc_cutoff
  )
  filtered_sig_any_df <- subset(
    filtered_df,
    !is.na(padj) & padj <= padj_cutoff & !is.na(log2FoldChange)
  )
  write.table(
    filtered_up_df,
    filtered_up_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    filtered_sig_any_df,
    filtered_sig_any_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  paste(
    status_label,
    "module wrote",
    nrow(filtered_df),
    "pathway rows and",
    nrow(filtered_up_df),
    "significant up rows at padj <=",
    padj_cutoff,
    "and log2FoldChange >=",
    log2fc_cutoff,
    ";",
    nrow(filtered_sig_any_df),
    "significant rows at padj <=",
    padj_cutoff
  )
}

run_gsea <- function(res_obj, gene_id_type, ontology, out_prefix) {
  ranked <- res_obj$log2FoldChange
  names(ranked) <- rownames(res_obj)
  ranked <- ranked[!is.na(ranked)]
  ranked <- sort(ranked, decreasing = TRUE)

  if (length(ranked) < 10) {
    return(list(success = FALSE, message = "Too few ranked genes for GSEA"))
  }

  gse_args <- list(
    geneList = ranked,
    ont = ontology,
    keyType = gene_id_type,
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = org.Hs.eg.db,
    pAdjustMethod = "BH"
  )

  # Use 10,000 permutations when supported by the installed
  # clusterProfiler/gseGO version.
  if ("nPerm" %in% names(formals(clusterProfiler::gseGO))) {
    gse_args$nPerm <- 10000
  }

  gse <- tryCatch(
    do.call(clusterProfiler::gseGO, gse_args),
    error = function(e) e
  )

  if (inherits(gse, "error")) {
    return(list(success = FALSE, message = conditionMessage(gse)))
  }

  gse_df <- as.data.frame(gse)
  if (nrow(gse_df) == 0) {
    return(list(success = FALSE, message = "No enriched terms returned"))
  }

  write.csv(gse_df, paste0(out_prefix, ".csv"), row.names = FALSE)
  list(success = TRUE, message = paste("Wrote", nrow(gse_df), "GSEA rows"))
}

validate_design <- function(metadata, design_formula, selected_samples) {
  mm <- model.matrix(design_formula, metadata)
  if (nrow(metadata) <= qr(mm)$rank) {
    stop(
      paste(
        "The selected samples do not provide enough replication to fit the DESeq2 model.",
        "This usually means too few samples per group or only one sample per condition.",
        "Selected samples:",
        paste(selected_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

derive_batch <- function(metadata) {
  if ("Batch" %in% colnames(metadata)) {
    return(metadata$Batch)
  }

  if ("Sequencing_pool" %in% colnames(metadata)) {
    pooled <- ifelse(
      metadata$Sequencing_pool == "GL069",
      "G69",
      ifelse(
        metadata$Sequencing_pool %in% c("GL072", "GL083"),
        "G72_83",
        ifelse(metadata$Sequencing_pool == "GL136", "G136", metadata$Sequencing_pool)
      )
    )
    return(pooled)
  }

  rep("single_batch", nrow(metadata))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
stop_if_missing(args, c("counts", "metadata", "outdir"))

counts_path <- args$counts
metadata_path <- args$metadata
outdir <- args$outdir
gene_id_col <- args$gene_id_col
length_col <- args$length_col
gene_id_type <- normalize_gene_id_type(args$gene_id_type)
filter_reads <- as.integer(args$filter_reads)
filter_min_samples <- as.integer(args$filter_min_samples)
reference_inducer <- args$reference_inducer
go_ontology <- toupper(args$go_ontology)
rif_modules <- load_rif_modules(args$rif_files, args$rif_names)
rif_up_padj <- as.numeric(args$rif_up_padj)
rif_up_log2fc <- as.numeric(args$rif_up_log2fc)
de_up_padj <- as.numeric(args$de_up_padj)
de_up_log2fc <- as.numeric(args$de_up_log2fc)
selected_samples <- trimws(unlist(strsplit(args$samples, ",", fixed = TRUE)))
selected_samples <- selected_samples[nzchar(selected_samples)]

if (!(gene_id_type %in% c("SYMBOL", "ENSEMBL"))) {
  stop("gene-id-type must be SYMBOL or ENSEMBL.", call. = FALSE)
}
if (is.na(rif_up_padj) || rif_up_padj < 0 || rif_up_padj > 1) {
  stop("--rif-up-padj must be a number between 0 and 1.", call. = FALSE)
}
if (is.na(rif_up_log2fc)) {
  stop("--rif-up-log2fc must be numeric.", call. = FALSE)
}
if (is.na(de_up_padj) || de_up_padj < 0 || de_up_padj > 1) {
  stop("--de-up-padj must be a number between 0 and 1.", call. = FALSE)
}
if (is.na(de_up_log2fc)) {
  stop("--de-up-log2fc must be numeric.", call. = FALSE)
}
de_up_suffix <- paste0(
  "padj",
  format_cutoff_label(de_up_padj),
  "_log2FC",
  format_cutoff_label(de_up_log2fc)
)

dir.create_safe(outdir)
dir.create_safe(file.path(outdir, "DE"))
dir.create_safe(file.path(outdir, "Tables"))
dir.create_safe(file.path(outdir, "Mean_norm_counts"))

counts_df <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read_metadata_file(metadata_path)

if (!(gene_id_col %in% colnames(counts_df))) {
  stop("The count file is missing gene ID column: ", gene_id_col, call. = FALSE)
}

if (!("sample" %in% colnames(metadata))) {
  stop("Metadata must contain a 'sample' column matching the count-matrix headers.", call. = FALSE)
}

required_meta_cols <- c("sample", "Genotype", "Group", "Inducer")
missing_meta <- setdiff(required_meta_cols, colnames(metadata))
if (length(missing_meta) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

for (col_name in required_meta_cols) {
  bad <- is.na(metadata[[col_name]]) | trimws(metadata[[col_name]]) == ""
  if (any(bad)) {
    stop("Metadata column '", col_name, "' contains blank values. Fill the template before running.", call. = FALSE)
  }
}

sample_columns <- setdiff(colnames(counts_df), c(gene_id_col, length_col))
if (length(selected_samples) > 0) {
  missing_selected <- setdiff(selected_samples, sample_columns)
  if (length(missing_selected) > 0) {
    stop("These selected samples are not in the count file: ",
         paste(missing_selected, collapse = ", "),
         call. = FALSE)
  }
  sample_columns <- selected_samples
}

metadata <- metadata[metadata$sample %in% sample_columns, , drop = FALSE]
metadata <- metadata[match(sample_columns, metadata$sample), , drop = FALSE]

if (any(is.na(metadata$sample))) {
  stop("Metadata does not cover all selected count columns.", call. = FALSE)
}

counts_mat <- as.matrix(counts_df[, sample_columns, drop = FALSE])
mode(counts_mat) <- "numeric"
rownames(counts_mat) <- counts_df[[gene_id_col]]

if (anyDuplicated(rownames(counts_mat))) {
  stop("The selected count matrix contains duplicated gene IDs.", call. = FALSE)
}

counts_mat <- remove_all_zero_rows(counts_mat, min_total_count = 0)
metadata$Read_counts <- colSums(counts_mat)
metadata$Sample_id <- metadata$sample
metadata$Group_gt_ind <- paste0(metadata$Genotype, metadata$Inducer)

metadata$Batch <- derive_batch(metadata)

factor_cols <- intersect(c("Genotype", "Group", "Inducer", "Group_gt_ind", "Batch"), colnames(metadata))
metadata[factor_cols] <- lapply(metadata[factor_cols], factor)
rownames(metadata) <- metadata$sample

design_formula <- if (nlevels(metadata$Batch) > 1) {
  ~ Batch + Group_gt_ind
} else {
  ~ Group_gt_ind
}

validate_design(metadata, design_formula, sample_columns)

gene_info <- data.frame(GeneID = rownames(counts_mat), stringsAsFactors = FALSE)
if (length_col %in% colnames(counts_df)) {
  length_map <- counts_df[, c(gene_id_col, length_col), drop = FALSE]
  rownames(length_map) <- length_map[[gene_id_col]]
  gene_info$Length <- length_map[gene_info$GeneID, length_col]
}
rownames(gene_info) <- gene_info$GeneID

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_mat),
  colData = metadata,
  design = design_formula
)

keep <- rowSums(counts(dds) >= filter_reads) >= filter_min_samples
dds <- dds[keep, ]
gene_info <- gene_info[rownames(dds), , drop = FALSE]

if (nrow(dds) == 0) {
  stop("No genes passed the count filter.", call. = FALSE)
}

dds <- estimateSizeFactors(dds)

write.table(
  counts(dds, normalized = TRUE),
  file.path(outdir, "Mean_norm_counts", "normalized_counts_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  metadata,
  file.path(outdir, "sample_metadata_used.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

summary_lines <- c(
  paste("Counts:", counts_path),
  paste("Metadata:", metadata_path),
  paste("Selected samples:", paste(sample_columns, collapse = ", ")),
  paste("Gene ID type:", gene_id_type),
  paste("DESeq2 design:", deparse(design_formula)),
  paste("Genes after count filter:", nrow(dds)),
  paste("DE significant-up cutoff: padj <=", de_up_padj, "and log2FoldChange >", de_up_log2fc),
  paste("RIF modules:", if (length(rif_modules) > 0) paste(names(rif_modules), collapse = ", ") else "none"),
  paste("RIF significant-up cutoff: padj <=", rif_up_padj, "and log2FoldChange >=", rif_up_log2fc)
)

gsea_notes <- character()

for (genotype in levels(metadata$Genotype)) {
  dds_sub <- dds[, dds$Genotype == genotype]
  if (ncol(dds_sub) < 2) {
    gsea_notes <- c(gsea_notes, paste("Skipped genotype", genotype, "because it has fewer than 2 samples."))
    next
  }

  dds_sub$Group_gt_ind <- droplevels(dds_sub$Group_gt_ind)
  dds_sub$Batch <- droplevels(dds_sub$Batch)

  ref_group <- paste0(genotype, reference_inducer)
  if (!(ref_group %in% levels(dds_sub$Group_gt_ind))) {
    gsea_notes <- c(
      gsea_notes,
      paste("Skipped genotype", genotype, "because reference group", ref_group, "is absent.")
    )
    next
  }

  dds_sub$Group_gt_ind <- relevel(dds_sub$Group_gt_ind, ref = ref_group)
  if (nlevels(dds_sub$Batch) > 1) {
    design(dds_sub) <- ~ Batch + Group_gt_ind
  } else {
    design(dds_sub) <- ~ Group_gt_ind
  }

  validate_design(as.data.frame(colData(dds_sub)), design(dds_sub), colnames(dds_sub))
  dds_sub <- DESeq(dds_sub)

  for (inducer in setdiff(levels(dds_sub$Inducer), reference_inducer)) {
    target_group <- paste0(genotype, inducer)
    result_name <- paste0("Group_gt_ind_", target_group, "_vs_", ref_group)
    if (!(result_name %in% resultsNames(dds_sub))) {
      gsea_notes <- c(gsea_notes, paste("Result name not found:", result_name))
      next
    }

    res <- results(dds_sub, list(c(result_name)))
    base_stub <- paste0(genotype, "_", inducer, "_vs_", reference_inducer)

    res_df <- write_result_table(
      res,
      file.path(outdir, "DE", paste0("DE_", base_stub, ".tsv")),
      gene_info,
      gene_id_type
    )
    subset_dds <- dds_sub[, dds_sub$Inducer %in% c(reference_inducer, inducer)]
    mean_norm_counts <- apply(counts(subset_dds, normalized = TRUE), 1, function(x) mean(x))
    write.table(
      x = as.data.frame(mean_norm_counts),
      file = file.path(outdir, "Mean_norm_counts", paste0("mean_norm_counts_", base_stub, ".txt")),
      sep = "\t"
    )

    gsea_status <- run_gsea(
      res,
      gene_id_type = gene_id_type,
      ontology = go_ontology,
      out_prefix = file.path(outdir, "Tables", paste0("GSEA_", base_stub, "_", go_ontology))
    )

    sig_ids <- rownames(subset(
      as.data.frame(res),
      !is.na(padj) & padj <= de_up_padj & !is.na(log2FoldChange) & log2FoldChange > de_up_log2fc
    ))
    write.table(
      data.frame(GeneID = sig_ids),
      file.path(outdir, "DE", paste0("DE_up_genes_", de_up_suffix, "_", base_stub, ".tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    rif_statuses <- character()
    for (rif_module in rif_modules) {
      rif_statuses <- c(
        rif_statuses,
        run_symbol_set_module(
          res_df = res_df,
          contrast_label = base_stub,
          gene_symbols = rif_module$symbols,
          outdir = outdir,
          module_stub = rif_module$name,
          status_label = rif_module$name,
          padj_cutoff = rif_up_padj,
          log2fc_cutoff = rif_up_log2fc
        )
      )
    }

    summary_lines <- c(
      summary_lines,
      paste("Contrast:", base_stub),
      paste("  significant up genes:", length(sig_ids)),
      paste("  GSEA status:", gsea_status$message)
    )
    if (length(rif_statuses) > 0) {
      summary_lines <- c(summary_lines, paste("  RIF status:", rif_statuses))
    }
  }
}

writeLines(c(summary_lines, gsea_notes), file.path(outdir, "run_summary.txt"))
cat(paste(c(summary_lines, gsea_notes), collapse = "\n"), "\n")
