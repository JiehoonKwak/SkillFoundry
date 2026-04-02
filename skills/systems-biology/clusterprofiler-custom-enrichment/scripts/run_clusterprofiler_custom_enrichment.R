#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

suppressWarnings(suppressMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("This script requires the 'jsonlite' R package.")
  }
}))

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1] %||% ""
script_file <- if (nzchar(script_arg)) sub("^--file=", "", script_arg) else file.path(getwd(), "run_clusterprofiler_custom_enrichment.R")
script_dir <- normalizePath(dirname(script_file), mustWork = FALSE)
skill_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
default_gene_file <- file.path(skill_dir, "examples", "toy_query_genes.txt")
default_term2gene <- file.path(skill_dir, "examples", "toy_term2gene.tsv")
default_term2name <- file.path(skill_dir, "examples", "toy_term2name.tsv")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list(
    genes = default_gene_file,
    term2gene = default_term2gene,
    term2name = default_term2name,
    lib_dir = Sys.getenv("BIOC_SKILL_R_LIB", unset = ""),
    install_missing = FALSE,
    min_size = 1L,
    max_size = 500L,
    out = NULL,
    describe_toy = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--genes") {
      parsed$genes <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--term2gene") {
      parsed$term2gene <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--term2name") {
      parsed$term2name <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--lib-dir") {
      parsed$lib_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--install-missing") {
      parsed$install_missing <- TRUE
      i <- i + 1L
    } else if (key == "--min-size") {
      parsed$min_size <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--max-size") {
      parsed$max_size <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--out") {
      parsed$out <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--describe-toy") {
      parsed$describe_toy <- TRUE
      i <- i + 1L
    } else if (key %in% c("-h", "--help")) {
      cat(
        paste(
          "Usage: Rscript run_clusterprofiler_custom_enrichment.R [options]",
          "",
          "Options:",
          "  --genes PATH           Query gene file with one gene per line",
          "  --term2gene PATH       TERM2GENE TSV with columns term,gene",
          "  --term2name PATH       TERM2NAME TSV with columns term,name",
          "  --lib-dir PATH         Optional R library path for Bioconductor installs",
          "  --install-missing      Install clusterProfiler into --lib-dir if missing",
          "  --min-size N           Minimum gene-set size (default 1)",
          "  --max-size N           Maximum gene-set size (default 500)",
          "  --describe-toy         Emit toy-input summary JSON without requiring clusterProfiler",
          "  --out PATH             Optional JSON output path",
          sep = "\n"
        )
      )
      quit(status = 0L)
    } else {
      stop(sprintf("Unknown argument: %s", key))
    }
  }
  parsed
}

configure_lib <- function(lib_dir) {
  if (!nzchar(lib_dir)) {
    return(invisible(NULL))
  }
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  .libPaths(unique(c(normalizePath(lib_dir), .libPaths())))
}

ensure_bioc_package <- function(pkg, lib_dir = "", install_missing = FALSE) {
  configure_lib(lib_dir)
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  if (!install_missing) {
    stop(sprintf(
      "Missing required package '%s'. Re-run with --install-missing --lib-dir <path> or set BIOC_SKILL_R_LIB.",
      pkg
    ))
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    stop("Missing required package 'BiocManager' needed to install Bioconductor dependencies.")
  }
  target_lib <- if (nzchar(lib_dir)) lib_dir else .libPaths()[1]
  dir.create(target_lib, recursive = TRUE, showWarnings = FALSE)
  BiocManager::install(pkg, lib = target_lib, ask = FALSE, update = FALSE, quiet = TRUE)
  configure_lib(target_lib)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Failed to install required package '%s'.", pkg))
  }
  invisible(TRUE)
}

read_genes <- function(path) {
  genes <- trimws(readLines(path, warn = FALSE))
  unique(genes[nzchar(genes)])
}

read_term_table <- function(path, expected) {
  table <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(expected %in% names(table)))
  table
}

describe_toy <- function(genes_path, term2gene_path, term2name_path) {
  genes <- read_genes(genes_path)
  term2gene <- read_term_table(term2gene_path, c("term", "gene"))
  term2name <- read_term_table(term2name_path, c("term", "name"))
  gene_sets <- sort(table(term2gene$term), decreasing = TRUE)
  list(
    genes_file = normalizePath(genes_path),
    term2gene_file = normalizePath(term2gene_path),
    term2name_file = normalizePath(term2name_path),
    query_gene_count = length(genes),
    term_count = length(unique(term2gene$term)),
    term_names = as.list(stats::setNames(term2name$name, term2name$term)),
    term_sizes = as.list(stats::setNames(as.integer(gene_sets), names(gene_sets)))
  )
}

run_clusterprofiler <- function(genes_path, term2gene_path, term2name_path, min_size, max_size) {
  genes <- read_genes(genes_path)
  term2gene <- read_term_table(term2gene_path, c("term", "gene"))
  term2name <- read_term_table(term2name_path, c("term", "name"))
  result <- clusterProfiler::enricher(
    gene = genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    minGSSize = min_size,
    maxGSSize = max_size,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  result_df <- if (is.null(result)) data.frame() else as.data.frame(result)
  rows <- vector("list", nrow(result_df))
  if (nrow(result_df) > 0) {
    result_df <- result_df[order(result_df$p.adjust, result_df$pvalue, decreasing = FALSE), ]
    for (idx in seq_len(nrow(result_df))) {
      row <- result_df[idx, ]
      rows[[idx]] <- list(
        term = as.character(row$ID),
        description = as.character(row$Description),
        gene_ratio = as.character(row$GeneRatio),
        bg_ratio = as.character(row$BgRatio),
        pvalue = unname(as.numeric(row$pvalue)),
        p_adjust = unname(as.numeric(row$p.adjust)),
        qvalue = if ("qvalue" %in% names(row)) unname(as.numeric(row$qvalue)) else NULL,
        gene_ids = strsplit(as.character(row$geneID), "/", fixed = TRUE)[[1]],
        count = as.integer(row$Count)
      )
    }
  }
  list(
    genes_file = normalizePath(genes_path),
    term2gene_file = normalizePath(term2gene_path),
    term2name_file = normalizePath(term2name_path),
    package = "clusterProfiler",
    package_version = as.character(utils::packageVersion("clusterProfiler")),
    query_gene_count = length(genes),
    term_count = length(unique(term2gene$term)),
    result_count = length(rows),
    results = rows
  )
}

write_json <- function(payload, out_path = NULL) {
  text <- jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (is.null(out_path)) {
    cat(text, "\n")
    return(invisible(NULL))
  }
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, out_path, useBytes = TRUE)
}

main <- function() {
  args <- parse_args()
  if (args$describe_toy) {
    payload <- describe_toy(args$genes, args$term2gene, args$term2name)
    write_json(payload, args$out)
    return(invisible(NULL))
  }
  ensure_bioc_package("clusterProfiler", lib_dir = args$lib_dir, install_missing = args$install_missing)
  payload <- run_clusterprofiler(args$genes, args$term2gene, args$term2name, args$min_size, args$max_size)
  write_json(payload, args$out)
}

main()
