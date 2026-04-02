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
script_file <- if (nzchar(script_arg)) sub("^--file=", "", script_arg) else file.path(getwd(), "run_fgsea_preranked.R")
script_dir <- normalizePath(dirname(script_file), mustWork = FALSE)
skill_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
default_stats <- file.path(skill_dir, "examples", "toy_ranked_stats.tsv")
default_pathways <- file.path(skill_dir, "examples", "toy_pathways.tsv")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list(
    stats = default_stats,
    pathways = default_pathways,
    lib_dir = Sys.getenv("BIOC_SKILL_R_LIB", unset = ""),
    install_missing = FALSE,
    min_size = 1L,
    max_size = 500L,
    seed = 7L,
    out = NULL,
    describe_toy = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--stats") {
      parsed$stats <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--pathways") {
      parsed$pathways <- args[[i + 1L]]
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
    } else if (key == "--seed") {
      parsed$seed <- as.integer(args[[i + 1L]])
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
          "Usage: Rscript run_fgsea_preranked.R [options]",
          "",
          "Options:",
          "  --stats PATH           Ranked stats TSV with columns gene,score",
          "  --pathways PATH        Pathway membership TSV with columns pathway,gene",
          "  --lib-dir PATH         Optional R library path for Bioconductor installs",
          "  --install-missing      Install fgsea into --lib-dir if missing",
          "  --min-size N           Minimum pathway size (default 1)",
          "  --max-size N           Maximum pathway size (default 500)",
          "  --seed N               RNG seed for deterministic runs (default 7)",
          "  --describe-toy         Emit toy-input summary JSON without requiring fgsea",
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

read_stats <- function(path) {
  table <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("gene", "score") %in% names(table)))
  stopifnot(!anyDuplicated(table$gene))
  stats <- table$score
  names(stats) <- table$gene
  stats[order(stats, decreasing = TRUE)]
}

read_pathways <- function(path) {
  table <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("pathway", "gene") %in% names(table)))
  split(table$gene, table$pathway)
}

describe_toy <- function(stats_path, pathways_path) {
  stats_table <- utils::read.delim(stats_path, sep = "\t", stringsAsFactors = FALSE)
  pathways_table <- utils::read.delim(pathways_path, sep = "\t", stringsAsFactors = FALSE)
  pathway_sizes <- sort(table(pathways_table$pathway), decreasing = TRUE)
  list(
    stats_file = normalizePath(stats_path),
    pathways_file = normalizePath(pathways_path),
    gene_count = nrow(stats_table),
    pathway_count = length(unique(pathways_table$pathway)),
    strongest_positive_gene = stats_table$gene[which.max(stats_table$score)],
    strongest_negative_gene = stats_table$gene[which.min(stats_table$score)],
    pathway_sizes = as.list(stats::setNames(as.integer(pathway_sizes), names(pathway_sizes)))
  )
}

run_fgsea <- function(stats_path, pathways_path, min_size, max_size, seed) {
  stats <- read_stats(stats_path)
  pathways <- read_pathways(pathways_path)
  set.seed(seed)
  result <- fgsea::fgsea(
    pathways = pathways,
    stats = stats,
    minSize = min_size,
    maxSize = max_size,
    eps = 0,
    BPPARAM = BiocParallel::SerialParam()
  )
  result <- result[order(result$padj, result$pval, decreasing = FALSE), ]
  rows <- vector("list", nrow(result))
  for (idx in seq_len(nrow(result))) {
    row <- result[idx, ]
    rows[[idx]] <- list(
      pathway = as.character(row$pathway),
      size = as.integer(row$size),
      ES = unname(as.numeric(row$ES)),
      NES = unname(as.numeric(row$NES)),
      pval = unname(as.numeric(row$pval)),
      padj = unname(as.numeric(row$padj)),
      leadingEdge = as.character(unlist(row$leadingEdge))
    )
  }
  list(
    stats_file = normalizePath(stats_path),
    pathways_file = normalizePath(pathways_path),
    package = "fgsea",
    package_version = as.character(utils::packageVersion("fgsea")),
    seed = as.integer(seed),
    min_size = as.integer(min_size),
    max_size = as.integer(max_size),
    gene_count = length(stats),
    pathway_count = length(pathways),
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
    payload <- describe_toy(args$stats, args$pathways)
    write_json(payload, args$out)
    return(invisible(NULL))
  }
  ensure_bioc_package("fgsea", lib_dir = args$lib_dir, install_missing = args$install_missing)
  payload <- run_fgsea(args$stats, args$pathways, args$min_size, args$max_size, args$seed)
  write_json(payload, args$out)
}

main()
