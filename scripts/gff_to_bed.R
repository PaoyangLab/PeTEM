#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-g", "--gff"), type = "character", help = "Input GFF3/GTF file"),
  make_option(c("-o", "--outdir"), type = "character", default = "annotation_bed", help = "Output directory"),
  make_option(c("--gene-feature"), type = "character", default = "gene", help = "Feature name used for gene rows (default: gene)"),
  make_option(c("--protein-coding-only"), type = "character", default = "y", help = "Filter gene rows to protein_coding entries using column 9 (y/n, default: y)")
)

opt <- parse_args(OptionParser(option_list = option_list))

die <- function(step, message) {
  stop(sprintf("[ERROR] Module Annotation | gff_to_bed.R | %s: %s", step, message), call. = FALSE)
}

if (is.null(opt$gff) || opt$gff == "") {
  cat("Usage: Rscript gff_to_bed.R --gff annotation.gff3 [--outdir annotation_bed] [--gene-feature gene] [--protein-coding-only y]\n", file = stderr())
  die("argument parsing", "missing required argument: --gff")
}

if (!file.exists(opt$gff)) {
  die("input validation", sprintf("input annotation file not found: %s", opt$gff))
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

parse_attr_value <- function(attr_text, key) {
  patterns <- c(
    sprintf("(^|;)%s=([^;]+)", key),
    sprintf("(^|;)%s[[:space:]]+\"([^\"]+)\"", key),
    sprintf("(^|;)%s[[:space:]]+([^;]+)", key)
  )

  for (pattern in patterns) {
    m <- regexec(pattern, attr_text, perl = TRUE)
    hit <- regmatches(attr_text, m)[[1]]
    if (length(hit) >= 3) {
      return(trimws(hit[3]))
    }
  }
  NULL
}

extract_attr_vector <- function(attr_vec, keys) {
  values <- rep(NA_character_, length(attr_vec))

  assign_matches <- function(pattern, group = "\\1") {
    pending <- is.na(values) | values == ""
    if (!any(pending)) {
      return(NULL)
    }
    matches <- grepl(pattern, attr_vec[pending], perl = TRUE)
    if (!any(matches)) {
      return(NULL)
    }
    extracted <- sub(pattern, group, attr_vec[pending][matches], perl = TRUE)
    values[pending][matches] <<- trimws(extracted)
  }

  for (key in keys) {
    key_escaped <- gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", key, perl = TRUE)
    assign_matches(sprintf(".*(?:^|;)%s=([^;]+).*", key_escaped))
    assign_matches(sprintf(".*(?:^|;)%s[[:space:]]+\"([^\"]+)\".*", key_escaped))
    assign_matches(sprintf(".*(?:^|;)%s[[:space:]]+([^;]+).*", key_escaped))
  }

  values
}

normalize_attr_values <- function(values, remove_transcript_suffix = FALSE) {
  values <- trimws(values)
  values <- sub("^[^:]+:", "", values)
  values <- gsub("[[:space:]]+", "_", values)
  if (remove_transcript_suffix) {
    values <- sub("\\.[0-9]+$", "", values)
  }
  values[values == ""] <- NA_character_
  values
}

extract_gene_names <- function(attr_vec) {
  normalize_attr_values(
    extract_attr_vector(attr_vec, c("gene_id", "locus_tag", "ID", "gene", "Parent", "transcript_id", "Name")),
    remove_transcript_suffix = FALSE
  )
}

extract_feature_names <- function(attr_vec, fallback_prefix) {
  values <- normalize_attr_values(
    extract_attr_vector(attr_vec, c("gene_id", "gene", "Parent", "transcript_id", "locus_tag", "Name", "ID")),
    remove_transcript_suffix = TRUE
  )
  missing <- is.na(values) | values == ""
  if (any(missing)) {
    values[missing] <- paste0(fallback_prefix, "_", seq_len(sum(missing)))
  }
  values
}

is_protein_coding <- function(attr_text) {
  lowered <- tolower(attr_text)
  grepl("protein_coding", lowered, fixed = TRUE)
}

to_bed6 <- function(df, name_values) {
  out <- data.frame(
    chr = df[[1]],
    start = pmax(0L, as.integer(df[[4]]) - 1L),
    end = as.integer(df[[5]]),
    name = name_values,
    score = 0,
    strand = df[[7]],
    stringsAsFactors = FALSE
  )
  out <- out[complete.cases(out[, c("chr", "start", "end", "name", "strand")]), , drop = FALSE]
  out <- out[order(out$chr, out$start, out$end, out$name), , drop = FALSE]
  out
}

write_bed <- function(df, filename) {
  out_path <- file.path(opt$outdir, filename)
  write.table(df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(sprintf("[INFO] Wrote %s", out_path))
}

step_name <- "read annotation"
gff <- tryCatch(
  read.delim(opt$gff, header = FALSE, sep = "\t", quote = "", comment.char = "#", stringsAsFactors = FALSE),
  error = function(e) die(step_name, sprintf("failed to read %s: %s", opt$gff, e$message))
)

if (ncol(gff) < 9) {
  die(step_name, sprintf("expected at least 9 columns in %s, found %d", opt$gff, ncol(gff)))
}

colnames(gff)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff$type_lower <- tolower(gff$type)

needed_types <- unique(c(
  tolower(opt$`gene-feature`),
  "cds",
  "exon",
  "five_prime_utr",
  "5utr",
  "five_prime_utr_primary_transcript",
  "three_prime_utr",
  "3utr",
  "three_prime_utr_primary_transcript"
))
gff <- gff[gff$type_lower %in% needed_types, , drop = FALSE]

step_name <- "parse gene names"
gene_rows <- gff[gff$type_lower == tolower(opt$`gene-feature`), , drop = FALSE]
if (tolower(opt$`protein-coding-only`) == "y") {
  gene_rows <- gene_rows[vapply(gene_rows$attributes, is_protein_coding, logical(1)), , drop = FALSE]
}

if (nrow(gene_rows) == 0) {
  die("gene extraction", "no gene rows matched the selected filters")
}

gene_rows$gene_name <- extract_gene_names(gene_rows$attributes)
missing_gene_names <- is.na(gene_rows$gene_name) | gene_rows$gene_name == ""
if (any(missing_gene_names)) {
  die("gene extraction", sprintf("failed to parse gene name from column 9 for %d gene row(s)", sum(missing_gene_names)))
}

step_name <- "extract feature beds"

feature_map <- list(
  CDS = c("cds"),
  exon = c("exon"),
  UTR5 = c("five_prime_utr", "5utr", "five_prime_utr_primary_transcript"),
  UTR3 = c("three_prime_utr", "3utr", "three_prime_utr_primary_transcript")
)

gene_bed <- to_bed6(gene_rows, gene_rows$gene_name)
write_bed(gene_bed, "gene.bed")

for (output_name in names(feature_map)) {
  feature_types <- feature_map[[output_name]]
  rows <- gff[gff$type_lower %in% feature_types, , drop = FALSE]
  if (nrow(rows) == 0) {
    message(sprintf("[INFO] No rows found for %s; skipping %s.bed", paste(feature_types, collapse = ", "), output_name))
    next
  }

  feature_names <- extract_feature_names(rows$attributes, tolower(output_name))
  if (tolower(opt$`protein-coding-only`) == "y") {
    keep <- feature_names %in% gene_bed$name
    rows <- rows[keep, , drop = FALSE]
    feature_names <- feature_names[keep]
    if (nrow(rows) == 0) {
      message(sprintf("[INFO] No protein-coding rows retained for %s; skipping %s.bed", paste(feature_types, collapse = ", "), output_name))
      next
    }
  }
  write_bed(to_bed6(rows, feature_names), sprintf("%s.bed", output_name))
}

message(sprintf("[DONE] Module Annotation | gff_to_bed.R | completed: outputs written to %s", opt$outdir))
