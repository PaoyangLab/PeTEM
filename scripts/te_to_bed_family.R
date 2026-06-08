#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input TE annotation file (GFF3/GTF-like tab-delimited file with >=9 columns)"),
  make_option(c("-o", "--outdir"), type = "character", default = "annotation_bed", help = "Output directory"),
  make_option(c("--feature-types"), type = "character", default = "transposable_element,transposable_element_gene,repeat_region", help = "Comma-separated TE feature types to keep"),
  make_option(c("--name-keys"), type = "character", default = "ID,Name,gene_id,transcript_id,locus_tag", help = "Comma-separated attribute keys used to extract the TE name"),
  make_option(c("--family-keys"), type = "character", default = "family,Family,classification,Classification,class,Class,superfamily,Superfamily,te_family,repeat_family", help = "Comma-separated attribute keys used to extract TE family")
)

opt <- parse_args(OptionParser(option_list = option_list))

die <- function(step, message) {
  stop(sprintf("[ERROR] Module Annotation | te_to_bed_family.R | %s: %s", step, message), call. = FALSE)
}

if (is.null(opt$input) || opt$input == "") {
  cat("Usage: Rscript te_to_bed_family.R --input te_annotation.gff3 [--outdir annotation_bed]\n", file = stderr())
  die("argument parsing", "missing required argument: --input")
}

if (!file.exists(opt$input)) {
  die("input validation", sprintf("input annotation file not found: %s", opt$input))
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

first_attr_match <- function(attr_text, keys) {
  for (key in keys) {
    value <- parse_attr_value(attr_text, key)
    if (!is.null(value) && nzchar(value)) {
      return(value)
    }
  }
  NA_character_
}

step_name <- "read annotation"
te_annot <- tryCatch(
  read.delim(opt$input, header = FALSE, sep = "\t", quote = "", comment.char = "#", stringsAsFactors = FALSE),
  error = function(e) die(step_name, sprintf("failed to read %s: %s", opt$input, e$message))
)

if (ncol(te_annot) < 9) {
  die(step_name, sprintf("expected at least 9 columns in %s, found %d", opt$input, ncol(te_annot)))
}

colnames(te_annot)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
te_annot$type_lower <- tolower(te_annot$type)

feature_types <- tolower(trimws(strsplit(opt$`feature-types`, ",", fixed = TRUE)[[1]]))
name_keys <- trimws(strsplit(opt$`name-keys`, ",", fixed = TRUE)[[1]])
family_keys <- trimws(strsplit(opt$`family-keys`, ",", fixed = TRUE)[[1]])

step_name <- "filter TE features"
te_rows <- te_annot[te_annot$type_lower %in% feature_types, , drop = FALSE]
if (nrow(te_rows) == 0) {
  die(step_name, sprintf("no TE rows matched feature types: %s", paste(feature_types, collapse = ", ")))
}

step_name <- "parse TE attributes"
te_rows$te_name <- vapply(te_rows$attributes, first_attr_match, character(1), keys = name_keys)
te_rows$te_family <- vapply(te_rows$attributes, first_attr_match, character(1), keys = family_keys)

if (any(is.na(te_rows$te_name) | te_rows$te_name == "")) {
  die(step_name, sprintf("failed to parse TE name from column 9 for %d row(s)", sum(is.na(te_rows$te_name) | te_rows$te_name == "")))
}

if (any(is.na(te_rows$te_family) | te_rows$te_family == "")) {
  die(step_name, sprintf("failed to parse TE family from column 9 for %d row(s)", sum(is.na(te_rows$te_family) | te_rows$te_family == "")))
}

step_name <- "write outputs"
te_bed <- data.frame(
  chr = te_rows$seqid,
  start = pmax(0L, as.integer(te_rows$start) - 1L),
  end = as.integer(te_rows$end),
  name = te_rows$te_name,
  score = 0,
  strand = te_rows$strand,
  stringsAsFactors = FALSE
)
te_bed <- te_bed[order(te_bed$chr, te_bed$start, te_bed$end, te_bed$name), , drop = FALSE]

te_family <- unique(data.frame(
  te_name = te_rows$te_name,
  te_family = te_rows$te_family,
  stringsAsFactors = FALSE
))
te_family <- te_family[order(te_family$te_name, te_family$te_family), , drop = FALSE]

te_bed_path <- file.path(opt$outdir, "TE.bed")
te_family_path <- file.path(opt$outdir, "TE_family.txt")

write.table(te_bed, file = te_bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(te_family, file = te_family_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message(sprintf("[INFO] Wrote %s", te_bed_path))
message(sprintf("[INFO] Wrote %s", te_family_path))
message(sprintf("[DONE] Module Annotation | te_to_bed_family.R | completed: outputs written to %s", opt$outdir))
