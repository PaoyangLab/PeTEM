#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--eg"), type="character", default="gene_expression.txt", help="Gene expression/count matrix"),
  make_option(c("--et"), type="character", default="TE_expression.txt", help="TE expression/count matrix"),
  make_option(c("-o", "--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--gene-out"), type="character", default="DEG.txt", help="Gene DE output name"),
  make_option(c("--te-out"), type="character", default="DETE.txt", help="TE DE output name")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

read_mat <- function(path) {
  df <- read.table(path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  as.data.frame(lapply(df, as.numeric), row.names=rownames(df), check.names=FALSE)
}

de_col <- function(cols) grepl("^(logFC_|PValue_|FDR_|DEG_)", cols)

sample_groups <- function(cols) {
  grp <- sub("([._-]rep)?[._-]?[0-9]+$", "", cols, ignore.case=TRUE)
  if (any(grp == "")) cols else grp
}

mean_by_group <- function(df, grp) {
  groups <- unique(grp)
  out <- sapply(groups, function(g) rowMeans(df[, grp == g, drop=FALSE], na.rm=TRUE))
  out <- as.data.frame(out, check.names=FALSE)
  rownames(out) <- rownames(df)
  out
}

fallback_pair <- function(avg, si, sj) {
  logfc <- log2((avg[[sj]] + 1e-6) / (avg[[si]] + 1e-6))
  data.frame(logFC=logfc, PValue=1, FDR=1, check.names=FALSE)
}

edgeR_pair <- function(df, grp, si, sj) {
  y <- edgeR::DGEList(counts=round(as.matrix(df)), group=factor(grp))
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y)
  tab <- edgeR::topTags(edgeR::exactTest(y, pair=c(si, sj)), n=Inf)$table
  tab[rownames(df), c("logFC", "PValue", "FDR"), drop=FALSE]
}

make_de <- function(infile, outfile) {
  raw <- read_mat(infile)
  if (any(grepl("^logFC_", colnames(raw)))) {
    write.table(raw, file=file.path(opt$outdir, outfile), sep="\t", quote=FALSE, col.names=NA)
    return(invisible(NULL))
  }
  df <- raw[, !de_col(colnames(raw)), drop=FALSE]
  grp <- sample_groups(colnames(df))
  avg <- mean_by_group(df, grp)
  groups <- colnames(avg)
  if (length(groups) < 2) stop("Need at least two stages/conditions: ", infile)

  use_edger <- requireNamespace("edgeR", quietly=TRUE) &&
    all(df >= 0, na.rm=TRUE) &&
    any(table(grp) > 1)

  out <- avg
  for (cmb in combn(groups, 2, simplify=FALSE)) {
    si <- cmb[1]; sj <- cmb[2]
    stat <- if (use_edger) edgeR_pair(df, grp, si, sj) else fallback_pair(avg, si, sj)
    out[[paste0("logFC_", sj, "_", si)]] <- stat$logFC
    out[[paste0("PValue_", sj, "_", si)]] <- stat$PValue
    out[[paste0("FDR_", sj, "_", si)]] <- stat$FDR
  }
  write.table(out, file=file.path(opt$outdir, outfile), sep="\t", quote=FALSE, col.names=NA)
}

make_de(opt$eg, opt$`gene-out`)
make_de(opt$et, opt$`te-out`)
