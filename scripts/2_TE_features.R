nearest_tss_distance <- function(te, gene) {
  gene_start <- as.numeric(gene[[2]])
  gene_end <- as.numeric(gene[[3]])
  gene_strand <- as.character(gene[[6]])
  if (any(!is.finite(gene_start) | !is.finite(gene_end) | gene_end <= gene_start) ||
      any(!(gene_strand %in% c("+", "-")))) {
    stop("Gene BED needs valid intervals and +/- strands to locate TSSs")
  }

  tss <- ifelse(gene_strand == "+", gene_start, gene_end - 1)
  tss_by_chr <- split(tss, as.character(gene[[1]]))
  distance <- rep(NA_real_, nrow(te))

  for (chrom in unique(te$chrom)) {
    rows <- which(te$chrom == chrom)
    positions <- sort(unique(tss_by_chr[[chrom]]))
    if (length(positions) == 0) {
      warning("No TSS on ", chrom, "; TE-TSS distance set to NA")
      next
    }

    first <- te$start[rows]
    last <- te$end[rows] - 1
    before_first <- findInterval(first - 1, positions)
    through_last <- findInterval(last, positions)
    within_te <- through_last > before_first

    left <- rep(Inf, length(rows))
    has_left <- before_first > 0
    left[has_left] <- first[has_left] - positions[before_first[has_left]]

    right <- rep(Inf, length(rows))
    next_idx <- through_last + 1
    has_right <- next_idx <= length(positions)
    right[has_right] <- positions[next_idx[has_right]] - last[has_right]

    distance[rows] <- ifelse(within_te, 0, pmin(left, right))
  }

  distance
}

run_te_features <- function(all_te, embedded_te, gene_bed, fasta, outdir,
                            group_ids=list(), output_prefix="OUTPUT_2", group_colors=NULL) {
  if (!nzchar(Sys.which("bedtools"))) stop("bedtools is required for --feature y")

  te <- data.frame(
    chrom=as.character(all_te[[1]]),
    start=as.numeric(all_te[[2]]),
    end=as.numeric(all_te[[3]]),
    id=as.character(all_te[[4]]),
    score=as.character(all_te[[5]]),
    strand=as.character(all_te[[6]]),
    stringsAsFactors=FALSE
  )
  te <- te[!duplicated(te[c("chrom", "start", "end", "id")]), , drop=FALSE]
  if (nrow(te) == 0 || any(!is.finite(te$start) | !is.finite(te$end) | te$start < 0 | te$end <= te$start)) {
    stop("All-TE BED contains invalid intervals")
  }

  te_key <- paste(te$chrom, te$start, te$end, te$id, sep="\t")
  embedded_key <- unique(paste(embedded_te[[1]], embedded_te[[2]], embedded_te[[3]], embedded_te[[4]], sep="\t"))
  if (any(!(embedded_key %in% te_key))) {
    stop("Embedded TE intervals do not match the all-TE BED; check chromosome names and coordinates")
  }
  te$embedded <- te_key %in% embedded_key

  fai <- read.delim(paste0(fasta, ".fai"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  if (ncol(fai) < 2 || anyDuplicated(fai[[1]])) stop("Invalid genome FASTA index")
  chrom_lengths <- setNames(as.numeric(fai[[2]]), as.character(fai[[1]]))
  missing_chrom <- setdiff(unique(te$chrom), names(chrom_lengths))
  if (length(missing_chrom) > 0) {
    stop("TE chromosomes absent from FASTA: ", paste(missing_chrom, collapse=", "))
  }
  if (any(te$end > chrom_lengths[te$chrom])) stop("TE intervals extend beyond FASTA chromosome lengths")

  gene <- read.delim(gene_bed, header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE)
  if (ncol(gene) < 6 || nrow(gene) == 0) stop("Gene BED must have at least six columns and one gene")
  distance <- nearest_tss_distance(te, gene)

  te_file <- tempfile("module2_te_", fileext=".bed")
  nuc_file <- tempfile("module2_nuc_", fileext=".txt")
  err_file <- tempfile("module2_nuc_", fileext=".err")
  on.exit(unlink(c(te_file, nuc_file, err_file)), add=TRUE)
  te_bed <- te[c("chrom", "start", "end", "id", "score", "strand")]
  te_bed$start <- sprintf("%.0f", te_bed$start)
  te_bed$end <- sprintf("%.0f", te_bed$end)
  write.table(te_bed,
              te_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  status <- suppressWarnings(system2("bedtools",
    args=c("nuc", "-fi", shQuote(fasta), "-bed", shQuote(te_file)),
    stdout=nuc_file, stderr=err_file))
  if (status != 0) stop("bedtools nuc failed: ", paste(readLines(err_file, warn=FALSE), collapse="; "))

  nuc <- read.delim(nuc_file, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE)
  if (nrow(nuc) != nrow(te) || !identical(as.character(nuc[[4]]), te$id)) {
    stop("bedtools nuc output does not match the all-TE BED")
  }
  count_base <- function(base) {
    col <- grep(paste0("_num_", base, "$"), names(nuc))
    if (length(col) != 1) stop("Cannot find ", base, " counts in bedtools nuc output")
    as.numeric(nuc[[col]])
  }
  a <- count_base("A")
  c <- count_base("C")
  g <- count_base("G")
  t <- count_base("T")
  called <- a + c + g + t

  per_te <- data.frame(
    TE_ID=te$id,
    Chromosome=te$chrom,
    Start=as.integer(te$start),
    End=as.integer(te$end),
    Embedded_TE=te$embedded,
    GC_content_pct=ifelse(called > 0, 100 * (c + g) / called, NA_real_),
    TE_length_bp=te$end - te$start,
    TE_TSS_distance_bp=distance,
    stringsAsFactors=FALSE
  )
  group_masks <- list("All TEs"=rep(TRUE, nrow(per_te)),
                      "Promoter-embedded TEs"=per_te$Embedded_TE)
  for (group in names(group_ids)) {
    missing_ids <- setdiff(group_ids[[group]], te$id)
    if (length(missing_ids) > 0) stop(group, " IDs absent from all-TE BED: ", paste(head(missing_ids), collapse=", "))
    group_masks[[group]] <- te$id %in% group_ids[[group]]
    per_te[[gsub("[^A-Za-z0-9]+", "_", group)]] <- group_masks[[group]]
  }
  write.table(per_te, file.path(outdir, paste0(output_prefix, "_TE_features_per_TE.txt")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")

  feature_labels <- c(
    GC_content_pct="GC content (%)",
    TE_length_bp="TE length (bp)",
    TE_TSS_distance_bp="TE-TSS distance (bp)"
  )
  rows <- list()
  for (feature in names(feature_labels)) {
    for (group in names(group_masks)) {
      values <- per_te[[feature]][group_masks[[group]]]
      values <- values[is.finite(values)]
      n <- length(values)
      sd_value <- if (n > 1) sd(values) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        Group=group, Feature=unname(feature_labels[[feature]]), N=n,
        Mean=if (n > 0) mean(values) else NA_real_,
        SD=sd_value, SEM=sd_value / sqrt(n),
        stringsAsFactors=FALSE
      )
    }
  }
  summary <- do.call(rbind, rows)
  write.table(summary, file.path(outdir, paste0(output_prefix, "_TE_features_summary.txt")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")

  plot_files <- c(
    GC_content_pct=paste0(output_prefix, "_TE_GC_content_bar.png"),
    TE_length_bp=paste0(output_prefix, "_TE_length_bar.png"),
    TE_TSS_distance_bp=paste0(output_prefix, "_TE_TSS_distance_bar.png")
  )
  if (is.null(group_colors)) {
    group_colors <- c("All TEs"="#A6A6A6", "Promoter-embedded TEs"="#CCD881",
                      "Positive TEs"="#638E77", "Negative TEs"="#8A79A5")
  }
  for (feature in names(feature_labels)) {
    plot_df <- summary[summary$Feature == feature_labels[[feature]], , drop=FALSE]
    plot_df$Group <- factor(plot_df$Group, levels=names(group_masks))
    plot_df$PlotMean <- ifelse(is.finite(plot_df$Mean), plot_df$Mean, 0)
    error_df <- plot_df[is.finite(plot_df$SEM), , drop=FALSE]
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x=Group, y=PlotMean, fill=Group)) +
      ggplot2::geom_col(width=0.62, colour="black", linewidth=0.4) +
      ggplot2::geom_errorbar(data=error_df,
        ggplot2::aes(ymin=pmax(0, Mean - SEM), ymax=Mean + SEM),
        width=0.18, linewidth=0.7) +
      ggplot2::geom_text(ggplot2::aes(y=PlotMean + ifelse(is.na(SEM), 0, SEM), label=paste0("n=", N)),
                         vjust=-0.6, size=4.5) +
      ggplot2::scale_fill_manual(values=group_colors) +
      ggplot2::scale_x_discrete(labels=c("All TEs"="All TEs", "Promoter-embedded TEs"="Promoter-\nembedded TEs",
                                         "Positive TEs"="Positive\nTEs", "Negative TEs"="Negative\nTEs")) +
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0, 0.2))) +
      ggplot2::labs(title=unname(feature_labels[[feature]]), x=NULL,
                    y=paste("Mean", feature_labels[[feature]])) +
      ggplot2::theme_classic(base_size=18) +
      ggplot2::theme(
        legend.position="none",
        axis.line.x.bottom=ggplot2::element_line(colour="black", linewidth=0.6),
        axis.line.y.left=ggplot2::element_line(colour="black", linewidth=0.6),
        axis.text.x=ggplot2::element_text(face="bold", size=if (length(group_masks) > 2) 12 else 15),
        plot.title=ggplot2::element_text(face="bold", hjust=0.5)
      )
    ggplot2::ggsave(file.path(outdir, plot_files[[feature]]), p,
                    width=6.5, height=5.5, dpi=400, bg="white")
  }
  invisible(per_te)
}

plot_te_family_composition <- function(per_te, family_file, outdir) {
  if (!requireNamespace("RColorBrewer", quietly=TRUE)) {
    stop("RColorBrewer is required for TE family composition")
  }
  raw_family <- read.table(family_file, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
  if (ncol(raw_family) == 2) {
    family <- raw_family[, 1:2]
  } else if (ncol(raw_family) >= 7) {
    family <- raw_family[, c(1, 7)]
  } else {
    stop("TE family file needs two columns (TE ID, family) or at least seven columns")
  }
  names(family) <- c("TE_ID", "TE_family")
  if (anyDuplicated(family$TE_ID)) stop("TE family file contains duplicate TE IDs")

  matched <- merge(per_te[, c("TE_ID", "Embedded_TE", "Positive_TEs", "Negative_TEs")],
                   family, by="TE_ID")
  if (nrow(matched) == 0) stop("No TE IDs match the TE family file")
  if (nrow(matched) < nrow(per_te)) {
    warning(nrow(per_te) - nrow(matched), " TEs lack family annotations and are excluded from the family plot")
  }

  groups <- list("All TEs"=rep(TRUE, nrow(matched)),
                 "Promoter-embedded TEs"=matched$Embedded_TE,
                 "Positive TEs"=matched$Positive_TEs,
                 "Negative TEs"=matched$Negative_TEs)
  families <- sort(unique(matched$TE_family))
  counts <- do.call(cbind, lapply(groups, function(mask) {
    as.numeric(table(factor(matched$TE_family[mask], levels=families)))
  }))
  dimnames(counts) <- list(families, names(groups))
  totals <- colSums(counts)
  if (any(totals[1:2] == 0)) stop("All-TE and embedded-TE groups must contain annotated families")
  percentages <- sweep(counts, 2, ifelse(totals > 0, totals, 1), "/") * 100
  enrichment <- log2(percentages[, 2] / percentages[, 1])
  family_order <- families[order(enrichment, decreasing=TRUE)]

  result <- data.frame("TE family"=family_order, check.names=FALSE)
  for (group in names(groups)) {
    result[[group]] <- counts[family_order, group]
    result[[paste0(group, " (%)")]] <- percentages[family_order, group]
  }
  for (group in names(groups)[-1]) {
    result[[paste0("Log2 enrichment (", group, " vs All TEs)")]] <-
      log2(percentages[family_order, group] / percentages[family_order, "All TEs"])
  }
  write.table(result, file.path(outdir, "OUTPUT_5_Promoter_embedded_TE_family.txt"),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")

  plot_data <- do.call(rbind, lapply(seq_along(groups), function(i) {
    data.frame(Group=names(groups)[i], TE_family=family_order,
               x=i, Percentage=percentages[family_order, i])
  }))
  plot_data$TE_family <- factor(plot_data$TE_family, levels=family_order)
  colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(family_order))
  names(colors) <- family_order
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x=x, y=Percentage, fill=TE_family)) +
    ggplot2::geom_col(width=0.45, colour="black", linewidth=0.3) +
    ggplot2::scale_x_continuous(breaks=1:4,
      labels=c("All TEs", "Promoter-\nembedded TEs", "Positive\nTEs", "Negative\nTEs"),
      limits=c(0.55, 4.45), expand=c(0, 0)) +
    ggplot2::scale_y_continuous(breaks=c(0, 25, 50, 75, 100), expand=c(0, 0)) +
    ggplot2::coord_cartesian(ylim=c(0, 100.5)) +
    ggplot2::scale_fill_manual(values=colors, name="TE Types") +
    petem_theme_classic() +
    ggplot2::labs(title="TE family composition", x=NULL, y="Percentage (%)") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=14, face="bold"),
                   axis.title.y=ggplot2::element_text(size=16),
                   axis.line.x.bottom=ggplot2::element_line(colour="black", linewidth=0.6),
                   axis.line.y.left=ggplot2::element_line(colour="black", linewidth=0.6),
                   legend.key.size=grid::unit(1.2, "lines"))
  ggplot2::ggsave(file.path(outdir, "OUTPUT_5_Promoter_embedded_TE_family_enrichment.png"),
                  p, width=10.5, height=6.3, dpi=400, bg="white")
  invisible(result)
}
