# Rscript 5_cross_condition_correlation.R --module0-dir module_0 --unexp "$unexp"

start_time <- Sys.time()

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
source(file.path(dirname(normalizePath(script_path)), "plot_defaults.R"))

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rlang))


#---- Option parser ----
raw_args <- commandArgs(trailingOnly=TRUE)
stage_pair_opt <- character(0)
stage_i <- which(raw_args %in% c("-stage", "--stage"))
if (length(stage_i) > 1) stop("Use -stage only once.")
if (length(stage_i) == 1) {
  i <- stage_i[1]
  if (length(raw_args) < i + 2) stop("-stage requires two values, e.g. -stage leaf root")
  stage_pair_opt <- raw_args[(i + 1):(i + 2)]
  raw_args <- raw_args[-c(i, i + 1, i + 2)]
}
option_list = list(
  make_option(c("--DEG"), type="character", help="Gene expression DEG file"),
  make_option(c("--DETE"), type="character", help="TE expression DETE file"),
  make_option(c("--module0-dir"), type="character", default=".", help="Directory containing Module 0 outputs such as Tab_* and TE_overlap_promoter.bed"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--unexp"), type="character", default="n", help="Include unexpressed TEs? (y/n)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args=raw_args)

include_unexp <- tolower(opt$unexp) == "y"
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
module0_file <- function(name) file.path(opt$`module0-dir`, name)
output_file <- function(name) file.path(opt$outdir, name)
if (is.null(opt$DEG) || opt$DEG == "") opt$DEG <- module0_file("DEG.txt")
if (is.null(opt$DETE) || opt$DETE == "") opt$DETE <- module0_file("DETE.txt")


#---- Helper functions ----
quadrant_counts <- function(df, xx, yy) {
  q1 <- sum(df[[xx]] > 0 & df[[yy]] > 0, na.rm=TRUE)
  q2 <- sum(df[[xx]] < 0 & df[[yy]] > 0, na.rm=TRUE)
  q3 <- sum(df[[xx]] < 0 & df[[yy]] < 0, na.rm=TRUE)
  q4 <- sum(df[[xx]] > 0 & df[[yy]] < 0, na.rm=TRUE)
  c(Q1=q1, Q2=q2, Q3=q3, Q4=q4)
}

lm_eqn_text <- function(df, yy, xx){
  m <- lm(df[[yy]] ~ df[[xx]], df)
  a <- format(coef(m)[1], digits=2)
  b <- format(coef(m)[2], digits=2)
  r2 <- format(summary(m)$r.squared, digits=3)
  eq <- paste0("y = ", a, " + ", b, " * x", ", R^2 = ", r2)
  return(eq)
}

plot_delta_scatter <- function(df, x, y, fname, title, xlab, ylab){
  counts <- quadrant_counts(df, x, y)
  eq_text <- paste0(lm_eqn_text(df, y, x))
  xlim <- max(abs(df[[x]]), na.rm=TRUE)
  ylim <- max(abs(df[[y]]), na.rm=TRUE)

  # ---- Spearman correlation ----
  cor_res <- cor.test(df[[x]], df[[y]], method="spearman", exact=FALSE)
  pval_text <- paste0(
    "Spearman rho = ", round(cor_res$estimate, 3),
    ", p = ", signif(cor_res$p.value, 3)
  )

  caption_text <- paste(eq_text, pval_text, sep="\n")

  png(file=fname, width=3000, height=2500, res=400)
  xvar <- x; yvar <- y
  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(color="#6E6E6E", alpha=0.45, size=1.6) +
    geom_hline(yintercept=0, linetype="dashed", color="gray70") +
    geom_vline(xintercept=0, linetype="dashed", color="gray70") +
    geom_smooth(method="lm", se=FALSE, col="black", formula = y ~ x) +
    scale_y_continuous(limits=c(-ylim, ylim)) +
    scale_x_continuous(limits=c(-xlim, xlim)) +
    petem_theme_bw() +
    annotate("text", x=xlim, y=ylim, label=paste("Q1:", counts["Q1"]), size=PETEM_ANNOTATION_TEXT_SIZE+2, hjust=1, fontface="bold") +
    annotate("text", x=-xlim, y=ylim, label=paste("Q2:", counts["Q2"]), size=PETEM_ANNOTATION_TEXT_SIZE+2, hjust=0, fontface="bold") +
    annotate("text", x=-xlim, y=-ylim, label=paste("Q3:", counts["Q3"]), size=PETEM_ANNOTATION_TEXT_SIZE+2, hjust=0, fontface="bold") +
    annotate("text", x=xlim, y=-ylim, label=paste("Q4:", counts["Q4"]), size=PETEM_ANNOTATION_TEXT_SIZE+2, hjust=1, fontface="bold") +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.caption=element_text(size=PETEM_AXIS_TEXT_SIZE+2, hjust=0)
    ) +
    labs(title=title, caption = caption_text, x=xlab, y=ylab)
  suppressMessages(print(p))

  dev.off()
}

#---- Read expression ----
gene_exp<-read.table(opt$DEG, row.names=1, header=T, sep="\t") 
TE_exp<-read.table(opt$DETE, row.names=1, header=T, sep="\t") 

#---- filter unexpressed ----
gene_stage_cols <- setdiff(colnames(gene_exp), grep("^(logFC_|PValue_|FDR_)", colnames(gene_exp), value=TRUE))
TE_stage_cols   <- setdiff(colnames(TE_exp), grep("^(logFC_|PValue_|FDR_)", colnames(TE_exp), value=TRUE))

if (length(gene_stage_cols) != 2 || length(TE_stage_cols) != 2) {
  stop("Module 5 requires exactly two expression stages/conditions in both DEG and DETE inputs.")
}
if (!setequal(gene_stage_cols, TE_stage_cols)) {
  stop("Module 5 requires DEG and DETE inputs to contain the same two expression stages/conditions.")
}

gene_exp <- gene_exp[rowSums(gene_exp[, gene_stage_cols, drop=FALSE]) != 0, , drop=FALSE]

if(!include_unexp){
  TE_exp <- TE_exp[rowSums(TE_exp[, TE_stage_cols, drop=FALSE]) != 0, , drop=FALSE]
}


#---- Rename columns ----
colnames(gene_exp) <- gsub("^logFC_", "dexp_", colnames(gene_exp))
colnames(TE_exp)   <- gsub("^logFC_", "dTEexp_", colnames(TE_exp))
colnames(gene_exp) <- gsub("^PValue_", "PV_g_", colnames(gene_exp))
colnames(TE_exp)   <- gsub("^PValue_", "PV_TE_", colnames(TE_exp))
colnames(gene_exp) <- gsub("^FDR_", "FDR_g_", colnames(gene_exp))
colnames(TE_exp)   <- gsub("^FDR_", "FDR_TE_", colnames(TE_exp))

#---- Read methylation ----
CG_TE  <- read.table(module0_file("Tab_TE_CG.txt"), header=TRUE, sep="\t")
CHG_TE <- read.table(module0_file("Tab_TE_CHG.txt"), header=TRUE, sep="\t")
CHH_TE <- read.table(module0_file("Tab_TE_CHH.txt"), header=TRUE, sep="\t")

methylation_stage_cols <- Reduce(intersect, list(colnames(CG_TE)[-1], colnames(CHG_TE)[-1], colnames(CHH_TE)[-1]))

# follow order of stages in delta expression files
fc_cols <- grep("^dexp_", colnames(gene_exp), value=TRUE)
if (length(fc_cols) == 0) {
  stop("No differential expression comparison columns detected in DEG file after renaming (expected source columns starting with 'logFC_').")
}
if (length(stage_pair_opt) == 0 && length(fc_cols) != 1) {
  stop("Module 5 requires exactly one logFC comparison in the DEG input.")
}

stage_pairs <- str_match(fc_cols, "^dexp_([^_]+)_([^_]+)$")
stage_pairs <- stage_pairs[,2:3, drop=FALSE]
valid_fc <- complete.cases(stage_pairs)
stage_pairs <- stage_pairs[valid_fc, , drop=FALSE]
fc_cols <- fc_cols[valid_fc]
if (nrow(stage_pairs) == 0) {
  stop("Unable to parse stage comparisons from logFC column names.")
}
if (length(stage_pair_opt) > 0) {
  keep <- stage_pairs[,1] == stage_pair_opt[1] & stage_pairs[,2] == stage_pair_opt[2]
  fc_cols <- fc_cols[keep]
  stage_pairs <- stage_pairs[keep, , drop=FALSE]
  if (length(fc_cols) != 1) stop("Requested -stage comparison not found in DEG/DETE.")
}
if (nrow(stage_pairs) != 1) {
  stop("Module 5 requires exactly one stage pair comparison.")
}

if (!all(stage_pairs[1, ] %in% gene_stage_cols)) {
  stop("Module 5 requires the logFC comparison to match the same two expression stages/conditions in DEG.")
}

te_fc_cols <- grep("^dTEexp_", colnames(TE_exp), value=TRUE)
if (length(stage_pair_opt) > 0) {
  te_fc_cols <- paste0("dTEexp_", stage_pair_opt[1], "_", stage_pair_opt[2])
}
if (length(te_fc_cols) != 1 || !te_fc_cols %in% colnames(TE_exp)) {
  stop("Module 5 requires exactly one logFC comparison in the DETE input.")
}

te_stage_pairs <- str_match(te_fc_cols, "^dTEexp_([^_]+)_([^_]+)$")
te_stage_pairs <- te_stage_pairs[,2:3, drop=FALSE]
te_stage_pairs <- te_stage_pairs[complete.cases(te_stage_pairs), , drop=FALSE]
if (nrow(te_stage_pairs) != 1 || !identical(as.vector(te_stage_pairs[1, ]), as.vector(stage_pairs[1, ]))) {
  stop("Module 5 requires DEG and DETE inputs to use the same single stage pair comparison.")
}

valid_pair <- apply(stage_pairs, 1, function(pair) all(pair %in% methylation_stage_cols))
if (!all(valid_pair)) {
  stop(sprintf("Missing methylation data for comparison(s): %s",
               paste(apply(stage_pairs[!valid_pair, , drop=FALSE], 1, paste, collapse=" vs "), collapse=", ")))
}
stage_pairs <- stage_pairs[valid_pair, , drop=FALSE]
fc_cols <- fc_cols[valid_pair]
if (nrow(stage_pairs) == 0) {
  stop("No stage comparisons have matching methylation columns in Tab_TE_* tables.")
}
if (length(methylation_stage_cols) != 2 || !setequal(methylation_stage_cols, gene_stage_cols)) {
  stop("Module 5 requires methylation tables to contain exactly the same two stages/conditions as DEG and DETE.")
}

# calculate delta methylation
for(k in seq_len(nrow(stage_pairs))){
  si <- stage_pairs[k,1]
  sj <- stage_pairs[k,2]
  
  CG_TE[[paste0("dTECG_", si, "_", sj)]]   <- (CG_TE[[si]] - CG_TE[[sj]]) * 100
  CHG_TE[[paste0("dTECHG_", si, "_", sj)]] <- (CHG_TE[[si]] - CHG_TE[[sj]]) * 100
  CHH_TE[[paste0("dTECHH_", si, "_", sj)]] <- (CHH_TE[[si]] - CHH_TE[[sj]]) * 100
}

#---- Merge with promoter overlaps ----
ins_promoter <- read.table(module0_file("TE_overlap_promoter.bed"), header=FALSE)
ins_promoter2 <- ins_promoter[,c("V4","V10")]

ins3 <- merge(gene_exp, ins_promoter2, by.x="row.names", by.y="V10")
colnames(ins3)[1] <- "gene_id"
ins4 <- merge(ins3, TE_exp, by.x="V4", by="row.names")
colnames(ins4)[colnames(ins4)=="V4"] <- "TE_id"

ins_CG  <- merge(ins4, CG_TE,  by.x="TE_id", by.y="ID")
ins_CHG <- merge(ins4, CHG_TE, by.x="TE_id", by.y="ID")
ins_CHH <- merge(ins4, CHH_TE, by.x="TE_id", by.y="ID")

#---- function for Q2/Q4 boxplots ----
plot_gene_TE_box <- function(df_subset, mC_type, si, sj, mode="Q2"){
  if(nrow(df_subset) == 0) return(NULL)
  
  # output gene-TE pair
  out_tbl <- df_subset[,c("gene_id","TE_id")]
  write.table(out_tbl,
              file=output_file(paste0("OUTPUT_5_",mode,"_gene_TE_pairs_", mC_type, "_", si, "_", sj,".txt")),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  # expression
  expr_mat <- df_subset[,c("gene_id", paste0(si,".x"), paste0(sj,".x"))]
  expr_mat_long <- melt(expr_mat, id.vars="gene_id", variable.name="stage", value.name="expr")
  expr_mat_long$expr <- log2(expr_mat_long$expr + 0.1)
  expr_mat_long$stage <- gsub("\\.x$", "", expr_mat_long$stage)
  expr_mat_long$stage <- factor(expr_mat_long$stage, levels=c(sj, si))
  
  # methylation
  meth_mat <- df_subset[,c("TE_id", si, sj)]
  meth_mat_long <- melt(meth_mat, id.vars="TE_id", variable.name="stage", value.name="methyl")
  meth_mat_long$methyl <- meth_mat_long$methyl * 100  
  meth_mat_long$stage <- factor(meth_mat_long$stage, levels=c(sj, si))
  
  expr_range <- range(expr_mat_long$expr, na.rm=TRUE)
  meth_range <- range(meth_mat_long$methyl, na.rm=TRUE)
  scale_factor <- diff(expr_range) / diff(meth_range)
  
  expr_means <- aggregate(expr ~ stage, expr_mat_long, mean)
  meth_means <- aggregate(methyl ~ stage, meth_mat_long, mean)
  
  p <- ggplot() +
    geom_boxplot(data=expr_mat_long, aes(x=stage, y=expr),
                 fill="#BFBFBF", width=0.35, position=position_nudge(x=-0.2)) +
    geom_point(data=expr_means, aes(x=stage, y=expr), shape=18, size=PETEM_ANNOTATION_TEXT_SIZE + 2, color="black",
               position=position_nudge(x=-0.2)) +
    geom_text(data=expr_means, aes(x=stage, y=expr, label=sprintf("%.1f", expr)),
              size=PETEM_ANNOTATION_TEXT_SIZE, vjust=-1.2, color="black", position=position_nudge(x=-0.2)) +
    geom_boxplot(data=meth_mat_long, aes(x=stage, y=methyl*scale_factor + expr_range[1] - meth_range[1]*scale_factor),
                 fill="#E5C1AF", width=0.35, position=position_nudge(x=0.2)) +
    geom_point(data=meth_means, aes(x=stage, y=methyl*scale_factor + expr_range[1] - meth_range[1]*scale_factor),
               shape=18, size=PETEM_ANNOTATION_TEXT_SIZE + 2, color="black", position=position_nudge(x=0.2)) +
    geom_text(data=meth_means, aes(x=stage, y=methyl*scale_factor + expr_range[1] - meth_range[1]*scale_factor,
                                   label=sprintf("%.1f", methyl)),
              size=PETEM_ANNOTATION_TEXT_SIZE+2, vjust=-1.2, color="black", position=position_nudge(x=0.2)) +
    scale_y_continuous(
      name="Expression (log2 RPKM)",
      sec.axis = sec_axis(
        ~ (. - expr_range[1] + meth_range[1]*scale_factor)/scale_factor,
        name=paste0(mC_type, " methylation level (%)")
      )
    ) +
    scale_x_discrete(name="Stage") +
    petem_theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y.right = element_text(angle=90)
    )
  
  ggsave(filename=output_file(paste0("OUTPUT_5_",mode,"_boxplot_", mC_type, "_", si, "_", sj,".png")), p, width=6, height=5)
  print(p)
}

#---- Scatter plots + Q2/Q4 boxplots ----
for(k in seq_along(fc_cols)){
  x_col <- fc_cols[k]
  si <- stage_pairs[k,1]
  sj <- stage_pairs[k,2]
  pv_col <- paste0("FDR_g_", si, "_", sj)

  # gene exp vs TE exp
  df_sub <- subset(ins4,
    abs(ins4[[x_col]]) > 1 & ins4[[pv_col]] < 0.05
  )
  if(nrow(df_sub) > 0){
    plot_delta_scatter(df_sub,
      x_col, paste0("dTEexp_", si, "_", sj),
      output_file(paste0("OUTPUT_5_geneexp_TEexp_change_", si, "_", sj, "_scatter.png")),
      paste0("TE and gene expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"Gene expression level (log2 RPKM FC)"),
      expression(Delta~"TE expression level (log2 RPKM FC)")
    )
  }

  # gene exp vs TE mC + Q2/Q4
  for(mC_type in c("CG","CHG","CHH")){
    df_TE <- get(paste0("ins_", mC_type))
    y_col <- paste0("dTE", mC_type, "_", si, "_", sj)
    df_sub <- subset(df_TE,
      abs(df_TE[[x_col]]) > 1 & df_TE[[pv_col]] < 0.05
    )
    if(nrow(df_sub) > 0){
      # scatter plot
      plot_delta_scatter(df_sub,
        x_col, y_col,
        output_file(paste0("OUTPUT_5_geneexp_TEm", mC_type, "_change_", si, "_", sj, "_scatter.png")),
        paste0("TE m", mC_type, " and gene expression changes\nbetween ", si, " and ", sj),
        expression(Delta~"Gene expression level (log2 RPKM FC)"),
        expression(Delta~paste("TE methylation level (%)"))
      )
      # Q2 boxplot
      df_Q2 <- subset(df_TE,
        abs(df_TE[[x_col]]) > 1 &
        df_TE[[pv_col]] < 0.05 &
        df_TE[[x_col]] < 0 &
        df_TE[[y_col]] > 0
      )
      plot_gene_TE_box(df_Q2, mC_type, si, sj, mode="Q2")
      # Q4 boxplot
      df_Q4 <- subset(df_TE,
        abs(df_TE[[x_col]]) > 1 &
        df_TE[[pv_col]] < 0.05 &
        df_TE[[x_col]] > 0 &
        df_TE[[y_col]] < 0
      )
      plot_gene_TE_box(df_Q4, mC_type, si, sj, mode="Q4")
    }
  }

  # TE exp vs TE mC
  for(mC_type in c("CG","CHG","CHH")){
    df_TE <- get(paste0("ins_", mC_type))
    x_TE_col <- paste0("dTEexp_", si, "_", sj)
    y_TE_col <- paste0("dTE", mC_type, "_", si, "_", sj)
    pv_TE_col <- paste0("FDR_TE_", si, "_", sj)

    df_sub <- subset(df_TE,
      abs(df_TE[[x_TE_col]]) > 1 & df_TE[[pv_TE_col]] < 0.05
    )
    if(nrow(df_sub) > 0){
      plot_delta_scatter(df_sub,
        x_TE_col, y_TE_col,
        output_file(paste0("OUTPUT_5_TEexp_TEm", mC_type, "_change_", si, "_", sj, "_scatter.png")),
        paste0("TE m", mC_type, " and TE expression changes\nbetween ", si, " and ", sj),
        expression(Delta~"TE expression level (log2 RPKM FC)"),
        expression(Delta~paste("TE methylation level (%)"))
      )
    }
  }
}

end_time <- Sys.time()
print(end_time-start_time)
