#!/usr/bin/env Rscript
# Rscript 4_single_condition_correlation.R --eg expression_gene.txt --et expression_TE.txt --unexp y/n \
#   --smooth 3 --ylim_CG 40 --ylim_CHG 5 --ylim_CHH 5 --ylim_TEexpTEmC_CH 8 --ylim_TEexpTEmC_CG 40


start_time <- Sys.time()

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
source(file.path(dirname(normalizePath(script_path)), "plot_defaults.R"))

suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

MODULE4_THEME_BASE_SIZE <- PETEM_THEME_BASE_SIZE + 2
MODULE4_CEX_AXIS <- PETEM_BASE_CEX_AXIS * 1.1
MODULE4_CEX_LABEL <- PETEM_BASE_CEX_LABEL * 1.1
MODULE4_CEX_NOTE <- PETEM_BASE_CEX_NOTE * 1.1
MODULE4_LINE_LABEL_SIZE <- 4.5
MODULE4_BAR_LABEL_SIZE <- 4.4
MODULE4_CG_COL <- "#CFA699"
MODULE4_CHG_COL <- "#71A061"
MODULE4_CHH_COL <- "#3871A6"
MODULE4_FRAME_WIDTH <- 1
MODULE4_GRID_COL <- "gray92"
plot.title = element_text(
    size = 24,
    face = "bold",
    hjust = 0.5
)

module4_theme_bw <- function() {
  petem_theme_bw(base_size = MODULE4_THEME_BASE_SIZE) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = MODULE4_FRAME_WIDTH),
      panel.grid.major = element_line(color = MODULE4_GRID_COL, linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

#---- Option parser ----
option_list = list(
  make_option(c("--eg"), type="character", help="Gene expression file"),
  make_option(c("--et"), type="character", help="TE expression file"),
  make_option(c("--module0-dir"), type="character", default=".", help="Directory containing Module 0 outputs such as Tab_* and TE_overlap_promoter.bed"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--unexp"), type="character", default="n", help="Include unexpressed TEs for sliding plots? (y/n)"),
  make_option(c("--smooth"), type="numeric", default=3, help="sliding window multiplier, 1-5 (default=3)"),
  make_option(c("--ylim_CG"), type="numeric", default=50, help="ylim for gene exp vs TE/promoter mC, CG (default=50)"),
  make_option(c("--ylim_CHG"), type="numeric", default=10, help="ylim for gene exp vs TE/promoter mC, CHG (default=10)"),
  make_option(c("--ylim_CHH"), type="numeric", default=10, help="ylim for gene exp vs TE/promoter mC, CHH (default=10)"),
  make_option(c("--ylim_TEexpTEmC_CH"), type="numeric", default=15, help="ylim for TE exp vs TE mC, CH (default=15)"),
  make_option(c("--ylim_TEexpTEmC_CG"), type="numeric", default=30, help="ylim for TE exp vs TE mC, CG (default=30)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (opt$smooth < 1 || opt$smooth > 5) stop("--smooth must be between 1 and 5.")

include_unexp <- tolower(opt$unexp) == "y"
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

module0_file <- function(name) file.path(opt$`module0-dir`, name)
output_file <- function(name) file.path(opt$outdir, name)

#---- Functions ----


sort_exp <- function(df, stage_exp){
  df_stage = df[complete.cases(df),]
  stage = df_stage[order(df_stage[,stage_exp]),]
  return(stage)
}

sliding <- function(vec){
  lth <- length(vec)
  step <- floor(lth / 100)
  step <- max(1, step)
  window <- lth - (100 - 1) * step
  window <- max(1, window)
  roll_mean <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm=TRUE)
  }
  return(as.vector(rollapply(zoo(vec), width=window, by=step, FUN=roll_mean, align="center")))
}

calc_window <- function(df){
  step <- max(1, floor(nrow(df)/100))
  window <- max(1, step*opt$smooth)
  window <- min(window, nrow(df))
  return(window)
}

collapse_methylation_rows <- function(df, id_col, value_cols) {
  if (nrow(df) == 0) return(df)
  mean_na <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm=TRUE)
  }
  out <- aggregate(df[, value_cols, drop=FALSE],
                   by=list(df[[id_col]]),
                   FUN=mean_na)
  colnames(out)[1] <- id_col
  out
}

expand_ylim <- function(requested, values, pad=1.08) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(requested)
  }
  max(requested, max(values, na.rm=TRUE) * pad)
}

run_cor_test <- function(x, y, method="pearson") {
  if (method == "spearman") {
    return(cor.test(x, y, method=method, exact=FALSE))
  }
  cor.test(x, y, method=method)
}

corr_mC <- function(df, exp, CG, CHG, CHH, method="pearson") {
  log2exp<-as.numeric(log2(df[[exp]]))
  cg_corr <- run_cor_test(log2exp, df[[CG]], method=method)
  chg_corr <- run_cor_test(log2exp, df[[CHG]], method=method)
  chh_corr <- run_cor_test(log2exp, df[[CHH]], method=method)
  data.frame(
    Corr = c(
      cg_corr$estimate,
      chg_corr$estimate,
      chh_corr$estimate
    ),
    Pvalue = c(
      cg_corr$p.value,
      chg_corr$p.value,
      chh_corr$p.value
    ),
    row.names = c("CG", "CHG", "CHH")
  )
}

corr_exp <- function(df, x, y, method="pearson") {
  corr <- run_cor_test(as.numeric(log2(df[[x]])), as.numeric(log2(df[[y]])), method=method)
  data.frame(Corr=corr$estimate, Pvalue=corr$p.value, row.names=paste(x, y, sep="_vs_"))
}

plot_corr_bar <- function(stage, gdf, gdfexp, TEdf, method="pearson") {
  cg_col <- MODULE4_CG_COL
  chg_col <- MODULE4_CHG_COL
  chh_col <- MODULE4_CHH_COL
  te_gray <- "#9E9E9E"

  # correlations
  TEexp_TEmC_corr <- corr_mC(TEdf, paste0(stage,"_TEexp"), paste0(stage,"_TE_CG"), paste0(stage,"_TE_CHG"), paste0(stage,"_TE_CHH"), method=method)
  exp_TEmC_corr   <- corr_mC(gdf, paste0(stage,"_exp"), paste0(stage,"_TE_CG"), paste0(stage,"_TE_CHG"), paste0(stage,"_TE_CHH"), method=method)
  exp_TEexp_corr  <- corr_exp(gdfexp, paste0(stage,"_exp"), paste0(stage,"_TEexp"), method=method)
  
  # combine
  corr_df <- data.frame(
    Comparison = c("TEexp vs TEmCG", "TEexp vs TEmCHG", "TEexp vs TEmCHH",
                   "geneexp vs TEexp", "geneexp vs TEmCG", "geneexp vs TEmCHG", "geneexp vs TEmCHH"),
    Corr = c(
      TEexp_TEmC_corr["CG","Corr"],
      TEexp_TEmC_corr["CHG","Corr"],
      TEexp_TEmC_corr["CHH","Corr"],
      exp_TEexp_corr[1,"Corr"],
      exp_TEmC_corr["CG","Corr"],
      exp_TEmC_corr["CHG","Corr"],
      exp_TEmC_corr["CHH","Corr"]
    ),
    Color = c(cg_col, chg_col, chh_col, te_gray, cg_col, chg_col, chh_col)
  )
  corr_df$Comparison <- factor(corr_df$Comparison, levels=rev(corr_df$Comparison))
  group_axis_labels <- c(
    "geneexp vs TEmCHH" = "",
    "geneexp vs TEmCHG" = "Gene expression\nvs TE methylation",
    "geneexp vs TEmCG" = "",
    "geneexp vs TEexp" = "Gene expression\nvs TE expression",
    "TEexp vs TEmCHH" = "",
    "TEexp vs TEmCHG" = "TE expression\nvs TE methylation",
    "TEexp vs TEmCG" = ""
  )
  te_row <- which(levels(corr_df$Comparison) == "geneexp vs TEexp")
  
  # plot
  xmin <- min(corr_df$Corr, na.rm=TRUE) - 0.1
  xmax <- max(corr_df$Corr, na.rm=TRUE) + 0.1
  right_labels <- data.frame(
    Comparison = factor(c("TEexp vs TEmCG", "TEexp vs TEmCHG", "TEexp vs TEmCHH"),
                        levels = levels(corr_df$Comparison)),
    label = c("CG", "CHG", "CHH"),
    color = c(cg_col, chg_col, chh_col),
    x = max(0.03, 0.04 * max(abs(c(xmin, xmax)))),
    stringsAsFactors = FALSE
  )
  outfile <- output_file(paste0("OUTPUT_4_",method,"_correlation_bar_",stage,".png"))
  png(file=outfile,width=2000,height=1500,res=400)
  print(
    ggplot(corr_df, aes(x=Corr,y=Comparison,fill=Color)) +
      geom_bar(stat="identity", width=0.55) +
      geom_hline(yintercept=c(te_row - 0.5, te_row + 0.5), color="black", linewidth=0.6) +
      geom_text(aes(label=round(Corr,2)),
                hjust=ifelse(corr_df$Corr>=0, -0.1, 1.1), size=MODULE4_BAR_LABEL_SIZE) +
      geom_text(data=right_labels, aes(x=x, y=Comparison, label=label, color=color),
                inherit.aes=FALSE, hjust=0, size=MODULE4_BAR_LABEL_SIZE + 0.2, fontface="bold") +
      scale_fill_identity() +
      scale_color_identity() +
      scale_y_discrete(labels = group_axis_labels) +
      xlim(xmin, xmax) +
      labs(  title = "Correlation strengths between\n TEs and genes",
      x=paste0(tools::toTitleCase(method)," correlation coefficient"), y=NULL) +
      module4_theme_bw() +
      theme(
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = PETEM_AXIS_TEXT_SIZE + 1.5,face = "bold",hjust = 0.5)
      )
  )
  dev.off()
}


#---- Read expression ----
gene_exp <- read.table(opt$eg, header=TRUE, sep="\t", row.names=1)
#gene_exp<-read.table("tttgene.txt", header=T, sep="\t", row.names=1) 
gene_exp <- gene_exp[rowSums(gene_exp) != 0, , drop=FALSE]  # drop unexpressed genes

TE_exp <- read.table(opt$et, header=TRUE, sep="\t", row.names=1)
#TE_exp<-read.table("tttTE.txt", header=T, sep="\t", row.names=1) 
#TE_exp <- TE_exp[rowSums(TE_exp) != 0, ]

#---- Determine stages ----
stages <- intersect(colnames(gene_exp), colnames(TE_exp))

#---- Read methylation ----
CG_TE <- read.table(module0_file("Tab_TE_CG.txt"), header=TRUE, sep="\t")
CHG_TE <- read.table(module0_file("Tab_TE_CHG.txt"), header=TRUE, sep="\t")
CHH_TE <- read.table(module0_file("Tab_TE_CHH.txt"), header=TRUE, sep="\t")

CG_gene <- read.table(module0_file("Tab_gene_CG.txt"), header=TRUE, sep="\t")
CHG_gene <- read.table(module0_file("Tab_gene_CHG.txt"), header=TRUE, sep="\t")
CHH_gene <- read.table(module0_file("Tab_gene_CHH.txt"), header=TRUE, sep="\t")

CG_promoter <- read.table(module0_file("Tab_promoter_CG.txt"), header=TRUE, sep="\t")
CHG_promoter <- read.table(module0_file("Tab_promoter_CHG.txt"), header=TRUE, sep="\t")
CHH_promoter <- read.table(module0_file("Tab_promoter_CHH.txt"), header=TRUE, sep="\t")

CG_promoterselves <- read.table(module0_file("Tab_promoterselves_CG.txt"), header=TRUE, sep="\t")
CHG_promoterselves <- read.table(module0_file("Tab_promoterselves_CHG.txt"), header=TRUE, sep="\t")
CHH_promoterselves <- read.table(module0_file("Tab_promoterselves_CHH.txt"), header=TRUE, sep="\t")
CG_promoterselves$ID  <- sub("_[0-9]+$", "", CG_promoterselves$ID)
CHG_promoterselves$ID <- sub("_[0-9]+$", "", CHG_promoterselves$ID)
CHH_promoterselves$ID <- sub("_[0-9]+$", "", CHH_promoterselves$ID)

methylation_stages <- Reduce(intersect, list(colnames(CG_TE)[-1], colnames(CHG_TE)[-1], colnames(CHH_TE)[-1],
                                             colnames(CG_promoterselves)[-1], colnames(CHG_promoterselves)[-1], colnames(CHH_promoterselves)[-1],
                                             colnames(CG_promoter)[-1], colnames(CHG_promoter)[-1], colnames(CHH_promoter)[-1]))
stages <- intersect(stages, methylation_stages)
if (length(stages) == 0) {
  stop("No overlapping stages between expression matrices and methylation tables.")
}
if (length(stages) != 2) {
  stop("Module 4 requires exactly two shared stages/conditions across expression and methylation inputs.")
}


ins_promoter <- read.table(module0_file("TE_overlap_promoter.bed"), header=FALSE)
ins_promoter2 <- ins_promoter[,c("V4","V10")]

#---- Loop over stages ----
for(stage in stages){

  missing_stage <- !all(c(stage %in% colnames(CG_TE),
                          stage %in% colnames(CHG_TE),
                          stage %in% colnames(CHH_TE),
                          stage %in% colnames(CG_promoterselves),
                          stage %in% colnames(CHG_promoterselves),
                          stage %in% colnames(CHH_promoterselves),
                          stage %in% colnames(CG_promoter),
                          stage %in% colnames(CHG_promoter),
                          stage %in% colnames(CHH_promoter)))
  if (missing_stage) {
    stop(sprintf("Stage %s is missing methylation columns required by module 4.", stage))
  }

  #---- Merge data ----
  insGene <- merge(ins_promoter2, gene_exp[,stage, drop=FALSE], by.x="V10", by.y="row.names", all.x=TRUE)
  colnames(insGene)[3] <- paste0(stage, "_exp")
  
  insGeneTE <- merge(insGene, TE_exp[,stage, drop=FALSE], by.x="V4", by.y="row.names", all.x=TRUE)
  colnames(insGeneTE)[4] <- paste0(stage, "_TEexp")
  

  # keep only current stage columns for TE methylation
  CG_TE_stage  <- CG_TE[, c("ID", stage), drop=FALSE]
  CHG_TE_stage <- CHG_TE[, c("ID", stage), drop=FALSE]
  CHH_TE_stage <- CHH_TE[, c("ID", stage), drop=FALSE]
  colnames(CG_TE_stage)[2]  <- paste0(stage, "_TE_CG")
  colnames(CHG_TE_stage)[2] <- paste0(stage, "_TE_CHG")
  colnames(CHH_TE_stage)[2] <- paste0(stage, "_TE_CHH")

  insGeneTEmC <- Reduce(function(x,y) merge(x,y,by.x="V4",by.y="ID",all.x=TRUE),
                        list(insGeneTE, CG_TE_stage, CHG_TE_stage, CHH_TE_stage))
  
  # Promoter with TE mC
  CG_promoterselves_stage  <- CG_promoterselves[, c("ID", stage), drop=FALSE]
  CHG_promoterselves_stage <- CHG_promoterselves[, c("ID", stage), drop=FALSE]
  CHH_promoterselves_stage <- CHH_promoterselves[, c("ID", stage), drop=FALSE]
  colnames(CG_promoterselves_stage)[2]  <- paste0(stage, "_promoterselves_CG")
  colnames(CHG_promoterselves_stage)[2] <- paste0(stage, "_promoterselves_CHG")
  colnames(CHH_promoterselves_stage)[2] <- paste0(stage, "_promoterselves_CHH")
  CG_promoterselves_stage <- collapse_methylation_rows(CG_promoterselves_stage, "ID", paste0(stage, "_promoterselves_CG"))
  CHG_promoterselves_stage <- collapse_methylation_rows(CHG_promoterselves_stage, "ID", paste0(stage, "_promoterselves_CHG"))
  CHH_promoterselves_stage <- collapse_methylation_rows(CHH_promoterselves_stage, "ID", paste0(stage, "_promoterselves_CHH"))
  
  insGenePromC <- Reduce(function(x,y) merge(x,y,by.x="V10",by.y="ID",all.x=TRUE),
                         list(insGeneTEmC, CG_promoterselves_stage, CHG_promoterselves_stage, CHH_promoterselves_stage))
 
  # Promoter without TE
  gene_exp_woTE <- gene_exp[!(rownames(gene_exp) %in% ins_promoter$V10), ,drop=FALSE]
  gene_exp_woTE$gene_id <- rownames(gene_exp_woTE)
  
  CG_promoter_stage  <- CG_promoter[, c("ID", stage), drop=FALSE]
  CHG_promoter_stage <- CHG_promoter[, c("ID", stage), drop=FALSE]
  CHH_promoter_stage <- CHH_promoter[, c("ID", stage), drop=FALSE]
  colnames(CG_promoter_stage)[2]  <- paste0(stage, "_promoter_CG")
  colnames(CHG_promoter_stage)[2] <- paste0(stage, "_promoter_CHG")
  colnames(CHH_promoter_stage)[2] <- paste0(stage, "_promoter_CHH")
  
  woTEGenePromC <- Reduce(function(x,y) merge(x,y,by.x="gene_id",by.y="ID",all.x=TRUE),
                          list(gene_exp_woTE, CG_promoter_stage, CHG_promoter_stage, CHH_promoter_stage))
  
  #---- Sorting & sliding ----
  TEdf <- sort_exp(insGeneTEmC[,c("V4","V10",paste0(stage,"_exp"),paste0(stage,"_TEexp"),
                                  paste0(stage,"_TE_CG"),paste0(stage,"_TE_CHG"),paste0(stage,"_TE_CHH"))],
                   paste0(stage,"_TEexp"))
  gdf <- sort_exp(insGeneTEmC[,c("V4","V10",paste0(stage,"_exp"),paste0(stage,"_TEexp"),
                                 paste0(stage,"_TE_CG"),paste0(stage,"_TE_CHG"),paste0(stage,"_TE_CHH"))],
                  paste0(stage,"_exp"))
  wo <- sort_exp(woTEGenePromC[,c("gene_id",paste0(stage),
                                  paste0(stage,"_promoter_CG"),paste0(stage,"_promoter_CHG"),paste0(stage,"_promoter_CHH"))],
                  paste0(stage))


  # remove unexpressed genes/TEs
  TEdf <- TEdf[(TEdf[[paste0(stage,"_exp")]] != 0) & (TEdf[[paste0(stage,"_TEexp")]] != 0), ]

  gdfexp <- gdf
  gdfexp  <- gdfexp[(gdfexp[[paste0(stage,"_exp")]] != 0) & (gdfexp[[paste0(stage,"_TEexp")]] != 0), ]

  wo   <- wo[(wo[[paste0(stage)]] != 0) , ]

  if(include_unexp){
    gdf  <- gdf[(gdf[[paste0(stage,"_exp")]] != 0), ]
  } else {
    gdf  <- gdf[(gdf[[paste0(stage,"_exp")]] != 0) & (gdf[[paste0(stage,"_TEexp")]] != 0), ]
  }

  pmC <- gdf[, c("V4", "V10", paste0(stage, "_exp")), drop=FALSE]
  pmC$.pair_order <- seq_len(nrow(pmC))
  pmC <- Reduce(function(x,y) merge(x,y,by.x="V10",by.y="ID",all.x=TRUE,sort=FALSE),
                list(pmC, CG_promoterselves_stage, CHG_promoterselves_stage, CHH_promoterselves_stage))
  pmC <- pmC[order(pmC$.pair_order), ]
  pmC$.pair_order <- NULL
  if (nrow(pmC) != nrow(gdf)) {
    stop("Promoter containing TEs rows must match TE-gene pairs rows.")
  }

  TE_CG <- sliding(TEdf[[paste0(stage,"_TE_CG")]]*100)
  TE_CHG <- sliding(TEdf[[paste0(stage,"_TE_CHG")]]*100)
  TE_CHH <- sliding(TEdf[[paste0(stage,"_TE_CHH")]]*100)
  
  gTE_CG <- sliding(gdf[[paste0(stage,"_TE_CG")]]*100)
  gTE_CHG <- sliding(gdf[[paste0(stage,"_TE_CHG")]]*100)
  gTE_CHH <- sliding(gdf[[paste0(stage,"_TE_CHH")]]*100)
  gTE_exp <- sliding(log2(gdfexp[[paste0(stage,"_TEexp")]]))
  
  pCG <- sliding(pmC[[paste0(stage,"_promoterselves_CG")]]*100)
  pCHG <- sliding(pmC[[paste0(stage,"_promoterselves_CHG")]]*100)
  pCHH <- sliding(pmC[[paste0(stage,"_promoterselves_CHH")]]*100)
  
  woCG <- sliding(wo[[paste0(stage,"_promoter_CG")]]*100)
  woCHG <- sliding(wo[[paste0(stage,"_promoter_CHG")]]*100)
  woCHH <- sliding(wo[[paste0(stage,"_promoter_CHH")]]*100)
  
  # window sizes
  TE_window <- calc_window(TEdf)
  g_window  <- calc_window(gdf)
  gexp_window  <- calc_window(gdfexp)
  pmC_window <- calc_window(pmC)
  wo_window <- calc_window(wo)


  #---- Plotting ----
  CG_col <- MODULE4_CG_COL
  CHG_col <- MODULE4_CHG_COL
  CHH_col <- MODULE4_CHH_COL
  TE_col <- "#B9653C"
  pro_wTE_col <- "#7F7F7F"
  pro_woTE_col <- "#BFBFBF"
  
# TE exp vs TE mC
te_ch_ylim <- expand_ylim(opt$ylim_TEexpTEmC_CH, c(TE_CHG, TE_CHH))
te_cg_ylim <- expand_ylim(opt$ylim_TEexpTEmC_CG, TE_CG)

# ----- pretty ticks -----
left_step <- 20
right_step <- 10

left_ticks <- seq(
  0,
  ceiling(te_cg_ylim / left_step) * left_step,
  by = left_step
)

right_ticks <- seq(
  0,
  ceiling(te_ch_ylim / right_step) * right_step,
  by = right_step
)

# map right ticks onto left axis coordinates
mapped_right_ticks <- right_ticks / max(right_ticks) * max(left_ticks)

png(
  file=output_file(paste0("OUTPUT_4_TEexp_TEmC_line_",stage,".png")),
  width=2500,
  height=2300,
  res=400
)

petem_base_par(list(
  mar=c(6,4.5,6,5)+0.1,
  bg="white",
  cex.axis=MODULE4_CEX_AXIS,
  cex.lab=MODULE4_CEX_LABEL,
  cex.main=MODULE4_CEX_LABEL
))

# ---- CH axis ----
plot(
  TE_CHG,
  lwd=PETEM_BASE_LINE_WIDTH,
  lty=1,
  col=CHG_col,
  type="n",
  axes=FALSE,
  ylim=c(0, max(right_ticks)),
  xlim=c(0,100),
  xlab=NA,
  ylab=NA,
  xaxt='n'
)

# aligned horizontal grid
abline(
  h=right_ticks,
  col=MODULE4_GRID_COL,
  lty=1,
  lwd=0.6
)

lines(
  TE_CHG,
  lwd=PETEM_BASE_LINE_WIDTH,
  lty=1,
  col=CHG_col
)

lines(
  TE_CHH,
  lwd=PETEM_BASE_LINE_WIDTH,
  lty=1,
  col=CHH_col
)

axis(
  4,
  at=right_ticks,
  labels=right_ticks,
  las=1,
  lwd=MODULE4_FRAME_WIDTH,
  font=2
)

mtext(
  expression(bold("TE CH methylation level (%)")),
  side=4,
  line=3.5,
  cex=MODULE4_CEX_LABEL
)

# ---- CG axis ----
par(new=TRUE)

plot(
  TE_CG,
  lwd=PETEM_BASE_LINE_WIDTH,
  lty=1,
  col=CG_col,
  type="l",
  axes=FALSE,
  ylim=c(0, max(left_ticks)),
  xlim=c(0,100),
  xlab=NA,
  ylab=NA,
  xaxt='n'
)

axis(
  2,
  at=left_ticks,
  labels=left_ticks,
  col="black",
  las=1,
  lwd=MODULE4_FRAME_WIDTH,
  font=2
)

mtext(
  expression(bold("TE CG methylation level (%)")),
  side=2,
  line=3,
  cex=MODULE4_CEX_LABEL
)

mtext(
  "Lowly expressed TEs      Highly expressed TEs",
  side=1,
  line=1,
  cex=MODULE4_CEX_LABEL,
  font=2
)

par(xpd=NA)

legend(
  "top",
  inset=c(0,-0.16),
  x.intersp=0.8,
  horiz=TRUE,
  legend=c("CG","CHG","CHH"),
  text.font=2,
  bty='n',
  lty=1,
  lwd=PETEM_BASE_LINE_WIDTH,
  col=c(CG_col,CHG_col,CHH_col),
  cex=MODULE4_CEX_AXIS
)

par(xpd=FALSE)

title(
  main = "Correlation between TE methylation \nand TE expression\n",
  cex.main = MODULE4_CEX_LABEL ,
  font.main = 2
)
box(lwd=MODULE4_FRAME_WIDTH)

mtext(
  paste0(
    "TEs: ",
    nrow(TEdf),
    ", window size: ",
    TE_window
  ),
  side=1,
  line=4,
  cex=MODULE4_CEX_NOTE+0.2
)

dev.off()
  
  # gene exp vs TE exp
  png(file=output_file(paste0("OUTPUT_4_geneexp_TEexp_line_",stage,".png")), width=2600,height=2200,res=400)
  petem_base_par(list(mar=c(5,4.5,4,5)+0.1, bg="white", cex.axis=MODULE4_CEX_AXIS+1, cex.lab=MODULE4_CEX_LABEL, cex.main=MODULE4_CEX_LABEL))
  plot(gTE_exp,lwd=PETEM_BASE_LINE_WIDTH,lty=1,col="gray50",type="n",axes=FALSE,xlim=c(0,100),xlab=NA,ylab=NA,xaxt='n')
  grid(nx=NA, ny=NULL, col=MODULE4_GRID_COL, lty=1, lwd=0.6)
  lines(gTE_exp,lwd=PETEM_BASE_LINE_WIDTH,lty=1,col="gray50")
  axis(2, ylim=c(0,1), col="black", las=1, lwd=MODULE4_FRAME_WIDTH, cex.axis=MODULE4_CEX_AXIS)
  mtext(expression(bold("TE expression (log2 RPKM)")),side=2,line=3,cex=MODULE4_CEX_LABEL)
  mtext("Lowly expressed genes      Highly expressed genes", side=1, line=1, cex=MODULE4_CEX_LABEL,font=2)
  box(lwd=MODULE4_FRAME_WIDTH)
  mtext(paste0("TE-gene pairs: ", nrow(gdfexp), ", window size: ", g_window), side=1, line=4, cex=MODULE4_CEX_NOTE+0.2)
  title(main = paste0("Correlation between TE expression\n", "and gene expression"),cex.main = MODULE4_CEX_LABEL+0.2 ,font.main = 2)

  dev.off()

  # gene exp vs TE/promoter mC

for(mtype in c("CG","CHG","CHH")){

    requested_ylim <- switch(
        mtype,
        CG = opt$ylim_CG,
        CHG = opt$ylim_CHG,
        CHH = opt$ylim_CHH
    )

    wo_vec <- get(paste0("wo", mtype))
    p_vec  <- get(paste0("p", mtype))
    te_vec <- get(paste0("gTE_", mtype))

    df_long <- bind_rows(
        data.frame(
            x = seq_along(wo_vec),
            type = "woTE",
            value = wo_vec
        ),
        data.frame(
            x = seq_along(p_vec),
            type = "wTE",
            value = p_vec
        ),
        data.frame(
            x = seq_along(te_vec),
            type = "TE",
            value = te_vec
        )
    )

    ylim_val <- expand_ylim(requested_ylim, df_long$value)

    line_colors <- c(
        "TE"   = CG_col,
        "wTE"  = pro_wTE_col,
        "woTE" = pro_woTE_col
    )

    line_labels <- c(
        "TE"   = "TE",
        "wTE"  = "Promoter containing TEs",
        "woTE" = "Promoter with no TEs"
    )

    

    rho_df <- data.frame(
        type = c("TE", "wTE", "woTE"),
        rho = c(
            cor(
                log2(gdf[[paste0(stage,"_exp")]]),
                gdf[[paste0(stage,"_TE_", mtype)]] * 100,
                method = "spearman",
                use = "complete.obs"
            ),
            cor(
                log2(pmC[[paste0(stage,"_exp")]]),
                pmC[[paste0(stage,"_promoterselves_", mtype)]] * 100,
                method = "spearman",
                use = "complete.obs"
            ),
            cor(
                log2(wo[[paste0(stage)]]),
                wo[[paste0(stage,"_promoter_", mtype)]] * 100,
                method = "spearman",
                use = "complete.obs"
            )
        ),
        stringsAsFactors = FALSE
    )

    rho_df$label <- sprintf("\u03c1=%.2f", rho_df$rho)

    line_end_df <- df_long %>%
        filter(is.finite(value)) %>%
        group_by(type) %>%
        filter(x == max(x)) %>%
        slice(1) %>%
        ungroup() %>%
        select(type, end_x = x, end_y = value)

    rho_df <- rho_df %>%
        left_join(line_end_df, by = "type") %>%
        filter(is.finite(rho), is.finite(end_x), is.finite(end_y))

    x_axis_max <- max(df_long$x, na.rm = TRUE) * 1.18
    if (nrow(rho_df) > 0) {
        rho_df$x <- rho_df$end_x + max(df_long$x, na.rm = TRUE) * 0.03
        rho_df$y <- rho_df$end_y + ylim_val * 0.04
        rho_df$y <- pmin(rho_df$y, ylim_val * 0.98)

        min_gap <- ylim_val * 0.07
        rho_df <- rho_df[order(rho_df$y), , drop = FALSE]

        if (nrow(rho_df) > 1) {
            for (i in 2:nrow(rho_df)) {
                if ((rho_df$y[i] - rho_df$y[i - 1]) < min_gap) {
                    rho_df$y[i] <- rho_df$y[i - 1] + min_gap
                }
            }
        }

        if (max(rho_df$y, na.rm = TRUE) > ylim_val * 0.98) {
            rho_df$y <- rho_df$y - (max(rho_df$y, na.rm = TRUE) - ylim_val * 0.98)
        }

        rho_df$y <- pmax(rho_df$y, ylim_val * 0.08)
    } else {
        rho_df$x <- numeric(0)
        rho_df$y <- numeric(0)
    }

    p <- ggplot(df_long, aes(x = x, y = value, color = type)) +
        geom_line(linewidth = 1.5) +
        geom_text(
            data = rho_df,
            aes(x = x, y = y, label = label, color = type),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 0.5,
            fontface = "bold",
            size = MODULE4_LINE_LABEL_SIZE + 1,
            show.legend = FALSE
        ) +
        scale_color_manual(values = line_colors, labels = line_labels, breaks = c( "TE", "wTE","woTE")) +
        coord_cartesian(ylim = c(0, ylim_val), xlim = c(0, x_axis_max)) +
        labs(    title = "Correlation between TE methylation\n and gene expression",
            y = paste0(mtype, " methylation level (%)"),
            x = "Lowly expressed genes             Highly expressed genes",
            caption = paste0(
                "TE-gene pairs: ", nrow(gdf),
                ", window size: ", g_window, "\n",
                "Promoter containing TEs: ", nrow(pmC),
                ", window size: ", pmC_window, "\n",
                "Promoter with no TEs: ", nrow(wo),
                ", window size: ", wo_window
            ),
            color = NULL
        ) +
        module4_theme_bw() +
        theme(
            legend.position = "top",
            legend.background = element_rect(fill = "white", color = NA),
            legend.box.background = element_rect(fill = "white", color = NA),
            legend.key = element_rect(fill = "white", color = NA),
            plot.margin = margin(t = 10, r = 10, b = 24, l = 10),
            plot.caption = element_text(
                hjust = 0.5,
                size = PETEM_AXIS_TEXT_SIZE + 1.5,
                lineheight = 1.1
            ),
            plot.title = element_text(size = PETEM_AXIS_TEXT_SIZE + 3,face = "bold",hjust = 0.5),
            plot.caption.position = "plot",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line = element_blank(),
            axis.text.y = element_text(size = PETEM_AXIS_TEXT_SIZE + 2)
        )

    ggsave(
        filename = output_file(paste0("OUTPUT_4_geneexp_TEm", mtype, "_line_", stage, ".png")),
        plot = p,
        width = 2300 / 400,
        height = 2500 / 400,
        dpi = 400,
        units = "in"
    )
}


  
  # Correlation tables
  # Pearson
  plot_corr_bar(stage, gdf, gdfexp, TEdf, method="pearson")
  # Spearman
  plot_corr_bar(stage, gdf, gdfexp, TEdf, method="spearman")
}

end_time <- Sys.time()
print(end_time-start_time)
