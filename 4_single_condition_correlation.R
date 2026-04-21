#!/usr/bin/env Rscript
# Rscript 4_single_condition_correlation.R --eg expression_gene.txt --et expression_TE.txt --unexp y/n \
#   --wd_num 156 --ylim_CG 40 --ylim_CHG 5 --ylim_CHH 5 --ylim_TEexpTEmC_CH 8 --ylim_TEexpTEmC_CG 40


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

#---- Option parser ----
option_list = list(
  make_option(c("--eg"), type="character", help="Gene expression file"),
  make_option(c("--et"), type="character", help="TE expression file"),
  make_option(c("--module0-dir"), type="character", default=".", help="Directory containing Module 0 outputs such as Tab_* and TE_overlap_promoter.bed"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--unexp"), type="character", default="n", help="Include unexpressed TEs for sliding plots? (y/n)"),
  make_option(c("--wd_num"), type="numeric", default=156, help="number of sliding window (default=156)"),
  make_option(c("--ylim_CG"), type="numeric", default=50, help="ylim for gene exp vs TE/promoter mC, CG (default=50)"),
  make_option(c("--ylim_CHG"), type="numeric", default=10, help="ylim for gene exp vs TE/promoter mC, CHG (default=10)"),
  make_option(c("--ylim_CHH"), type="numeric", default=10, help="ylim for gene exp vs TE/promoter mC, CHH (default=10)"),
  make_option(c("--ylim_TEexpTEmC_CH"), type="numeric", default=15, help="ylim for TE exp vs TE mC, CH (default=15)"),
  make_option(c("--ylim_TEexpTEmC_CG"), type="numeric", default=30, help="ylim for TE exp vs TE mC, CG (default=30)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

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
  step <- floor(lth/opt$wd_num)
  window <- lth-(opt$wd_num-1)*step
  return(as.vector(rollapply(zoo(vec), width=window, by=step, FUN=mean, align="center")))
}

calc_window <- function(df){
  step <- floor(nrow(df)/opt$wd_num)
  window <- nrow(df)-(opt$wd_num-1)*step
  return(window)
}

corr_mC <- function(df, exp, CG, CHG, CHH, method="pearson") {
  log2exp<-as.numeric(log2(df[[exp]]))
  data.frame(
    Corr = c(
      cor.test(log2exp, df[[CG]], method=method)$estimate,
      cor.test(log2exp, df[[CHG]], method=method)$estimate,
      cor.test(log2exp, df[[CHH]], method=method)$estimate
    ),
    Pvalue = c(
      cor.test(log2exp, df[[CG]], method=method)$p.value,
      cor.test(log2exp, df[[CHG]], method=method)$p.value,
      cor.test(log2exp, df[[CHH]], method=method)$p.value
    ),
    row.names = c("CG", "CHG", "CHH")
  )
}

corr_exp <- function(df, x, y, method="pearson") {
  corr <- cor.test(as.numeric(log2(df[[x]])), as.numeric(log2(df[[y]])), method=method)
  data.frame(Corr=corr$estimate, Pvalue=corr$p.value, row.names=paste(x, y, sep="_vs_"))
}

plot_corr_bar <- function(stage, gdf, gdfexp, TEdf, method="pearson") {
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
    Color = c(rep("#D3A355",3), "#2B5A78", rep("#A4BE78",3))
  )
  corr_df$Comparison <- factor(corr_df$Comparison, levels=rev(corr_df$Comparison))
  
  # plot
  xmin <- min(corr_df$Corr, na.rm=TRUE) - 0.1
  xmax <- max(corr_df$Corr, na.rm=TRUE) + 0.1
  outfile <- output_file(paste0("OUTPUT_4_",method,"_correlation_bar_",stage,".png"))
  png(file=outfile,width=2000,height=1500,res=400)
  print(
    ggplot(corr_df, aes(x=Corr,y=Comparison,fill=Color)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=round(Corr,2)),
                hjust=ifelse(corr_df$Corr>=0, -0.1, 1.1), size=4) +
      scale_fill_identity() +
      xlim(xmin, xmax) +
      labs(x=paste0(method," correlation coefficient"), y=NULL) +
      petem_theme_bw() +
      theme(
        panel.border=element_rect(colour="black", fill=NA, linewidth=1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
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
if (length(colnames(gene_exp)) != 2 || length(colnames(TE_exp)) != 2) {
  stop("Module 4 requires exactly two expression stages/conditions in both gene and TE expression matrices.")
}
if (length(stages) != 2) {
  stop("Module 4 requires exactly two shared stages/conditions between gene and TE expression matrices.")
}

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
  pmC <- sort_exp(insGenePromC[,c("V4","V10",paste0(stage,"_exp"),paste0(stage,"_TEexp"),
                                  paste0(stage,"_promoterselves_CG"),paste0(stage,"_promoterselves_CHG"),paste0(stage,"_promoterselves_CHH"))],
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
    pmC  <- pmC[(pmC[[paste0(stage,"_exp")]] != 0), ]
  } else {
    gdf  <- gdf[(gdf[[paste0(stage,"_exp")]] != 0) & (gdf[[paste0(stage,"_TEexp")]] != 0), ]
    pmC  <- pmC[(pmC[[paste0(stage,"_exp")]] != 0) & (pmC[[paste0(stage,"_TEexp")]] != 0), ]
  }

  pmC  <- pmC[, c(2,3,5,6,7)]
  pmC <- unique(pmC)

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
  #pmC_window <- calc_window(pmC)
  wo_window <- calc_window(wo)


  #---- Plotting ----
  CG_col <- "#CFA699"
  CHG_col <- "#71A061"
  CHH_col <- "#3871A6"
  TE_col <- "#B9653C"
  pro_wTE_col <- "#7F7F7F"
  pro_woTE_col <- "#BFBFBF"
  
  # TE exp vs TE mC
  png(file=output_file(paste0("OUTPUT_4_TEexp_TEmC_line_",stage,".png")), width=2600,height=2200,res=400)
  petem_base_par(list(mar=c(6,4.5,6,5)+0.1))
  plot(TE_CHG,lwd=PETEM_BASE_LINE_WIDTH,lty=1,col=CHG_col,type="l",axes=FALSE,
       ylim=c(0,opt$ylim_TEexpTEmC_CH),xlim=c(0,opt$wd_num),xlab=NA,ylab=NA,xaxt='n')
  lines(TE_CHH,lwd=PETEM_BASE_LINE_WIDTH,lty=1,col=CHH_col)
  axis(4,las=1)
  mtext(expression(bold("TE mCH(%)")), side=4, line=3.5, cex=PETEM_BASE_CEX_LABEL)
  box()
  par(new=TRUE)
  plot(TE_CG,lwd=PETEM_BASE_LINE_WIDTH,lty=1,col=CG_col,type="l",axes=FALSE,ylim=c(0,opt$ylim_TEexpTEmC_CG),xlim=c(0,opt$wd_num),xlab=NA,ylab=NA,xaxt='n')
  axis(2, col="black", las=1)
  mtext(expression(bold("TE mCG(%)")), side=2, line=3, cex=PETEM_BASE_CEX_LABEL)
  mtext("Lowly expressed TEs      Highly expressed TEs", side=1, line=1, cex=PETEM_BASE_CEX_LABEL, font=2)
  par(xpd=NA)
  legend("top", inset=c(0,-0.08), x.intersp=0.8, horiz=FALSE,
         legend=c("CG","CHG","CHH"), text.font=2, bty='n', lty=1, lwd=PETEM_BASE_LINE_WIDTH, col=c(CG_col,CHG_col,CHH_col), cex=PETEM_BASE_CEX_AXIS)
  par(xpd=FALSE)
  grid(nx=NA, ny=NULL, col="gray70", lty=3, lwd=1)
  mtext(paste0("TEs: ", nrow(TEdf), ", window size: ", TE_window), side=1, line=4, cex=PETEM_BASE_CEX_NOTE)
  dev.off()
  
  # gene exp vs TE exp
  png(file=output_file(paste0("OUTPUT_4_geneexp_TEexp_line_",stage,".png")), width=2600,height=2200,res=400)
  petem_base_par(list(mar=c(5,4.5,4,5)+0.1))
  plot(gTE_exp,lwd=PETEM_BASE_LINE_WIDTH,lty=1,col="gray50",type="l",axes=FALSE,xlim=c(0,opt$wd_num),xlab=NA,ylab=NA,xaxt='n')
  axis(2, ylim=c(0,1),col="black",las=1)
  mtext(expression(bold("TE expression (log2 RPKM)")),side=2,line=3,cex=PETEM_BASE_CEX_LABEL)
  mtext("Lowly expressed genes      Highly expressed genes", side=1, line=1, cex=PETEM_BASE_CEX_LABEL,font=2)
  grid(nx=NA,ny=NULL,col="gray70",lty=3,lwd=1)
  box()
  mtext(paste0("TE-gene pairs: ", nrow(gdfexp), ", window size: ", gexp_window), side=1, line=4, cex=PETEM_BASE_CEX_NOTE)
  dev.off()
  
  # gene exp vs TE/promoter mC

for(mtype in c("CG","CHG","CHH")){

    # Select y-axis upper limit.
    ylim_val <- switch(mtype, CG=opt$ylim_CG, CHG=opt$ylim_CHG, CHH=opt$ylim_CHH)

    df <- data.frame(
        x = seq_along(get(paste0("wo", mtype))),
        woTE = get(paste0("wo", mtype)),
        wTE  = get(paste0("p", mtype)),
        TE   = get(paste0("gTE_", mtype))
    )

    df_long <- df %>%
        pivot_longer(cols=c("woTE","wTE","TE"),
                     names_to="type", values_to="value")

    line_colors <- c(
        "TE"   = CG_col,
        "wTE"  = pro_wTE_col,
        "woTE" = pro_woTE_col
    )

    line_labels <- c(
        "TE"   = "TE",
        "wTE"  = "Promoters w TEs",
        "woTE" = "Promoters w/o TEs"
    )

    p <- ggplot(df_long, aes(x=x, y=value, color=type)) +
        geom_line(linewidth=1.5) +
        scale_color_manual(values=line_colors, labels=line_labels) +
        coord_cartesian(ylim=c(0, ylim_val), xlim=c(0, opt$wd_num)) +
        labs(y="Methylation level (%)",
             x="Lowly expressed genes             Highly expressed genes",
             color=NULL) +
        petem_theme_bw() +
        theme(
            legend.position = "top",
            panel.grid.minor = element_blank()
        ) +
        annotate("text", 
                 x = Inf, y = -Inf, hjust=1.05, vjust=-1.5, size=PETEM_ANNOTATION_TEXT_SIZE,
                 label = paste0(
                     "TE-gene pairs: ", nrow(gdf),
                     ", window size: ", g_window,
                     "; Promoters w/o TEs: ", nrow(wo),
                     ", window size: ", wo_window
                 ))

    ggsave(
        filename = output_file(paste0("OUTPUT_4_geneexp_TEm",mtype,"_line_",stage,".png")),
        plot = p,
        width = 2600/400, height = 2200/400, dpi=400, units="in"
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
