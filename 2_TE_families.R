#!/usr/bin/env Rscript

# Rscript 2_TE_families.R -a TE.bed -i TE_overlap_promoter.bed -T TE_family.txt 

start_time <- Sys.time()

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
source(file.path(dirname(normalizePath(script_path)), "plot_defaults.R"))

# ==== Load packages ====
library(optparse)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)


# Options -------------
option_list <- list(
  make_option(c("-a", "--all"),   type="character", help="All TE annotation file (bed format)"),
  make_option(c("-i", "--ins"),   type="character", help="Annotation file of TEs overlapping with promoters (bed format)"),
  make_option(c("-T", "--TE"),    type="character", help="TE families file (bed format)"),
  make_option(c("-o", "--outdir"), type="character", default=".", help="Output directory")
)

args <- parse_args(OptionParser(option_list=option_list))

if (is.null(args$all) | is.null(args$ins) | is.null(args$TE)) {
  cat("Missing required arguments\n\n")
  print_help(OptionParser(option_list=option_list))
  quit(status=1)
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

read_te_family_table <- function(path) {
  raw_df <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
  if (ncol(raw_df) == 2) {
    out <- raw_df[, 1:2]
  } else if (ncol(raw_df) >= 7) {
    out <- raw_df[, c(1, 7)]
  } else {
    stop("Unsupported TE family format. Expected either 2 columns (TE_id, family) or >=7 columns with TE family in column 7.")
  }
  colnames(out) <- c("TE_id", "TE_family")
  out
}


# Read input files ----------
df1 <- read.table(args$all, sep="\t", header=FALSE)
df2 <- read.table(args$ins, sep="\t", header=FALSE)
TE_families <- read_te_family_table(args$TE)


df1<-df1[,c(1:6)]
df1<-df1[!duplicated(df1),]

df2<-df2[,c(1:6)]
df2<-df2[!duplicated(df2),]

# Count TE families ----------
df1m <- merge(df1, TE_families, by.x="V4", by.y="TE_id")
df2m <- merge(df2, TE_families, by.x="V4", by.y="TE_id")

c1 <- as.data.frame(table(df1m$TE_family))
c2 <- as.data.frame(table(df2m$TE_family))

df_TE <- merge(c1, c2, by="Var1", all=TRUE)
df_TE[is.na(df_TE)] <- 0
colnames(df_TE)<-c("TE", "df1", "df2")

print(df_TE)

# totals
total_df1 <- sum(df_TE$df1)
total_df2 <- sum(df_TE$df2)

# percentage
df_TE$pc_df1 <- df_TE$df1 * 100 / total_df1
df_TE$pc_df2 <- df_TE$df2 * 100 / total_df2

# enrichment
df_TE$enrich <- log2((df_TE$df2 / total_df2) / (df_TE$df1 / total_df1))

# Fisher exact test
p_list <- numeric(nrow(df_TE))
for(i in seq_len(nrow(df_TE))){
  a <- df_TE$df2[i]
  b <- df_TE$df1[i]
  c <- total_df2 - a
  d <- total_df1 - b
  contingency <- matrix(c(a,c,b,d), nrow=2, byrow=TRUE)
  p_list[i] <- fisher.test(contingency)$p.value
}

df_TE$pvalue_num <- p_list
out_df_TE<-df_TE

df_TE$pvalue <- format.pval(p_list, digits=3, scientific=TRUE)

df_TE <- df_TE[order(df_TE$enrich, decreasing = TRUE), , drop = FALSE]
df_TE$TE <- factor(df_TE$TE, levels = df_TE$TE)

colnames(out_df_TE)<-c("TE family","All TEs", "Promoter-embedded TEs", "All TEs (%)", "Promoter-embedded TEs (%)", "Log2 enrichment", "pvalue")
write.table(out_df_TE, file=file.path(args$outdir, "OUTPUT_2_Promoter_embedded_TE_family.txt"), sep="\t", quote=F, row.names=F)

# labels
df_TE$signif <- ifelse(df_TE$pvalue_num < 0.001, "***",
                ifelse(df_TE$pvalue_num < 0.01, "**",
                ifelse(df_TE$pvalue_num < 0.05, "*", "ns")))
df_TE$text <- paste0("Log2 FC=", sprintf("%.2f", df_TE$enrich), ", ", df_TE$signif)
df_TE$text_y <- 100 - (cumsum(df_TE$pc_df2) - df_TE$pc_df2 / 2)

df_long <- data.frame(
  TE = rep(df_TE$TE, times = 2),
  type = rep(c("All TEs", "Promoter-\nembedded TEs"), each = nrow(df_TE)),
  percentage = c(df_TE$pc_df1, df_TE$pc_df2),
  stringsAsFactors = FALSE
)

set3_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(df_long$TE)))
names(set3_colors) <- levels(df_TE$TE)

df_pos <- df_TE[df_TE$enrich > 0, ]
df_pos <- df_pos[order(df_pos$pvalue_num, -df_pos$enrich), , drop = FALSE][1:3, , drop = FALSE]
df_neg <- df_TE[df_TE$enrich < 0, ]
df_neg <- df_neg[order(df_neg$pvalue_num, df_neg$enrich), , drop = FALSE][1:3, , drop = FALSE]
df_label <- rbind(df_pos, df_neg)

df_label <- df_label[order(df_label$text_y), , drop = FALSE]
df_label$y_start <- df_label$text_y
df_label$text_y <- seq(25, 75, length.out = nrow(df_label))
df_label$y_end <- df_label$text_y

# Plot
png(file=file.path(args$outdir, "OUTPUT_2_Promoter_embedded_TE_family_enrichment.png"), width=4500, height=2200, res=400)

ggplot(df_long, aes(x = type, y = percentage, alluvium = TE)) +
  geom_alluvium(aes(fill = TE), width = 0.3, alpha = 0.6) +
  geom_stratum(aes(stratum = TE, fill = TE), width = 0.3) +
  geom_segment(data = df_label,
               aes(x = 2.12, xend = 2.5, y = y_start, yend = y_end),
               color = "gray30", linewidth = 0.7, inherit.aes = FALSE) +
  scale_x_discrete(limits = c("All TEs", "Promoter-\nembedded TEs", "", "", "")) +
  scale_y_continuous(limits = c(0,100), expand=c(0,0)) +
  scale_fill_manual(values = set3_colors) +
  petem_theme_classic() +
  labs(title = "Enriched families of promoter-embedded TEs", x="", y="Percentage (%)", fill="TE Types") +
  theme(
    axis.text.x = element_text(size = PETEM_AXIS_TEXT_SIZE + 2),
    axis.text.y = element_text(size = PETEM_AXIS_TEXT_SIZE + 2),
    axis.title.y = element_text(size = PETEM_AXIS_TITLE_SIZE + 2),
    plot.title = element_text(size = PETEM_PLOT_TITLE_SIZE + 2),
    legend.title = element_text(size = PETEM_LEGEND_TITLE_SIZE + 2),
    legend.text = element_text(size = PETEM_LEGEND_TEXT_SIZE + 2),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color='black'),
    axis.line.y.left = element_line(color='black'),
    legend.key.size = unit(1.4, "lines")
  ) +
  geom_text(data = df_label,
            aes(x = 2.6, y = y_end, label = text),
            hjust = 0, size = PETEM_ANNOTATION_TEXT_SIZE + 1, inherit.aes = FALSE)

dev.off()

end_time <- Sys.time()
print(end_time-start_time)
