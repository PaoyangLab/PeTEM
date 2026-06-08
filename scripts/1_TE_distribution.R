# Rscript 1_TE_distribution.R
start_time <- Sys.time()

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
source(file.path(dirname(normalizePath(script_path)), "plot_defaults.R"))

args <- commandArgs(trailingOnly = TRUE)

# Bed files
ins_gene_file <- args[1]
ins_CDS_file <- args[2]
ins_5UTR_file <- args[3]
ins_exon_file <- args[4]
ins_intron_file <- args[5]
ins_3UTR_file <- args[6]
ins_promoter_file <- args[7]
ins_IGR_file <- args[8]
TE_file <- args[9]
gene_file <- args[10]
CDS_file <- args[11]
UTR5_file <- args[12]
exon_file <- args[13]
intron_file <- args[14]
UTR3_file <- args[15]
promoter_file <- args[16]
IGR_file <- args[17]
genome_file <- args[18]


#------------------
# Load data
ins_gene <- read.table(ins_gene_file)
ins_CDS <- read.table(ins_CDS_file)
ins_5UTR <- read.table(ins_5UTR_file)
ins_exon <- read.table(ins_exon_file)
ins_intron <- read.table(ins_intron_file)
ins_3UTR <- read.table(ins_3UTR_file)
ins_promoter <- read.table(ins_promoter_file)
ins_IGR <- read.table(ins_IGR_file)
TEnew <- read.table(TE_file)
TEnew$V7 <- TEnew$V3 - TEnew$V2

TE_ins <- function(df) {
  df2 <- df[!duplicated(df[, c("V1", "V2", "V3")]), ]
  return(sum(df2$V3 - df2$V2) )
}

TE_insertion <- c(TE_ins(ins_promoter), TE_ins(ins_gene), TE_ins(ins_5UTR), TE_ins(ins_CDS),
                  TE_ins(ins_exon), TE_ins(ins_intron), TE_ins(ins_3UTR), TE_ins(ins_IGR))
TE_insertion2 <- TE_insertion*100/sum(as.numeric(TEnew$V7))

# Load genomic regions
gene <- read.table(gene_file)
CDS <- read.table(CDS_file)
UTR5 <- read.table(UTR5_file)
exon <- read.table(exon_file)
intron <- read.table(intron_file)
UTR3 <- read.table(UTR3_file)
promoter <- read.table(promoter_file)
IGR <- read.table(IGR_file)
genome <- read.table(genome_file)

genome_percent <- function(df){
  df$V7 <- (df$V3 - df$V2) 
  sum(df$V7)*100/sum(genome$V3)
}

whole_genome <- c(genome_percent(promoter), genome_percent(gene), genome_percent(UTR5),
                  genome_percent(CDS), genome_percent(exon), genome_percent(intron),
                  genome_percent(UTR3), genome_percent(IGR))
TEinsert_enrich <- log2(TE_insertion2 / whole_genome)

regions <- c("Promoter","Gene","5'UTR","CDS","Exon","Intron","3'UTR","IGR")
colors <- c("Promoter"="#BD7634","Gene"="#000000","5'UTR"="#666666","CDS"="#333333",
            "Exon"="#999999","Intron"="#BBBBBB","3'UTR"="#DDDDDD","IGR"="#D8C5AF")

library(ggplot2)
suppressPackageStartupMessages(library(ggbreak))

plot_bar <- function(values, ylab, title, yuplim, y_limit = NA){
  df <- data.frame(Region=factor(regions, levels=regions), Value=values)
  df$label_y <- ifelse(df$Value < 0, 0, df$Value)
  df$vjust <- ifelse(df$Value < 0, -0.5, -0.5)

  p <- ggplot(df, aes(x=Region, y=Value, fill=Region)) +
    geom_bar(stat="identity", color="black") +
    geom_text(aes(x=Region, y=label_y, label=sprintf("%.2f", Value), vjust=vjust),
              size=PETEM_ANNOTATION_TEXT_SIZE + 1.8, color="firebrick4", fontface="bold") +
    scale_fill_manual(values=colors) +
    petem_theme_classic() +
    theme(axis.text.x = element_text(angle=40, hjust=1, size = PETEM_AXIS_TEXT_SIZE + 3.5),
          axis.text.y = element_text(size = PETEM_AXIS_TEXT_SIZE + 3.5),
          axis.title.x = element_text(size = PETEM_AXIS_TITLE_SIZE + 3.5),
          axis.title.y = element_text(size = PETEM_AXIS_TITLE_SIZE + 3.5),
          plot.title = element_text(size = PETEM_PLOT_TITLE_SIZE + 3.5),
          axis.ticks.x = element_blank(),
          legend.position="none") +
    labs(x="", y=ylab, title=title)
  if (is.na(y_limit)) {
    p <- p + expand_limits(y = max(df$Value) * yuplim)
  } else {
    low <- min(0, min(df$Value))
    p <- p + scale_y_continuous(limits = c(low, y_limit))
  }
  return(p)
}

#------------------
png("OUTPUT_1_TE_distribution_enrichment.png", width=2000, height=1800, res=300)
print(plot_bar(TEinsert_enrich, expression(bold("Enrichment (Log2)")), "Enrichment of TE distribution", 1.8))
dev.off()

png("OUTPUT_1_TE_distribution_percentage.png", width=2000, height=1800, res=300)
print(plot_bar(TE_insertion2, "Percentage of TEs (%)", "Distribution of TEs across genomic features", 1.2, y_limit = 100))
dev.off()

end_time <- Sys.time()
print(end_time-start_time)
