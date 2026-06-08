#!/usr/bin/env bash
set -euo pipefail

MODULE_NAME="Module 3"
SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
CURRENT_STEP="argument parsing"
LOG="LOG_3_2_TE_impact_distance_plot.log"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

die() {
    echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: $*" >&2
    exit 1
}

usage_missing() {
    echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: missing required argument(s): $*" >&2
    usage
}

require_file() {
    local flag=$1
    local path=$2
    [[ -f $path ]] || die "input file for ${flag} not found: ${path}"
}

trap 'rc=$?; echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: command failed with exit code ${rc}" >&2; exit ${rc}' ERR

# ====================
# usage
# ====================
usage() {
    echo "Usage: bash 3_2_TE_impact_distance_plot.sh -g gene.bed -t TE.bed -eg expression_gene.txt -et expression_TE.txt -lim 15000 -tick 5000 -WD 200 -unexp y/n [-types CG,CHG,CHH,raw] [-font 20]" >&2
    exit 1
}

echo "[$(date)] Pipeline started" | tee "$LOG"
start_allall=$(date +%s)

# ====================
# parse args
# ====================
LIMIT=15000
MAJOR_TICK=5000
WD=200
INCLUDE_UNEXPRESSED_TE="n"
TYPES="CG,CHG,CHH"
BASE_FONT_SIZE=20

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -g) GENE_BED="$2"; shift 2 ;;
    -t) TE_BED="$2"; shift 2 ;;
    -eg) GENE_EXP="$2"; shift 2 ;;
    -et) TE_EXP="$2"; shift 2 ;;
    -lim) LIMIT="$2"; shift 2 ;;
    -tick) MAJOR_TICK="$2"; shift 2 ;;
    -WD) WD="$2"; shift 2 ;;
    -unexp) INCLUDE_UNEXPRESSED_TE="$2"; shift 2 ;;
    -types) TYPES="$2"; shift 2 ;;
    -font) BASE_FONT_SIZE="$2"; shift 2 ;;
    *) die "unknown option: ${key}" ;;
  esac
done

missing_args=()
[[ -n ${GENE_BED:-} ]] || missing_args+=("-g gene.bed")
[[ -n ${TE_BED:-} ]] || missing_args+=("-t TE.bed")
[[ -n ${GENE_EXP:-} ]] || missing_args+=("-eg expression_gene.txt")
[[ -n ${TE_EXP:-} ]] || missing_args+=("-et expression_TE.txt")
(( ${#missing_args[@]} == 0 )) || usage_missing "${missing_args[*]}"

CURRENT_STEP="input validation"
require_file "-g" "$GENE_BED"
require_file "-t" "$TE_BED"
require_file "-eg" "$GENE_EXP"
require_file "-et" "$TE_EXP"
[[ -d pre_step3 ]] || die "required directory not found: pre_step3. Run module 3 preprocessing first."

CURRENT_STEP="step 1 - prepare sorted BED files"
sort -k1,1 -k2,2n "$GENE_BED" > gene_sort.bed
sort -k1,1 -k2,2n "$TE_BED" > TE_sort.bed

echo "[$(date)] Input files parsed" | tee -a "$LOG"

# ====================
# step 1: find nearby TEs for each gene
# ====================
CURRENT_STEP="step 1 - find nearby TEs"
echo "[$(date)] Step 1. Preprocessing: find nearby TEs for each gene" | tee -a "$LOG"

Rscript - "$GENE_EXP" "$TE_EXP" "$INCLUDE_UNEXPRESSED_TE" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
gene_exp <- args[1]
TE_exp <- args[2]
include_unexp <- args[3]

gene_bed <- read.table("gene_sort.bed", header=FALSE, stringsAsFactors=FALSE)
TE_bed <- read.table("TE_sort.bed", header=FALSE, stringsAsFactors=FALSE)
gene_exp <- read.table(gene_exp, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
TE_exp <- read.table(TE_exp, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

gene_exp <- gene_exp[rowSums(gene_exp, na.rm=TRUE) != 0, , drop=FALSE]
if(tolower(include_unexp) != "y"){
  TE_exp <- TE_exp[rowSums(TE_exp, na.rm=TRUE) != 0, , drop=FALSE]
}

gene_bed2 <- gene_bed[gene_bed$V4 %in% row.names(gene_exp), ]
TE_bed2 <- TE_bed[TE_bed$V4 %in% row.names(TE_exp), ]

write.table(gene_bed2, "expressed_gene.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(TE_bed2, "expressed_TE.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
EOF

bedtools closest -a expressed_gene.bed -b expressed_TE.bed -id -d -D a > expgene_closest_expTE.bed

# ====================
# step 2: split high/low genes and generate upstream/downstream BED files
# ====================
CURRENT_STEP="step 2 - build expression strata and adjacent regions"
echo "[$(date)] Step 2. Split highly and lowly expressed genes and build neighboring TE regions" | tee -a "$LOG"

Rscript - "$GENE_EXP" "$LIMIT" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
gene_exp_file <- args[1]
limit <- as.numeric(args[2])

gene_bed <- read.table("gene_sort.bed", header=FALSE, stringsAsFactors=FALSE)
gene_exp <- read.table(gene_exp_file, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

adjacent <- function(df, up=TRUE, limit=limit){
  df2 <- df
  if(up){
    df2[df2$V6 == "+", 3] <- df2[df2$V6 == "+", 2] - 1
    df2[df2$V6 == "+", 2] <- df2[df2$V6 == "+", 3] - limit
    df2[df2$V6 == "-", 2] <- df2[df2$V6 == "-", 3] + 1
    df2[df2$V6 == "-", 3] <- df2[df2$V6 == "-", 2] + limit
  } else {
    df2[df2$V6 == "-", 3] <- df2[df2$V6 == "-", 2] - 1
    df2[df2$V6 == "-", 2] <- df2[df2$V6 == "-", 3] - limit
    df2[df2$V6 == "+", 2] <- df2[df2$V6 == "+", 3] + 1
    df2[df2$V6 == "+", 3] <- df2[df2$V6 == "+", 2] + limit
  }
  df2$V2[df2$V2 < 1] <- 1
  df2[5] <- 0
  return(df2)
}

clo_TE <- read.table("expgene_closest_expTE.bed", sep="\t", header=FALSE)
clo_TE2 <- clo_TE[clo_TE$V10 != ".", ]
clo_TE2 <- clo_TE2[clo_TE2$V13 >= (-limit), ]

gene_exp2 <- gene_exp[row.names(gene_exp) %in% clo_TE2$V4, , drop=FALSE]
stages <- colnames(gene_exp2)
sink("OUTPUT_3_TE_impact_distance_gene_TE_number.txt")

for(stage in stages){
  vals <- gene_exp2[, stage, drop=FALSE]
  vals <- vals[vals[, 1] > 0, , drop=FALSE]

  sorted <- vals[order(vals[, 1]), , drop=FALSE]
  qn <- floor(nrow(sorted) / 4)
  if(qn < 1){
    cat(stage, ": skipped; fewer than 4 expressed genes with nearby TEs\n\n")
    next
  }

  low <- head(sorted, qn)
  high <- tail(sorted, qn)

  low_bed <- merge(gene_bed, low, by.x="V4", by.y="row.names", all.y=TRUE)
  high_bed <- merge(gene_bed, high, by.x="V4", by.y="row.names", all.y=TRUE)
  low_bed <- low_bed[, c(2,3,4,1,7,6)]
  high_bed <- high_bed[, c(2,3,4,1,7,6)]

  write.table(high_bed, paste0("high_", stage, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(low_bed, paste0("low_", stage, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

  low_up <- adjacent(low_bed, up=TRUE, limit=limit)
  high_up <- adjacent(high_bed, up=TRUE, limit=limit)
  low_down <- adjacent(low_bed, up=FALSE, limit=limit)
  high_down <- adjacent(high_bed, up=FALSE, limit=limit)

  write.table(low_up, paste0("low_", stage, "_up.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(high_up, paste0("high_", stage, "_up.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(low_down, paste0("low_", stage, "_down.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(high_down, paste0("high_", stage, "_down.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

  cat(stage, ":", "\n")
  cat("  Highly expressed gene number:", nrow(high_bed), "\n")
  cat("  Lowly expressed gene number:", nrow(low_bed), "\n")

  high_genes <- high_bed$V4
  low_genes <- low_bed$V4

  high_te <- length(unique(clo_TE2$V10[clo_TE2$V4 %in% high_genes]))
  low_te <- length(unique(clo_TE2$V10[clo_TE2$V4 %in% low_genes]))

  cat("  TEs nearby highly expressed genes:", high_te, "\n")
  cat("  TEs nearby lowly expressed genes:", low_te, "\n\n")
}
sink()
EOF

# ====================
# step 3: intersect TE regions with methylation tables
# ====================
start_step4=$(date +%s)
CURRENT_STEP="step 3 - intersect methylation tables"
echo "[$(date)] Step 3. Intersect TE methylation data" | tee -a "$LOG"

IFS=',' read -r -a METH_TYPES <<< "$TYPES"
stages=($(ls *_up.bed 2>/dev/null | cut -d'_' -f2 | sort -u))

for stage in "${stages[@]}"; do
(
    start=$(date +%s)
    echo "[INFO] Processing stage $stage" | tee -a "$LOG"

    missing_types=()
    for type in "${METH_TYPES[@]}"; do
        if [[ ! -f "pre_step3/pre3_${stage}_${type}.bed" ]]; then
            missing_types+=("$type")
        fi
    done
    if [[ ${#missing_types[@]} -gt 0 ]]; then
        echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: missing methylation files for stage ${stage}: ${missing_types[*]}" | tee -a "$LOG"
        exit 1
    fi

    for expr in low high; do
        for dir in up down; do
            input="${expr}_${stage}_${dir}.bed"
            wTE="wTE_${expr}_${stage}_${dir}.bed"
            bedtools intersect -a "$input" -b expressed_TE.bed > "$wTE"

            for type in "${METH_TYPES[@]}"; do
                bedtools intersect -a "$wTE" -b "pre_step3/pre3_${stage}_${type}.bed" -wa -wb > "wTE_${expr}_${stage}_${dir}_${type}.bed"
            done

            rm -f "$wTE"
        done
    done

    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec" | tee -a "$LOG"
)&
done
wait

end_step4=$(date +%s)
echo "[$(date)] Step 3 finished in $((end_step4-start_step4)) sec" | tee -a "$LOG"

# ====================
# step 4: plot module 3 outputs
# ====================
start_step5=$(date +%s)
CURRENT_STEP="step 4 - generate plots"
echo "[$(date)] Step 4. Calculate regional methylation and generate plots" | tee -a "$LOG"

shopt -s nullglob
stage_files=(wTE_low_*_up_*.bed)
shopt -u nullglob

stages=()
declare -A seen_stage=()
for file in "${stage_files[@]}"; do
    stage=$(basename "$file" | cut -d'_' -f3)
    [[ -z "$stage" ]] && continue
    if [[ -z ${seen_stage[$stage]:-} ]]; then
        stages+=("$stage")
        seen_stage[$stage]=1
    fi
done

if [[ ${#stages[@]} -eq 0 ]]; then
    echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: no TE methylation overlap files found for plotting." | tee -a "$LOG"
    exit 1
fi

for stage in "${stages[@]}"; do
{
    start=$(date +%s)
    echo "[INFO] Processing stage $stage" | tee -a "$LOG"

    Rscript - "$LIMIT" "$MAJOR_TICK" "$WD" "$stage" "$SCRIPT_DIR" "$TYPES" "$BASE_FONT_SIZE" <<'EOF'
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))

args <- commandArgs(trailingOnly=TRUE)
limit <- as.numeric(args[1])
major_tick <- as.numeric(args[2])
WD <- as.numeric(args[3])
stage <- args[4]
script_dir <- args[5]
types <- trimws(strsplit(args[6], ",")[[1]])
base_font_size <- as.numeric(args[7])
source(file.path(script_dir, "plot_defaults.R"))

safe_read <- function(path){
  if(!file.exists(path) || file.info(path)$size == 0){
    return(NULL)
  }
  x <- tryCatch(read.table(path, sep="\t", header=FALSE, stringsAsFactors=FALSE), error=function(e) NULL)
  if(is.null(x) || nrow(x) == 0){
    return(NULL)
  }
  return(x)
}

for_up <- function(mydf, mybed){
  if(is.null(mydf) || is.null(mybed)) return(data.frame())
  mydf2 <- mydf[, c("V4", "V8", "V11", "V12")]
  mybed2 <- mybed[, c("V2", "V3", "V4", "V6")]
  mymy <- merge(mydf2, mybed2, by.x=c("V4", "V12"), by.y=c("V4", "V6"))
  rev <- mymy[mymy$V12 == "-", ]
  fow <- mymy[mymy$V12 == "+", ]
  rev$dist <- rev$V8 - rev$V3
  fow$dist <- fow$V2 - fow$V8
  rbind(rev, fow)
}

for_down <- function(mydf, mybed){
  if(is.null(mydf) || is.null(mybed)) return(data.frame())
  mydf2 <- mydf[, c("V4", "V8", "V11", "V12")]
  mybed2 <- mybed[, c("V2", "V3", "V4", "V6")]
  mymy <- merge(mydf2, mybed2, by.x=c("V4", "V12"), by.y=c("V4", "V6"))
  rev <- mymy[mymy$V12 == "-", ]
  fow <- mymy[mymy$V12 == "+", ]
  rev$dist <- rev$V2 - rev$V8
  fow$dist <- fow$V8 - fow$V3
  rbind(rev, fow)
}

find_border_ttest <- function(df, region, limit){
  sub <- df[df$region == region, ]
  border <- NA
  first_valid <- TRUE

  for(i in 0:limit){
    if(region == "upstream"){
      range_data <- sub[sub$dist >= -i & sub$dist <= 0, ]
    } else {
      range_data <- sub[sub$dist <= i & sub$dist >= 0, ]
    }

    low_vals <- range_data$V11[range_data$expr == "Lowly expressed genes"]
    high_vals <- range_data$V11[range_data$expr == "Highly expressed genes"]
    low_vals <- low_vals[is.finite(low_vals)]
    high_vals <- high_vals[is.finite(high_vals)]

    if(length(low_vals) < 5 || length(high_vals) < 5){
      next
    }

    mean_low <- mean(low_vals, na.rm=TRUE)
    mean_high <- mean(high_vals, na.rm=TRUE)
    test <- tryCatch(t.test(low_vals, high_vals), error=function(e) NULL)
    if(is.null(test)) next

    if(mean_low > mean_high && test$p.value < 0.05){
      border <- i
      first_valid <- FALSE
    } else {
      if(first_valid) return(NA)
      return(border)
    }
  }

  return(border)
}

linedf <- function(df){
  if(is.null(df) || nrow(df) == 0){
    return(data.frame(distance=numeric(), mC=numeric(), lower=numeric(), upper=numeric()))
  }

  df2 <- df[, c("V4", "V11", "dist")]
  df2$V11 <- as.numeric(df2$V11)
  out <- data.frame()

  for(i in 0:RAN2){
    left <- i * SS + 0.5
    right <- i * SS + WD
    sub <- df2[abs(df2$dist) > left & abs(df2$dist) <= right, ]
    sub <- sub[is.finite(sub$V11), ]

    if(nrow(sub) > 0){
      gene_mean <- tapply(sub$V11, sub$V4, mean, na.rm=TRUE)
      gene_mean <- gene_mean[is.finite(gene_mean)]

      if(length(gene_mean) >= 2){
        m <- mean(gene_mean, na.rm=TRUE)
        if(sd(gene_mean, na.rm=TRUE) == 0){
          lower <- m
          upper <- m
        } else {
          ci <- tryCatch(t.test(gene_mean)$conf.int, error=function(e) c(NA, NA))
          lower <- ci[1]
          upper <- ci[2]
        }
      } else if(length(gene_mean) == 1){
        m <- mean(gene_mean, na.rm=TRUE)
        lower <- m
        upper <- m
      } else {
        m <- NA
        lower <- NA
        upper <- NA
      }
    } else {
      m <- NA
      lower <- NA
      upper <- NA
    }

    center <- (left + right - 0.5) / 2
    out <- rbind(out, data.frame(distance=center, mC=m, lower=lower, upper=upper))
  }

  out
}

exprs <- c("low", "high")
dirs <- c("up", "down")

RAN <- limit
SS <- 4
Start <- WD / 2
RAN2 <- floor((RAN - WD / 2) / SS)

line_cols <- c(
  "Highly expressed genes" = "#B44A53",
  "Lowly expressed genes" = "#509ABC"
)

ci_cols <- c(
  "Highly expressed genes" = "#DEB9B7",
  "Lowly expressed genes" = "#ACD8EB"
)

for(type in types){
  message("Plotting stage=", stage, ", type=", type)

  data_list <- list()
  missing_any <- FALSE

  for(expr in exprs){
    for(dir in dirs){
      mydf_path <- paste0("wTE_", expr, "_", stage, "_", dir, "_", type, ".bed")
      mybed_path <- paste0(expr, "_", stage, ".txt")
      mydf <- safe_read(mydf_path)
      mybed <- safe_read(mybed_path)

      if(is.null(mydf) || is.null(mybed)){
        warning("Skipping empty/missing file: ", mydf_path)
        missing_any <- TRUE
        data_list[[paste0(expr, "_", dir)]] <- data.frame()
        next
      }

      if(ncol(mydf) < 12){
        warning("Skipping file with fewer than 12 columns: ", mydf_path)
        missing_any <- TRUE
        data_list[[paste0(expr, "_", dir)]] <- data.frame()
        next
      }

      if(dir == "up"){
        df <- for_up(mydf, mybed)
        if(nrow(df) > 0) df$dist <- 0 - df$dist
      } else {
        df <- for_down(mydf, mybed)
      }

      if(nrow(df) > 0){
        df$V11 <- as.numeric(df$V11)
        # If values are proportions, convert to percent. If already raw percentage/count-like values, keep them.
        if(all(is.na(df$V11) | (df$V11 >= 0 & df$V11 <= 1))){
          df$V11 <- df$V11 * 100
        }
      }

      data_list[[paste0(expr, "_", dir)]] <- df
    }
  }

  low_up <- data_list[["low_up"]]
  low_down <- data_list[["low_down"]]
  high_up <- data_list[["high_up"]]
  high_down <- data_list[["high_down"]]

  if(any(vapply(list(low_up, low_down, high_up, high_down), nrow, numeric(1)) == 0)){
    warning("Skipping plot for stage=", stage, ", type=", type, " because one or more required input sets are empty.")
    next
  }

  point_low <- rbind(low_up, low_down)
  point_low$expr <- "Lowly expressed genes"

  point_high <- rbind(high_up, high_down)
  point_high$expr <- "Highly expressed genes"

  point_all <- rbind(point_low, point_high)
  point_all <- point_all[is.finite(point_all$V11) & is.finite(point_all$dist), ]

  gap <- limit / 10

  up_df_point <- point_all[point_all$dist < 0, ]
  up_df_point$dist_shift <- up_df_point$dist

  down_df_point <- point_all[point_all$dist > 0, ]
  down_df_point$dist_shift <- down_df_point$dist + gap

  df_plot <- rbind(up_df_point, down_df_point)
  df_plot$region <- ifelse(df_plot$dist < 0, "upstream", "downstream")
  df_plot$expr <- factor(df_plot$expr, levels=c("Highly expressed genes", "Lowly expressed genes"))

  up_border <- find_border_ttest(df_plot, "upstream", limit)
  down_border <- find_border_ttest(df_plot, "downstream", limit)

  up_title <- if(!is.na(up_border)) paste0("Upstream border: ", -up_border, " bp") else "Upstream border: NA"
  dn_title <- if(!is.na(down_border)) paste0("Downstream border: ", down_border, " bp") else "Downstream border: NA"

  df_up_low <- linedf(low_up)
  df_down_low <- linedf(low_down)
  df_up_high <- linedf(high_up)
  df_down_high <- linedf(high_down)

  df_up_low$distance <- 0 - df_up_low$distance
  df_up_high$distance <- 0 - df_up_high$distance

  low_line <- rbind(df_up_low[, c("distance", "mC")], df_down_low[, c("distance", "mC")])
  low_line$expr <- "Lowly expressed genes"

  high_line <- rbind(df_up_high[, c("distance", "mC")], df_down_high[, c("distance", "mC")])
  high_line$expr <- "Highly expressed genes"

  line_all <- rbind(low_line, high_line)
  line_all <- line_all[complete.cases(line_all), ]

  low_CI_plot <- rbind(df_up_low, df_down_low)
  low_CI_plot$expr <- "Lowly expressed genes"

  high_CI_plot <- rbind(df_up_high, df_down_high)
  high_CI_plot$expr <- "Highly expressed genes"

  up_df_line <- line_all[line_all$distance < 0, ]
  up_df_line$distance_shift <- up_df_line$distance

  down_df_line <- line_all[line_all$distance > 0, ]
  down_df_line$distance_shift <- down_df_line$distance + gap

  up_low_CI_plot <- low_CI_plot[low_CI_plot$distance < 0 & complete.cases(low_CI_plot[, c("lower", "upper")]), ]
  up_low_CI_plot$distance_shift <- up_low_CI_plot$distance

  up_high_CI_plot <- high_CI_plot[high_CI_plot$distance < 0 & complete.cases(high_CI_plot[, c("lower", "upper")]), ]
  up_high_CI_plot$distance_shift <- up_high_CI_plot$distance

  down_low_CI_plot <- low_CI_plot[low_CI_plot$distance > 0 & complete.cases(low_CI_plot[, c("lower", "upper")]), ]
  down_low_CI_plot$distance_shift <- down_low_CI_plot$distance + gap

  down_high_CI_plot <- high_CI_plot[high_CI_plot$distance > 0 & complete.cases(high_CI_plot[, c("lower", "upper")]), ]
  down_high_CI_plot$distance_shift <- down_high_CI_plot$distance + gap

  breaks <- c(
    seq(-limit, -100, major_tick),
    0,
    gap,
    seq(gap + major_tick, gap + limit, major_tick)
  )

  labels <- c(
    seq(-limit, -100, major_tick),
    "TSS",
    "TTS",
    seq(major_tick, limit, major_tick)
  )

  vline_pos <- c()
  if(!is.na(up_border)) vline_pos <- c(vline_pos, -up_border)
  if(!is.na(down_border)) vline_pos <- c(vline_pos, down_border + gap)

  y_lab <- if(tolower(type) == "raw"){
    "TE raw methylation signal"
  } else {
    paste0("TE ", type, " methylation level (%)")
  }

  png(
    file=paste0("OUTPUT_3_TE_impact_distance_", stage, "_", type, ".png"),
    width=5000,
    height=2200,
    res=400
  )

  p <- ggplot() +
    geom_point(
      data=df_plot,
      mapping=aes(x=dist_shift, y=V11, color=expr),
      size=0.01,
      alpha=0.1,
      show.legend=FALSE
    ) +
    geom_ribbon(
      data=up_low_CI_plot,
      mapping=aes(x=distance_shift, ymin=lower, ymax=upper),
      fill=ci_cols["Lowly expressed genes"],
      alpha=0.27,
      show.legend=FALSE
    ) +
    geom_ribbon(
      data=down_low_CI_plot,
      mapping=aes(x=distance_shift, ymin=lower, ymax=upper),
      fill=ci_cols["Lowly expressed genes"],
      alpha=0.27,
      show.legend=FALSE
    ) +
    geom_ribbon(
      data=up_high_CI_plot,
      mapping=aes(x=distance_shift, ymin=lower, ymax=upper),
      fill=ci_cols["Highly expressed genes"],
      alpha=0.27,
      show.legend=FALSE
    ) +
    geom_ribbon(
      data=down_high_CI_plot,
      mapping=aes(x=distance_shift, ymin=lower, ymax=upper),
      fill=ci_cols["Highly expressed genes"],
      alpha=0.27,
      show.legend=FALSE
    ) +
    geom_line(
      data=up_df_line,
      mapping=aes(x=distance_shift, y=mC, color=expr, group=expr),
      linewidth=0.9
    ) +
    geom_line(
      data=down_df_line,
      mapping=aes(x=distance_shift, y=mC, color=expr, group=expr),
      linewidth=0.9
    ) +
    geom_hline(yintercept=c(0, 100), color="grey40", linewidth=0.3) +
    geom_vline(
      xintercept=vline_pos,
      color="grey50",
      linewidth=0.3,
      linetype=2
    ) +
    petem_theme_bw() +
    scale_color_manual(
      values=line_cols,
      guide=guide_legend(
        override.aes=list(
          linewidth=1.4,
          alpha=1,
          fill=NA
        )
      )
    ) +
    scale_x_continuous(breaks=breaks, labels=labels) +
    scale_y_continuous(breaks=c(0, 50, 100)) +
    ggtitle(paste0(up_title, "; ", dn_title)) +
    theme(
      legend.position="top",
      legend.title=element_blank(),
      legend.background=element_rect(fill="white", color=NA),
      legend.box.background=element_rect(fill="white", color=NA),
      legend.key=element_rect(fill="white", color=NA),
      legend.key.size=unit(1.2, "lines"),
      text=element_text(size=base_font_size),
      plot.title=element_text(size=base_font_size + 2, face="bold", hjust=0.5),
      axis.title=element_text(size=base_font_size + 2),
      axis.text=element_text(size=base_font_size),
      legend.text=element_text(size=base_font_size),
      panel.grid.minor=element_blank()
    ) +
    labs(
      x="Distance to gene (bp)",
      y=y_lab,
      color=NULL
    )

  print(p)
  dev.off()
}
EOF

    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec" | tee -a "$LOG"
}&
done
wait

end_step5=$(date +%s)
echo "[$(date)] Step 4 finished in $((end_step5-start_step5)) sec" | tee -a "$LOG"

rm -f wTE_*.bed low_* high_* expgene_closest_expTE.bed expressed_gene.bed expressed_TE.bed gene_sort.bed TE_sort.bed

end_allall=$(date +%s)
echo "[$(date)] Pipeline finished in $((end_allall-start_allall)) sec" | tee -a "$LOG"
