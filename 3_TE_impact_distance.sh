#!/usr/bin/env bash
set -euo pipefail

MODULE_NAME="Module 3"
SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
CURRENT_STEP="argument parsing"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Combined module 3: methylation preprocessing + TE impact distance plotting
# Usage:
#   bash 3_TE_impact_distance.sh \
#     -m sample1.CGmap.gz sample2.CGmap.gz ... \
#     -g gene.bed -t TE.bed \
#     -eg gene_expression.txt -et TE_expression.txt \
#     [-d 15000] [-tick 5000] [-p 10] [-w 100] [-l poly] [-unexp y|n] [-c y|n]

die() {
  echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: $*" >&2
  exit 1
}

usage() {
  echo "Usage: bash 3_TE_impact_distance.sh -m sample1.CGmap.gz [sample2.CGmap.gz ...] -g gene.bed -t TE.bed -eg expression_gene.txt -et expression_TE.txt [-d 15000] [-tick 5000] [-p 10] [-w 100] [-l poly] [-unexp y|n] [-c y|n]" >&2
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

collect_shared_stages() {
  python3 - "$1" "$2" <<'PYTHON_HELPER'
import sys
import pandas as pd

gene_exp = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
te_exp = pd.read_csv(sys.argv[2], sep="\t", index_col=0)
stages = sorted(set(gene_exp.columns).intersection(te_exp.columns))
for stage in stages:
    print(stage)
PYTHON_HELPER
}

preprocess_complete_for_stages() {
  local pre_dir=$1
  shift
  local stages=("$@")

  [[ -d "$pre_dir" ]] || return 1
  (( ${#stages[@]} > 0 )) || return 1

  local stage
  local context
  for stage in "${stages[@]}"; do
    for context in CG CHG CHH; do
      local path="${pre_dir}/pre3_${stage}_${context}.bed"
      [[ -s "$path" ]] || return 1
    done
  done
}

trap 'rc=$?; echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: command failed with exit code ${rc}" >&2; exit ${rc}' ERR

METH_FILES=()
GENE_BED=""
TE_BED=""
GENE_EXP=""
TE_EXP=""
LIMIT=15000
MAJOR_TICK=5000
TOP_BOTTOM_PERCENT=10
WINDOW_SIZE=100
LINE_MODE="raw"
INCLUDE_UNEXPRESSED_TE="n"
RUN_CONTROL_PLOT="n"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m)
      shift
      while [[ $# -gt 0 && ! $1 =~ ^- ]]; do
        METH_FILES+=("$1"); shift
      done
      ;;
    -g) GENE_BED="$2"; shift 2 ;;
    -t) TE_BED="$2"; shift 2 ;;
    -eg) GENE_EXP="$2"; shift 2 ;;
    -et) TE_EXP="$2"; shift 2 ;;
    -d|-lim) LIMIT="$2"; shift 2 ;;
    -tick) MAJOR_TICK="$2"; shift 2 ;;
    -p) TOP_BOTTOM_PERCENT="$2"; shift 2 ;;
    -w|-WD) WINDOW_SIZE="$2"; shift 2 ;;
    -l) LINE_MODE="$2"; shift 2 ;;
    -unexp) INCLUDE_UNEXPRESSED_TE="$2"; shift 2 ;;
    -c) RUN_CONTROL_PLOT="$2"; shift 2 ;;
    *) die "unknown option: $1" ;;
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
[[ $TOP_BOTTOM_PERCENT =~ ^[0-9]+([.][0-9]+)?$ ]] || die "invalid value for -p: ${TOP_BOTTOM_PERCENT}"
[[ $WINDOW_SIZE =~ ^[0-9]+$ ]] || die "invalid value for -w: ${WINDOW_SIZE}"
[[ $LIMIT =~ ^[0-9]+$ ]] || die "invalid value for -d: ${LIMIT}"
awk "BEGIN {exit !(${TOP_BOTTOM_PERCENT} > 0 && ${TOP_BOTTOM_PERCENT} < 50)}" || die "value for -p must be greater than 0 and less than 50: ${TOP_BOTTOM_PERCENT}"
(( WINDOW_SIZE > 0 )) || die "value for -w must be greater than 0: ${WINDOW_SIZE}"
(( LIMIT > 0 )) || die "value for -d must be greater than 0: ${LIMIT}"
[[ $LINE_MODE == "raw" || $LINE_MODE == "poly" ]] || die "invalid value for -l: ${LINE_MODE}. Supported values: poly"
for f in "${METH_FILES[@]}"; do
  require_file "-m" "$f"
done

CURRENT_STEP="step A0 - inspect preprocessing state"
mapfile -t SHARED_STAGES < <(collect_shared_stages "$GENE_EXP" "$TE_EXP")
(( ${#SHARED_STAGES[@]} > 0 )) || die "no shared stages found between $GENE_EXP and $TE_EXP"

SKIP_PRE=0
if preprocess_complete_for_stages "pre_step3" "${SHARED_STAGES[@]}"; then
  SKIP_PRE=1
  echo "[INFO] Found complete pre_step3 outputs for stages: ${SHARED_STAGES[*]}; module 3 will skip preprocessing."
fi

if [[ $SKIP_PRE -eq 0 && ${#METH_FILES[@]} -eq 0 ]]; then
  usage_missing "-m sample1.CGmap.gz [sample2.CGmap.gz ...]"
fi

LOG="LOG_3_TE_impact_distance.log"
echo "[`date`] Combined TE impact distance pipeline started" | tee "$LOG"
start_all=$(date +%s)

#####################################
# Part A: methylation preprocessing (from 3_1)
#####################################
if [[ $SKIP_PRE -eq 0 ]]; then
  CURRENT_STEP="step A1 - unzip and filter CGmap files"
  echo "[`date`] [A] Preprocess methylation: unzip + filter" | tee -a "$LOG"
  rm -rf pre_step3
  mkdir -p pre_step3

  for f in "${METH_FILES[@]}"; do
  (
    start=$(date +%s)
    base=$(basename "$f" .CGmap.gz)
    gunzip -c "$f" | awk '$8>=4 {print $1"\t"$3"\t"$2"\t"$4"\t"$6}' > "pre_step3/${base}.txt"
    end=$(date +%s)
    echo "[INFO] Preprocessed $f in $((end-start)) sec" | tee -a "$LOG"
  )&
  done
  wait
  echo "[INFO] All replicates done." | tee -a "$LOG"

  CURRENT_STEP="step A2 - compute per-stage average methylation"
  echo "[`date`] [A] Calculate average mC of each site at each stage" | tee -a "$LOG"
  stages=($(printf '%s\n' "${METH_FILES[@]}" | xargs -n1 basename | cut -d'_' -f1 | sort -u))
  for stage in "${stages[@]}"; do
  (
    start=$(date +%s)
    echo "[INFO] Processing stage $stage" | tee -a "$LOG"
    python3 - <<EOF
import pandas as pd, glob
stage = "${stage}"
files = glob.glob(f"pre_step3/{stage}_*.txt")
if not files:
    raise SystemExit(f"No files found for stage {stage}")

dfs = []
for f in files:
    df = pd.read_csv(f, sep="\t", header=None,
                     names=["chr", "site", "nt", "CNN", "mC"])
    dfs.append(df)

combined = pd.concat(dfs, axis=0, ignore_index=True)
m = combined.groupby(["chr", "site", "nt", "CNN"], as_index=False)["mC"].mean()
m = m.rename(columns={"mC": f"{stage}_mC"})

m["strand"] = m.nt.replace({"C": "+", "G": "-"})
m["name"] = ["site_" + str(i + 1) for i in range(len(m))]

for type in ["CG", "CHG", "CHH"]:
    sub = m[m.CNN == type][["chr", "site", "site", "name", f"{stage}_mC", "strand"]].dropna()
    sub.to_csv(f"pre_step3/pre3_{stage}_{type}.bed", sep="\t", index=False, header=False)
EOF
    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec" | tee -a "$LOG"
  )&
  done
  wait
  echo "[`date`] [A] Methylation preprocessing finished" | tee -a "$LOG"
else
  echo "[`date`] [A] Reusing existing pre_step3 files and skipping methylation preprocessing." | tee -a "$LOG"
  stages=("${SHARED_STAGES[@]}")
fi

#####################################
# Part B: plotting pipeline (from 3_2)
#####################################
CURRENT_STEP="step B0 - prepare sorted BED files"
echo "[`date`] [B] Input parsing" | tee -a "$LOG"
sort -k1,1 -k2,2n "$GENE_BED" > gene_sort.bed
sort -k1,1 -k2,2n "$TE_BED" > TE_sort.bed

CURRENT_STEP="step B1 - find nearby TEs"
echo "[`date`] [B] Step 1. Find nearby TEs for each gene" | tee -a "$LOG"
Rscript - "$GENE_EXP" "$TE_EXP"  "$INCLUDE_UNEXPRESSED_TE" <<'EOF' 

args <- commandArgs(trailingOnly=TRUE)
gene_exp <- args[1]
TE_exp   <- args[2]
include_unexp <- args[3]

gene_bed <- read.table("gene_sort.bed", header=F, stringsAsFactors=F)
TE_bed   <- read.table("TE_sort.bed", header=F, stringsAsFactors=F)
gene_exp    <- read.table(gene_exp, header=T, stringsAsFactors=F)
TE_exp    <- read.table(TE_exp, header=T, stringsAsFactors=F)

gene_exp <- gene_exp[rowSums(gene_exp) != 0, , drop=FALSE]
if(tolower(include_unexp) != "y"){
    TE_exp <- TE_exp[rowSums(TE_exp) != 0, , drop=FALSE]
}

gene_bed2 <- gene_bed[gene_bed$V4 %in% row.names(gene_exp), ]
TE_bed2   <- TE_bed[TE_bed$V4 %in% row.names(TE_exp), ]

write.table(gene_bed2, "expressed_gene.bed", sep="\t", quote=F, col.names=F, row.names=F)
write.table(TE_bed2, "expressed_TE.bed", sep="\t", quote=F, col.names=F, row.names=F)
EOF

bedtools closest -a expressed_gene.bed -b expressed_TE.bed -id -d -D a > expgene_closest_expTE.bed

CURRENT_STEP="step B2 - build expression strata and adjacent regions"
echo "[`date`] [B] Step 2. Split high/low genes and build adjacent regions" | tee -a "$LOG"
Rscript - "$GENE_EXP" "$LIMIT" "$TOP_BOTTOM_PERCENT" <<'EOF' 
args <- commandArgs(trailingOnly=TRUE)

gene_exp <- args[1]
limit <- as.numeric(args[2])
top_bottom_percent <- as.numeric(args[3])

gene_bed <- read.table("gene_sort.bed", header=F, stringsAsFactors=F)
gene_exp    <- read.table(gene_exp, header=T, stringsAsFactors=F)

# Function for adjacent regions
adjacent <- function(df, up=TRUE, limit=limit){
  df2 <- df
  if(up){
    df2[df2$V6=="+",3] <- df2[df2$V6=="+",2] - 1
    df2[df2$V6=="+",2] <- df2[df2$V6=="+",3] - (limit)
    df2[df2$V6=="-",2] <- df2[df2$V6=="-",3] + 1
    df2[df2$V6=="-",3] <- df2[df2$V6=="-",2] + (limit)
  } else {
    df2[df2$V6=="-",3] <- df2[df2$V6=="-",2] - 1
    df2[df2$V6=="-",2] <- df2[df2$V6=="-",3] - (limit)
    df2[df2$V6=="+",2] <- df2[df2$V6=="+",3] + 1
    df2[df2$V6=="+",3] <- df2[df2$V6=="+",2] + (limit)
  }
  df2$V2[df2$V2 < 1] <- 1
  df2[5] <- 0
  return(df2)
}

# Focus on TEs located within adjacent regions of genes
clo_TE <- read.table("expgene_closest_expTE.bed", sep="\t", header=F)
clo_TE2 <- clo_TE[(clo_TE$V10)!= ".", ]
clo_TE2 <- clo_TE2[clo_TE2$V13 >= (-limit), ]

gene_exp2 <- gene_exp[row.names(gene_exp) %in% clo_TE2$V4, ]

# Highly and lowly expressed genes for each stage
stages <- colnames(gene_exp2)
sink("OUTPUT_3_TE_impact_distance_gene_TE_number.txt")

for(stage in stages){
  vals <- gene_exp2[, stage, drop=FALSE]
  vals <- vals[vals[,1]>0, , drop=FALSE]
  if(nrow(vals) == 0){
    next
  }
  
  sorted <- vals[order(vals[,1]), , drop=FALSE]
  selected_n <- max(1, ceiling(nrow(sorted) * top_bottom_percent / 100))
  low  <- head(sorted, selected_n)
  high <- tail(sorted, selected_n)
  
  low_bed  <- merge(gene_bed, low, by.x= "V4", by.y= "row.names", all.y=T)
  high_bed <- merge(gene_bed, high, by.x= "V4", by.y= "row.names", all.y=T)
  low_bed <- low_bed[,c(2,3,4,1,7,6)]
  high_bed <- high_bed[,c(2,3,4,1,7,6)]

  write.table(high_bed, paste0("high_", stage, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(low_bed,  paste0("low_", stage, ".txt"),  row.names=F, col.names=F, quote=F, sep="\t")

  low_up   <- adjacent(low_bed, up=TRUE,  limit=limit)
  high_up  <- adjacent(high_bed, up=TRUE,  limit=limit)
  low_down <- adjacent(low_bed, up=FALSE, limit=limit)
  high_down<- adjacent(high_bed, up=FALSE, limit=limit)

  write.table(low_up,   paste0("low_",  stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
  write.table(high_up,  paste0("high_", stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
  write.table(low_down, paste0("low_",  stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(high_down,paste0("high_", stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")

  cat( stage, ":","\n")
  cat("  Highly expressed gene number:", nrow(high_bed), "\n")
  cat("  Lowly expressed gene number:", nrow(low_bed), "\n")
  
  high_genes <- high_bed$V4
  low_genes  <- low_bed$V4

  high_te <- length(unique(clo_TE2$V10[clo_TE2$V4 %in% high_genes]))
  low_te  <- length(unique(clo_TE2$V10[clo_TE2$V4 %in% low_genes]))
  
  cat("  TEs nearby highly expressed genes:", high_te, "\n")
  cat("  TEs nearby lowly expressed genes:", low_te,  "\n\n")
}
sink()
EOF

CURRENT_STEP="step B3 - intersect methylation tables"
echo "[`date`] [B] Step 3. Intersect TE methylation data" | tee -a "$LOG"
stages=($(ls *_up.bed | cut -d'_' -f2 | sort -u))
for stage in "${stages[@]}"; do
(
  start=$(date +%s)
  echo "[INFO] Processing stage $stage"   | tee -a "$LOG"

  missing_types=()
  for type in CG CHG CHH; do
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

          # 1. with TE
          wTE="wTE_${expr}_${stage}_${dir}.bed"
          bedtools intersect -a "$input" -b expressed_TE.bed > "$wTE"

          # 2. without TE
          woTE="woTE_${expr}_${stage}_${dir}.bed"
          bedtools intersect -a "$input" -b expressed_TE.bed -v > "$woTE"

          for type in CG CHG CHH; do
              bedtools intersect -a "$wTE" -b "pre_step3/pre3_${stage}_${type}.bed" -wa -wb > "wTE_${expr}_${stage}_${dir}_${type}.bed"
              bedtools intersect -a "$woTE" -b "pre_step3/pre3_${stage}_${type}.bed" -wa -wb > "woTE_${expr}_${stage}_${dir}_${type}.bed"
          done

          rm -f "$wTE" "$woTE"
      done
  done

  end=$(date +%s)
  echo "[INFO] Stage $stage done in $((end-start)) sec"   | tee -a "$LOG"
)&
done
wait

CURRENT_STEP="step B4 - generate plots"
echo "[`date`] [B] Step 4. Plotting" | tee -a "$LOG"
shopt -s nullglob
stage_files=(wTE_low_*_up_CG.bed)
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
else
for stage in "${stages[@]}"; do
{
    start=$(date +%s)
    echo "[INFO] Processing stage $stage"  | tee -a "$LOG"

    Rscript - "$LIMIT" "$MAJOR_TICK" "$WINDOW_SIZE" "$stage" "$SCRIPT_DIR" "$RUN_CONTROL_PLOT" "$LINE_MODE" <<'EOF'
library(ggplot2)
suppressPackageStartupMessages(library(gplots))

args <- commandArgs(trailingOnly=TRUE)
limit <- as.numeric(args[1])
major_tick <- as.numeric(args[2])
window_size <- as.numeric(args[3])
stage <- args[4]
script_dir <- args[5]
run_control_plot <- tolower(args[6])
line_mode <- tolower(args[7])
source(file.path(script_dir, "plot_defaults.R"))

# Functions
for_up <- function(mydf, mybed){
  mydf2 <- mydf[,c("V4","V8","V11","V12")]
  mybed2 <- mybed[,c("V2","V3","V4","V6")]
  mymy <- merge(mydf2, mybed2, by.x=c("V4","V12"), by.y=c("V4","V6"))
  rev <- mymy[mymy$V12=="-",]
  fow <- mymy[mymy$V12=="+",]
  rev$dist <- rev$V8 - rev$V3
  fow$dist <- fow$V2 - fow$V8
  return(rbind(rev,fow))
}

for_down <- function(mydf, mybed){
  mydf2 <- mydf[,c("V4","V8","V11","V12")]
  mybed2 <- mybed[,c("V2","V3","V4","V6")]
  mymy <- merge(mydf2, mybed2, by.x=c("V4","V12"), by.y=c("V4","V6"))
  rev <- mymy[mymy$V12=="-",]
  fow <- mymy[mymy$V12=="+",]
  rev$dist <- rev$V2 - rev$V8
  fow$dist <- fow$V8 - fow$V3
  return(rbind(rev,fow))
}

linedf <- function(df){
  df2 <- df[,c("V11","dist")]
  df3 <- aggregate(df2$V11, list(df2$dist), mean)
  colnames(df3) <- c("dist","V11")
  df_list <- c()
  for(i in 0:RAN2){
    val <- mean(df3[(abs(df3$dist) > (i * SS + 0.5)) & (abs(df3$dist) <= (i * SS + window_size)), ]$V11, na.rm=TRUE)
    df_list <- c(df_list,val)
  }
  if(length(df_list) == 0){
    return(df_list)
  }

  # Use a centered moving average across adjacent 100 bp bins.
  smoothed <- stats::filter(df_list, rep(1/3, 3), sides = 2)
  smoothed[is.na(smoothed)] <- df_list[is.na(smoothed)]
  as.numeric(smoothed)
}

types <- c("CG","CHG","CHH")
exprs <- c("low","high")
dirs <- c("up","down")

RAN <- limit
SS <- window_size
Start <- window_size / 2
RAN2 <- floor((RAN - window_size / 2) / SS)

plot_prefixes <- c("wTE")
if(run_control_plot == "y"){
  plot_prefixes <- c(plot_prefixes, "woTE")
}

for(prefix in plot_prefixes){
for(type in types){
  for(expr in exprs){
    for(dir in dirs){
      mydf <- read.table(paste0(prefix,"_",expr,"_",stage,"_",dir,"_",type,".bed"), sep="\t")
      mybed <- read.table(paste0(expr,"_",stage,".txt"), sep="\t")
      if(dir=="up"){
        df <- for_up(mydf,mybed)
        df$dist <- 0 - df$dist 
      } else {
        df <- for_down(mydf,mybed)
      }
      df$V11 <- df$V11*100
      assign(paste0(expr,"_",dir),df)
    }
  }

  point_low <- rbind(get("low_up"), get("low_down"))
  point_low$expr <- "Lowly expressed genes"
  point_high <- rbind(get("high_up"), get("high_down"))
  point_high$expr <- "Highly expressed genes"
  point_all <- rbind(point_low,point_high)

  df_up_low   <- linedf(get("low_up"))
  df_down_low <- linedf(get("low_down"))
  df_up_high  <- linedf(get("high_up"))
  df_down_high<- linedf(get("high_down"))

  up_low_line   <- data.frame(distance=c(seq(-Start, -RAN, -SS)), mC=df_up_low)
  up_high_line  <- data.frame(distance=c(seq(-Start, -RAN, -SS)), mC=df_up_high)
  down_low_line <- data.frame(distance=c(seq(Start, RAN, SS)), mC=df_down_low)
  down_high_line<- data.frame(distance=c(seq(Start, RAN, SS)), mC=df_down_high)

  low_line  <- rbind(up_low_line, down_low_line); low_line$expr <- "Lowly expressed genes"
  high_line <- rbind(up_high_line, down_high_line); high_line$expr <- "Highly expressed genes"
  line_all  <- rbind(low_line,high_line)
  line_all <- line_all[complete.cases(line_all),]

  # Check border
  check_border <- cbind(low_line, high_line)
  check_up <- check_border[check_border$distance < 0, ]
  check_dn <- check_border[check_border$distance > 0, ]
  border_up <- head(check_up[check_up[,2] < check_up[,5], ], 1)
  border_dn <- head(check_dn[check_dn[,2] < check_dn[,5], ], 1)

  up_title <- if(nrow(border_up)>0) paste0("Upstream border: ", border_up$distance[1], " bp") else ""
  dn_title <- if(nrow(border_dn)>0) paste0("Downstream border: ", border_dn$distance[1], " bp") else ""

  
  gap <- limit/10

  up_df_point <- point_all[point_all$dist < 0,]
  up_df_line  <- line_all[line_all$distance < 0,]
  up_df_point$dist_shift <- up_df_point$dist
  up_df_line$distance_shift <- up_df_line$distance

  down_df_point <- point_all[point_all$dist > 0,]
  down_df_line  <- line_all[line_all$distance > 0,]
  down_df_point$dist_shift <- down_df_point$dist + gap
  down_df_line$distance_shift <- down_df_line$distance + gap

  df_point_all <- rbind(up_df_point, down_df_point)
  
  left_breaks <- if(major_tick < limit) seq(-limit, -major_tick, major_tick) else -limit
  right_breaks <- if(major_tick < limit) seq(gap + major_tick, gap + limit, major_tick) else gap + limit
  left_labels <- if(major_tick < limit) seq(-limit, -major_tick, major_tick) else -limit
  right_labels <- if(major_tick < limit) seq(major_tick, limit, major_tick) else limit
  breaks <- c(left_breaks, 0, gap, right_breaks)
  labels <- c(left_labels, "TSS", "TTS", right_labels)

  output_file <- if(prefix == "woTE"){
    paste0("OUTPUT_3_TE_impact_distance_control_",stage,"_",type,".png")
  } else {
    paste0("OUTPUT_3_TE_impact_distance_",stage,"_",type,".png")
  }
  y_label <- if(prefix == "woTE") "Control methylation level (%)" else "TE methylation level (%)"

  png(file=output_file, width=5000, height=2000, res=400)

  p <- ggplot() +
  geom_point(df_point_all, mapping=aes(x=dist_shift,y=V11,color=expr), size=0.01, alpha=0.03)

  if(line_mode == "poly"){
    p <- p +
      geom_smooth(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                  method="lm", formula=y ~ poly(x, 3), se=TRUE, linewidth=0.8, alpha=0.18) +
      geom_smooth(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                  method="lm", formula=y ~ poly(x, 3), se=TRUE, linewidth=0.8, alpha=0.18)
  } else {
    p <- p +
      geom_line(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr)) +
      geom_line(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr))
  }

  p <- p +
  petem_theme_bw() +
  scale_color_manual(values=c("#B44A53","#509ABC")) +
  scale_x_continuous(breaks=breaks, labels=labels) +
  ggtitle(paste0(up_title, "; ", dn_title)) +
  theme(legend.position="top",
    panel.grid.minor = element_blank()) +
  labs(x="Distance to gene (bp)", y=y_label, color=NULL)

  print(p)

  dev.off()
}
}

EOF

    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec"  | tee -a "$LOG"
}&
done
wait
fi

rm wTE_*.bed woTE_*.bed low_* high_* expgene_closest_expTE.bed expressed_gene.bed expressed_TE.bed

end_all=$(date +%s)
echo "[`date`] Combined pipeline finished in $((end_all-start_all)) sec" | tee -a "$LOG"
