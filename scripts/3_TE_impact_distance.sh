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
#     -g gene.bed -t TE.bed -i TE_overlap_promoter.bed \
#     -eg gene_expression.txt -et TE_expression.txt \
#     [-d 15000] [-tick auto] [-p 10] [-w 100] [-l raw|linear|poly|loess] [-CI y|n] [-border y|n] [-unexp y|n] [-nTE y|n]

die() {
  echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: $*" >&2
  exit 1
}

usage() {
  echo "Usage: bash 3_TE_impact_distance.sh -m sample1.CGmap.gz [sample2.CGmap.gz ...] -g gene.bed -t TE.bed -i TE_overlap_promoter.bed -eg expression_gene.txt -et expression_TE.txt [-d 15000] [-tick auto] [-p 10] [-w 100] [-l raw|linear|poly|loess] [-CI y|n] [-border y|n] [-unexp y|n] [-nTE y|n]" >&2
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
import re
import sys
import pandas as pd

gene_exp = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
te_exp = pd.read_csv(sys.argv[2], sep="\t", index_col=0)
de_prefix = re.compile(r"^(logFC_|PValue_|FDR_|dexp_|dTEexp_|PV_)")
stages = sorted(
    c for c in set(gene_exp.columns).intersection(te_exp.columns)
    if not de_prefix.match(c)
)
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

plot_cache_complete_for_stages() {
  local cache_dir=$1
  local run_control_plot=$2
  shift 2
  local stages=("$@")

  [[ -d "$cache_dir" ]] || return 1
  (( ${#stages[@]} > 0 )) || return 1

  local stage expr dir context prefix
  for stage in "${stages[@]}"; do
    for expr in low high; do
      [[ -s "${cache_dir}/${expr}_${stage}.txt" ]] || return 1
      [[ -s "${cache_dir}/ctrl_${expr}_${stage}.txt" ]] || return 1
    done
    for prefix in wTE woTE; do
      [[ $prefix == "woTE" && $run_control_plot != "y" ]] && continue
      for expr in low high; do
        for dir in up down; do
          for context in CG CHG CHH; do
            [[ -s "${cache_dir}/${prefix}_${expr}_${stage}_${dir}_${context}.bed" ]] || return 1
          done
        done
      done
    done
  done
}

clear_plot_temp() {
  rm -f wTE_*.bed woTE_*.bed low_*.txt high_*.txt ctrl_low_*.txt ctrl_high_*.txt low_*_up.bed low_*_down.bed high_*_up.bed high_*_down.bed ctrl_low_*_up.bed ctrl_low_*_down.bed ctrl_high_*_up.bed ctrl_high_*_down.bed
}

restore_plot_cache() {
  local cache_dir=$1
  clear_plot_temp
  cp -f "${cache_dir}"/*.txt "${cache_dir}"/*.bed .
}

save_plot_cache() {
  local cache_dir=$1
  rm -rf "$cache_dir"
  mkdir -p "$cache_dir"
  cp low_*.txt high_*.txt ctrl_low_*.txt ctrl_high_*.txt wTE_*.bed woTE_*.bed "$cache_dir"/ 2>/dev/null || true
}

make_plot_cache_key() {
  local pre_dir=$1
  shift
  python3 - "$GENE_BED" "$TE_BED" "$PROMOTER_EMBEDDED_BED" "$GENE_EXP" "$TE_EXP" \
    "$LIMIT" "$TOP_BOTTOM_PERCENT" "$INCLUDE_UNEXPRESSED_TE" "$RUN_CONTROL_PLOT" "$pre_dir" "$@" <<'PYTHON_HELPER'
import hashlib
import json
import os
import sys

gene_bed, te_bed, promoter_bed, gene_exp, te_exp, limit, pct, unexp, nte, pre_dir, *stages = sys.argv[1:]

def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def meta(path, content_hash=False):
    st = os.stat(path)
    out = {
        "path": os.path.abspath(path),
        "size": st.st_size,
        "mtime_ns": st.st_mtime_ns,
    }
    if content_hash:
        out["sha256"] = file_hash(path)
        out.pop("mtime_ns")
    return out

data = {
    "cache_version": "module3_closest_gene_filter_v3",
    "inputs": [meta(p, content_hash=True) for p in [gene_bed, te_bed, promoter_bed, gene_exp, te_exp]],
    "params": {
        "limit": limit,
        "top_bottom_percent": pct,
        "include_unexpressed_te": unexp,
        "run_control_plot": nte,
    },
    "stages": stages,
    "pre3": [],
}

for stage in stages:
    for context in ("CG", "CHG", "CHH"):
        data["pre3"].append(meta(os.path.join(pre_dir, f"pre3_{stage}_{context}.bed")))

print(hashlib.sha256(json.dumps(data, sort_keys=True).encode("utf-8")).hexdigest())
PYTHON_HELPER
}

trap 'rc=$?; echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: command failed with exit code ${rc}" >&2; exit ${rc}' ERR

METH_FILES=()
GENE_BED=""
TE_BED=""
PROMOTER_EMBEDDED_BED=""
GENE_EXP=""
TE_EXP=""
LIMIT=4000
MAJOR_TICK=0
TOP_BOTTOM_PERCENT=10
WINDOW_SIZE=100
LINE_MODE="raw"
SHOW_CI="n"
SHOW_BORDER="n"
INCLUDE_UNEXPRESSED_TE="n"
RUN_CONTROL_PLOT="y"

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
    -i) PROMOTER_EMBEDDED_BED="$2"; shift 2 ;;
    -eg) GENE_EXP="$2"; shift 2 ;;
    -et) TE_EXP="$2"; shift 2 ;;
    -d|-lim) LIMIT="$2"; shift 2 ;;
    -tick) MAJOR_TICK="$2"; shift 2 ;;
    -p) TOP_BOTTOM_PERCENT="$2"; shift 2 ;;
    -w|-WD) WINDOW_SIZE="$2"; shift 2 ;;
    -l) LINE_MODE="$2"; shift 2 ;;
    -CI) SHOW_CI="$2"; shift 2 ;;
    -border) SHOW_BORDER="$2"; shift 2 ;;
    -unexp) INCLUDE_UNEXPRESSED_TE="$2"; shift 2 ;;
    -nTE) RUN_CONTROL_PLOT="$2"; shift 2 ;;
    *) die "unknown option: $1" ;;
  esac
done

missing_args=()
[[ -n ${GENE_BED:-} ]] || missing_args+=("-g gene.bed")
[[ -n ${TE_BED:-} ]] || missing_args+=("-t TE.bed")
[[ -n ${PROMOTER_EMBEDDED_BED:-} ]] || missing_args+=("-i TE_overlap_promoter.bed")
[[ -n ${GENE_EXP:-} ]] || missing_args+=("-eg expression_gene.txt")
[[ -n ${TE_EXP:-} ]] || missing_args+=("-et expression_TE.txt")
(( ${#missing_args[@]} == 0 )) || usage_missing "${missing_args[*]}"

CURRENT_STEP="input validation"
require_file "-g" "$GENE_BED"
require_file "-t" "$TE_BED"
require_file "-i" "$PROMOTER_EMBEDDED_BED"
require_file "-eg" "$GENE_EXP"
require_file "-et" "$TE_EXP"
[[ $TOP_BOTTOM_PERCENT =~ ^[0-9]+([.][0-9]+)?$ ]] || die "invalid value for -p: ${TOP_BOTTOM_PERCENT}"
[[ $WINDOW_SIZE =~ ^[0-9]+$ ]] || die "invalid value for -w: ${WINDOW_SIZE}"
[[ $LIMIT =~ ^[0-9]+$ ]] || die "invalid value for -d: ${LIMIT}"
[[ $MAJOR_TICK =~ ^[0-9]+$ ]] || die "invalid value for -tick: ${MAJOR_TICK}"
awk "BEGIN {exit !(${TOP_BOTTOM_PERCENT} > 0 && ${TOP_BOTTOM_PERCENT} < 50)}" || die "value for -p must be greater than 0 and less than 50: ${TOP_BOTTOM_PERCENT}"
(( WINDOW_SIZE > 0 )) || die "value for -w must be greater than 0: ${WINDOW_SIZE}"
(( LIMIT > 0 )) || die "value for -d must be greater than 0: ${LIMIT}"
[[ $SHOW_CI == "y" || $SHOW_CI == "n" ]] || die "invalid value for -CI: ${SHOW_CI}. Supported values: y, n"
[[ $SHOW_BORDER == "y" || $SHOW_BORDER == "n" ]] || die "invalid value for -border: ${SHOW_BORDER}. Supported values: y, n"
[[ $LINE_MODE == "raw" || $LINE_MODE == "linear" || $LINE_MODE == "poly" || $LINE_MODE == "loess" ]] || die "invalid value for -l: ${LINE_MODE}. Supported values: raw, linear, poly, loess"
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
PLOT_CACHE_DIR="pre_step3/module3_plot_cache"
PLOT_CACHE_KEY_FILE="pre_step3/module3_plot_cache.key"
PLOT_CACHE_KEY=$(make_plot_cache_key "pre_step3" "${stages[@]}")
SKIP_PLOT_PRE=0
if [[ -f "$PLOT_CACHE_KEY_FILE" ]] &&
   [[ $(<"$PLOT_CACHE_KEY_FILE") == "$PLOT_CACHE_KEY" ]] &&
   plot_cache_complete_for_stages "$PLOT_CACHE_DIR" "$RUN_CONTROL_PLOT" "${stages[@]}"; then
  SKIP_PLOT_PRE=1
  CURRENT_STEP="step B0 - restore plot cache"
  echo "[`date`] [B] Reusing cached module 3 plot inputs for identical data-affecting parameters." | tee -a "$LOG"
  restore_plot_cache "$PLOT_CACHE_DIR"
fi

if [[ $SKIP_PLOT_PRE -eq 0 ]]; then
clear_plot_temp
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

de_prefix <- "^(logFC_|PValue_|FDR_|dexp_|dTEexp_|PV_)"
gene_stage_cols <- setdiff(colnames(gene_exp), grep(de_prefix, colnames(gene_exp), value=TRUE))
TE_stage_cols <- setdiff(colnames(TE_exp), grep(de_prefix, colnames(TE_exp), value=TRUE))
shared_stage_cols <- intersect(gene_stage_cols, TE_stage_cols)
if(length(shared_stage_cols) == 0){
    stop("No shared expression stage columns found between gene and TE expression files.")
}

if(tolower(include_unexp) != "y"){
    gene_exp <- gene_exp[rowSums(gene_exp[, shared_stage_cols, drop=FALSE], na.rm=TRUE) != 0, , drop=FALSE]
    TE_exp <- TE_exp[rowSums(TE_exp[, shared_stage_cols, drop=FALSE], na.rm=TRUE) != 0, , drop=FALSE]
}

gene_bed2 <- gene_bed[gene_bed$V4 %in% row.names(gene_exp), ]
TE_bed2   <- TE_bed[TE_bed$V4 %in% row.names(TE_exp), ]

write.table(gene_bed2, "expressed_gene.bed", sep="\t", quote=F, col.names=F, row.names=F)
write.table(TE_bed2, "expressed_TE.bed", sep="\t", quote=F, col.names=F, row.names=F)
EOF

bedtools closest -a expressed_gene.bed -b expressed_TE.bed -id -d -D a > expgene_closest_expTE.bed

CURRENT_STEP="step B2 - build expression strata and adjacent regions"
echo "[`date`] [B] Step 2. Split high/low genes and build adjacent regions" | tee -a "$LOG"
Rscript - "$GENE_EXP" "$LIMIT" "$TOP_BOTTOM_PERCENT" "$INCLUDE_UNEXPRESSED_TE" "${stages[@]}" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)

gene_exp <- args[1]
limit <- as.numeric(args[2])
top_bottom_percent <- as.numeric(args[3])
include_unexp <- args[4]
allowed_stages <- args[-c(1:4)]

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
  df2$V2 <- pmax(0, as.numeric(df2$V2))
  df2$V3 <- pmax(0, as.numeric(df2$V3))
  df2 <- df2[df2$V3 > df2$V2, , drop=FALSE]
  df2[5] <- 0
  return(df2)
}

clo_TE <- read.table("expgene_closest_expTE.bed", sep="\t", header=FALSE, stringsAsFactors=FALSE)
if(ncol(clo_TE) < 13){
  stop("expgene_closest_expTE.bed does not contain the expected closest TE distance column.")
}
clo_TE2 <- clo_TE[clo_TE$V10 != "." & clo_TE$V13 >= (-limit), , drop=FALSE]
nearby_gene_ids <- unique(clo_TE2$V4)
nearby_gene_ids <- nearby_gene_ids[!is.na(nearby_gene_ids) & nearby_gene_ids != "."]
gene_ids <- gene_bed$V4

plot_gene_ids <- intersect(gene_ids, intersect(row.names(gene_exp), nearby_gene_ids))
control_gene_ids <- setdiff(intersect(gene_ids, row.names(gene_exp)), nearby_gene_ids)

gene_exp_plot <- gene_exp[row.names(gene_exp) %in% plot_gene_ids, , drop = FALSE]
gene_exp_control <- gene_exp[row.names(gene_exp) %in% control_gene_ids, , drop = FALSE]

# Highly and lowly expressed genes for each stage
stages <- intersect(colnames(gene_exp_plot), allowed_stages)
if(length(stages) == 0){
  stop("No overlapping expression stages matched the methylation stage names.")
}
sink("OUTPUT_3_TE_impact_distance_gene_TE_number.txt")

for(stage in stages){
  vals_plot <- gene_exp_plot[, stage, drop=FALSE]
  vals_control <- gene_exp_control[, stage, drop=FALSE]
  if(tolower(include_unexp) != "y"){
    vals_plot <- vals_plot[vals_plot[,1]>0, , drop=FALSE]
    vals_control <- vals_control[vals_control[,1]>0, , drop=FALSE]
  }
  if(nrow(vals_plot) == 0){
    next
  }

  sorted_plot <- vals_plot[order(vals_plot[,1]), , drop=FALSE]
  selected_plot_n <- max(1, ceiling(nrow(sorted_plot) * top_bottom_percent / 100))
  low_plot  <- head(sorted_plot, selected_plot_n)
  high_plot <- tail(sorted_plot, selected_plot_n)

  low_plot_bed  <- merge(gene_bed, low_plot, by.x= "V4", by.y= "row.names", all.y=T)
  high_plot_bed <- merge(gene_bed, high_plot, by.x= "V4", by.y= "row.names", all.y=T)
  low_plot_bed <- low_plot_bed[,c(2,3,4,1,7,6)]
  high_plot_bed <- high_plot_bed[,c(2,3,4,1,7,6)]

  write.table(high_plot_bed, paste0("high_", stage, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(low_plot_bed,  paste0("low_", stage, ".txt"),  row.names=F, col.names=F, quote=F, sep="\t")

  low_up   <- adjacent(low_plot_bed, up=TRUE,  limit=limit)
  high_up  <- adjacent(high_plot_bed, up=TRUE,  limit=limit)
  low_down <- adjacent(low_plot_bed, up=FALSE, limit=limit)
  high_down<- adjacent(high_plot_bed, up=FALSE, limit=limit)

  write.table(low_up,   paste0("low_",  stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
  write.table(high_up,  paste0("high_", stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
  write.table(low_down, paste0("low_",  stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(high_down,paste0("high_", stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")

  if(nrow(vals_control) > 0){
    sorted_control <- vals_control[order(vals_control[,1]), , drop=FALSE]
    selected_control_n <- max(1, ceiling(nrow(sorted_control) * top_bottom_percent / 100))
    low_control  <- head(sorted_control, selected_control_n)
    high_control <- tail(sorted_control, selected_control_n)

    low_control_bed  <- merge(gene_bed, low_control, by.x= "V4", by.y= "row.names", all.y=T)
    high_control_bed <- merge(gene_bed, high_control, by.x= "V4", by.y= "row.names", all.y=T)
    low_control_bed <- low_control_bed[,c(2,3,4,1,7,6)]
    high_control_bed <- high_control_bed[,c(2,3,4,1,7,6)]

    write.table(high_control_bed, paste0("ctrl_high_", stage, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
    write.table(low_control_bed,  paste0("ctrl_low_", stage, ".txt"),  row.names=F, col.names=F, quote=F, sep="\t")

    ctrl_low_up   <- adjacent(low_control_bed, up=TRUE,  limit=limit)
    ctrl_high_up  <- adjacent(high_control_bed, up=TRUE,  limit=limit)
    ctrl_low_down <- adjacent(low_control_bed, up=FALSE, limit=limit)
    ctrl_high_down<- adjacent(high_control_bed, up=FALSE, limit=limit)

    write.table(ctrl_low_up,   paste0("ctrl_low_",  stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
    write.table(ctrl_high_up,  paste0("ctrl_high_", stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
    write.table(ctrl_low_down, paste0("ctrl_low_",  stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")
    write.table(ctrl_high_down,paste0("ctrl_high_", stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")
  }

  cat( stage, ":","\n")
  cat("  TE-nearby highly expressed gene number:", nrow(high_plot_bed), "\n")
  cat("  TE-nearby lowly expressed gene number:", nrow(low_plot_bed), "\n")
  cat("  No nearby TE expressed gene number:", nrow(vals_control), "\n\n")
}
sink()
EOF

CURRENT_STEP="step B3 - intersect methylation tables"
echo "[`date`] [B] Step 3. Intersect TE methylation data" | tee -a "$LOG"
shopt -s nullglob
stage_files=(high_*.txt)
shopt -u nullglob
stages=()
for file in "${stage_files[@]}"; do
  stage=$(basename "$file")
  stage=${stage#high_}
  stage=${stage%.txt}
  [[ -n $stage ]] && stages+=("$stage")
done
pids=()
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
          control_input="ctrl_${expr}_${stage}_${dir}.bed"

          # 1. with TE
          wTE="wTE_${expr}_${stage}_${dir}.bed"
          bedtools intersect -a "$input" -b expressed_TE.bed > "$wTE"

          # 2. without TE
          woTE="woTE_${expr}_${stage}_${dir}.bed"
          if [[ ! -f "$control_input" ]]; then
              echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: missing control gene region BED for stage ${stage}: ${control_input}" | tee -a "$LOG"
              exit 1
          fi
          bedtools intersect -a "$control_input" -b expressed_TE.bed -v > "$woTE"

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
pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "$pid"
done
save_plot_cache "$PLOT_CACHE_DIR"
printf '%s\n' "$PLOT_CACHE_KEY" > "$PLOT_CACHE_KEY_FILE"
fi

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

    Rscript - "$LIMIT" "$MAJOR_TICK" "$WINDOW_SIZE" "$stage" "$SCRIPT_DIR" "$RUN_CONTROL_PLOT" "$LINE_MODE" "$SHOW_CI" "$TOP_BOTTOM_PERCENT" "$SHOW_BORDER" <<'EOF'
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
show_ci <- tolower(args[8])
top_bottom_percent <- as.numeric(args[9])
show_border <- tolower(args[10])
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
  # smoothed <- stats::filter(df_list, rep(1/3, 3), sides = 2)
  # smoothed[is.na(smoothed)] <- df_list[is.na(smoothed)]
  # as.numeric(smoothed)
  as.numeric(df_list)
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
      bed_prefix <- if(prefix == "woTE") "ctrl_" else ""
      mybed <- read.table(paste0(bed_prefix, expr,"_",stage,".txt"), sep="\t")
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

  low_label <- paste0("Bottom ", format(top_bottom_percent, trim = TRUE), "% lowly expressed genes")
  high_label <- paste0("Top ", format(top_bottom_percent, trim = TRUE), "% highly expressed genes")

  df_up_low   <- linedf(get("low_up"))
  df_down_low <- linedf(get("low_down"))
  df_up_high  <- linedf(get("high_up"))
  df_down_high<- linedf(get("high_down"))

  up_low_line   <- data.frame(distance=c(seq(-Start, -RAN, -SS)), mC=df_up_low)
  up_high_line  <- data.frame(distance=c(seq(-Start, -RAN, -SS)), mC=df_up_high)
  down_low_line <- data.frame(distance=c(seq(Start, RAN, SS)), mC=df_down_low)
  down_high_line<- data.frame(distance=c(seq(Start, RAN, SS)), mC=df_down_high)

  low_line  <- rbind(up_low_line, down_low_line); low_line$expr <- low_label
  high_line <- rbind(up_high_line, down_high_line); high_line$expr <- high_label
  line_all  <- rbind(low_line,high_line)
  line_all <- line_all[complete.cases(line_all),]

  # Check border
  check_border <- cbind(low_line, high_line)
  check_up <- check_border[check_border$distance < 0, ]
  check_dn <- check_border[check_border$distance > 0, ]
  border_up <- head(check_up[check_up[,2] < check_up[,5], ], 1)
  border_dn <- head(check_dn[check_dn[,2] < check_dn[,5], ], 1)

  format_border <- function(x){
    x <- as.numeric(x)
    ifelse(abs(x) <= 50, 0, x)
  }
  border_up_value <- if(nrow(border_up)>0) format_border(border_up$distance[1]) else NA
  border_dn_value <- if(nrow(border_dn)>0) format_border(border_dn$distance[1]) else NA
  up_title <- if(!is.na(border_up_value)) paste0("Upstream border: ", border_up_value, " bp") else ""
  dn_title <- if(!is.na(border_dn_value)) paste0("Downstream border: ", border_dn_value, " bp") else ""
  border_subtitle <- if(show_border == "y") paste0(up_title, "; ", dn_title) else NULL
  border_lines <- data.frame(x = numeric(0))

  
  gap <- limit/10
  up_df_line  <- line_all[line_all$distance < 0,]
  up_df_line$distance_shift <- up_df_line$distance

  down_df_line  <- line_all[line_all$distance > 0,]
  down_df_line$distance_shift <- down_df_line$distance + gap

  if (major_tick <= 0) {
    tick_values <- unique(round(seq(limit / 4, limit, by = limit / 4)))
  } else {
    tick_values <- unique(seq(major_tick, limit, by = major_tick))
  }
  left_breaks <- -rev(tick_values)
  right_breaks <- gap + tick_values
  left_labels <- -rev(tick_values)
  right_labels <- tick_values
  breaks <- c(left_breaks, 0, gap, right_breaks)
  labels <- c(left_labels, "TSS", "TES", right_labels)

  if(show_border == "y"){
    border_x <- c()
    if(nrow(border_up) > 0){
      border_x <- c(border_x, border_up_value)
    }
    if(nrow(border_dn) > 0){
      border_x <- c(border_x, border_dn_value + gap)
    }
    border_lines <- data.frame(x = border_x)
  }

  output_file <- if(prefix == "woTE"){
    paste0("OUTPUT_3_nonTE_impact_distance_",stage,"_",type,".png")
  } else {
    paste0("OUTPUT_3_TE_impact_distance_",stage,"_",type,".png")
  }
  y_label <- if(prefix == "woTE"){
    paste0("non-TE ", type, " methylation level (%)")
  } else {
    paste0("TE ", type, " methylation level (%)")
  }

  plot_prefix_label <- if(prefix == "woTE") "Non-TE" else "TE"
  plot_title <- paste0(plot_prefix_label, " ", type, " methylation around genes with different expression level")

  if(prefix == "woTE"){
    ymax <- max(line_all$mC, na.rm = TRUE)
    if(ymax < 5){
      y_limit <- 5
      y_breaks <- c(0,1,2,3,4,5)
    } else if(ymax < 20){
      y_limit  <- 20
      y_breaks <- c(0,5,10,15,20)
    } else if(ymax < 50){
      y_limit  <- 50
      y_breaks <- c(0,10,20,30,40,50)
    } else {
      y_limit  <- 100
      y_breaks <- c(0,20,40,60,80,100)
    }  
  } else {
    y_limit  <- 100
    y_breaks <- c(0,20,40,60,80,100)
  }

  png(file=output_file, width=3400, height=2000, res=400)

  p <- ggplot()

  if(line_mode == "poly"){
    if(show_ci == "y"){
      p <- p +
        geom_smooth(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,fill=expr,group=expr),
                    method="lm", formula=y ~ poly(x, 3), se=TRUE, linewidth=1.1, alpha=0.25) +
        geom_smooth(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,fill=expr,group=expr),
                    method="lm", formula=y ~ poly(x, 3), se=TRUE, linewidth=1.1, alpha=0.25)
    } else {
      p <- p +
        geom_smooth(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                    method="lm", formula=y ~ poly(x, 3), se=FALSE, linewidth=1.1) +
        geom_smooth(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                    method="lm", formula=y ~ poly(x, 3), se=FALSE, linewidth=1.1)
    }
  } else if(line_mode == "linear"){
    if(show_ci == "y"){
      p <- p +
        geom_smooth(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,fill=expr,group=expr),
                    method="lm", formula=y ~ x, se=TRUE, linewidth=1.1, alpha=0.25) +
        geom_smooth(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,fill=expr,group=expr),
                    method="lm", formula=y ~ x, se=TRUE, linewidth=1.1, alpha=0.25)
    } else {
      p <- p +
        geom_smooth(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                    method="lm", formula=y ~ x, se=FALSE, linewidth=1.1) +
        geom_smooth(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                    method="lm", formula=y ~ x, se=FALSE, linewidth=1.1)
    }
  } else if(line_mode == "loess"){
    if(show_ci == "y"){
      p <- p +
        geom_smooth(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,fill=expr,group=expr),
                    method="loess", se=TRUE, linewidth=1.1, alpha=0.25, span=0.4) +
        geom_smooth(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,fill=expr,group=expr),
                    method="loess", se=TRUE, linewidth=1.1, alpha=0.25, span=0.4)
    } else {
      p <- p +
        geom_smooth(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                    method="loess", se=FALSE, linewidth=1.1, span=0.4) +
        geom_smooth(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr),
                    method="loess", se=FALSE, linewidth=1.1, span=0.4)
    }
  } else {
    p <- p +
      geom_line(data=up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr), linewidth=1.1) +
      geom_line(data=down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr), linewidth=1.1)
  }

  p <- p +
  petem_theme_bw() +
  scale_color_manual(values=setNames(c("#B44A53","#509ABC"), c(high_label, low_label)),
                     breaks=c(high_label, low_label)) +
  scale_fill_manual(values=setNames(c("#F8CDD5","#CFE6FA"), c(high_label, low_label)),
                    breaks=c(high_label, low_label), guide = "none") +
  scale_x_continuous(breaks=breaks, labels=labels) +
  scale_y_continuous(limits=c(0, y_limit), breaks=y_breaks) +
  theme(legend.position="top",
    legend.text = element_text(size = PETEM_LEGEND_TEXT_SIZE + 3),
    panel.grid.minor = element_blank(),
    plot.title=element_text(size=PETEM_LEGEND_TEXT_SIZE + 4,face="bold",hjust=0.5)) +
  labs(title=plot_title, subtitle=border_subtitle, x="Distance to gene (bp)", y=y_label, color=NULL)

  if(nrow(border_lines) > 0){
    p <- p + geom_vline(data = border_lines, aes(xintercept = x),
                        inherit.aes = FALSE, linetype = "dashed",
                        linewidth = 0.8, color = "gray35")
  }

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

rm -f wTE_*.bed woTE_*.bed low_* high_* ctrl_low_* ctrl_high_* expgene_closest_expTE.bed expressed_gene.bed expressed_TE.bed

end_all=$(date +%s)
echo "[`date`] Combined pipeline finished in $((end_all-start_all)) sec" | tee -a "$LOG"
