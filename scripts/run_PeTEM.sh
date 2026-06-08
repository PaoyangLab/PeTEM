#!/usr/bin/env bash
set -euo pipefail

MODULE_NAME="Pipeline Runner"
SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
CURRENT_STEP="initialization"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ENV_CHECK_SCRIPT="$SCRIPT_DIR/env_check.sh"

die() {
  echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: $*" >&2
  exit 1
}

trap 'rc=$?; echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: command failed with exit code ${rc}" >&2; exit ${rc}' ERR

abs_path() {
  python3 - "$1" <<'PYTHON_HELPER'
import os, sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PYTHON_HELPER
}

ask_yes_no() {
  local __var=$1
  local prompt=$2
  local default=${3:-y}
  local prompt_suffix
  if [[ ${default,,} == y ]]; then
    prompt_suffix=" (Y/n): "
    default=y
  else
    prompt_suffix=" (y/N): "
    default=n
  fi
  while true; do
    read -r -p "${prompt}${prompt_suffix}" reply || exit 1
    reply=${reply,,}
    if [[ -z $reply ]]; then
      reply=$default
    fi
    case $reply in
      y|yes)
        printf -v "${__var}" "y"
        return 0
        ;;
      n|no)
        printf -v "${__var}" "n"
        return 0
        ;;
      *)
        echo "Please answer y or n."
        ;;
    esac
  done
}

ask_file() {
  local __var=$1
  local prompt=$2
  local expected=${3:-}
  local module_label=${4:-General}
  while true; do
    read -r -p "$prompt" value || exit 1
    if [[ -z $value ]]; then
      echo "[ERROR] ${module_label} | ${SCRIPT_NAME} | ${CURRENT_STEP}: input is required."
      continue
    fi
    local resolved
    resolved=$(abs_path "$value")
    if [[ ! -f $resolved ]]; then
      echo "[ERROR] ${module_label} | ${SCRIPT_NAME} | ${CURRENT_STEP}: input file not found: $resolved"
      continue
    fi
    if [[ -n $expected ]]; then
      shopt -s nocasematch
      local IFS='|'
      read -r -a __suffixes <<< "$expected"
      local matched=false
      for suffix in "${__suffixes[@]}"; do
        if [[ $resolved == *${suffix} ]]; then
          matched=true
          break
        fi
      done
      shopt -u nocasematch
      if [[ $matched == false ]]; then
        echo "[ERROR] ${module_label} | ${SCRIPT_NAME} | ${CURRENT_STEP}: expected a file ending with ${expected}, got: $resolved"
        continue
      fi
    fi
    printf -v "${__var}" "%s" "$resolved"
    return 0
  done
}

ask_optional_list() {
  local __array=$1
  local -n __ref=$__array
  local prompt=$2
  local expected_suffix=${3:-}
  local module_label=${4:-General}
  while true; do
    local -a values
    read -r -a values -p "$prompt" || exit 1
    if [[ ${#values[@]} -eq 0 ]]; then
      echo "[ERROR] ${module_label} | ${SCRIPT_NAME} | ${CURRENT_STEP}: please provide at least one file."
      continue
    fi
    local -a resolved_list=()
    local valid=true
    for v in "${values[@]}"; do
      local resolved
      resolved=$(abs_path "$v")
      if [[ ! -f $resolved ]]; then
        echo "[ERROR] ${module_label} | ${SCRIPT_NAME} | ${CURRENT_STEP}: input file not found: $resolved"
        valid=false
        break
      fi
      if [[ -n $expected_suffix ]]; then
        shopt -s nocasematch
        local IFS='|'
        read -r -a __suffixes <<< "$expected_suffix"
        local matched=false
        for suffix in "${__suffixes[@]}"; do
          if [[ $resolved == *${suffix} ]]; then
            matched=true
            break
          fi
        done
        shopt -u nocasematch
        if [[ $matched == false ]]; then
          echo "[ERROR] ${module_label} | ${SCRIPT_NAME} | ${CURRENT_STEP}: expected files ending with ${expected_suffix}, got: $resolved"
          valid=false
          break
        fi
      fi
      resolved_list+=("$resolved")
    done
    if [[ $valid == true ]]; then
      __ref=("${resolved_list[@]}")
      return 0
    else
      echo "[ERROR] ${module_label} | ${SCRIPT_NAME} | ${CURRENT_STEP}: please re-enter the file list."
    fi
  done
}

ask_with_default() {
  local __var=$1
  local prompt=$2
  local default_value=$3
  while true; do
    read -r -p "$prompt" value || exit 1
    if [[ -z $value ]]; then
      printf -v "${__var}" "%s" "$default_value"
      return 0
    fi
    if [[ $value =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
      printf -v "${__var}" "%s" "$value"
      return 0
    else
      echo "Please provide a numeric value or press Enter for the default ($default_value)."
    fi
  done
}

run_environment_check() {
  [[ -x "$ENV_CHECK_SCRIPT" ]] || {
    die "environment check script not found or not executable: $ENV_CHECK_SCRIPT"
  }
  bash "$ENV_CHECK_SCRIPT"
}

module3_preprocess_complete() {
  local workdir=$1
  local gene_exp_file=$2
  local te_exp_file=$3

  [[ -d "$workdir/pre_step3" ]] || return 1
  [[ -f "$gene_exp_file" ]] || return 1
  [[ -f "$te_exp_file" ]] || return 1

  python3 - "$workdir" "$gene_exp_file" "$te_exp_file" <<'PYTHON_HELPER'
import os
import sys
import pandas as pd

workdir, gene_exp_path, te_exp_path = sys.argv[1:4]
pre_dir = os.path.join(workdir, "pre_step3")

gene_exp = pd.read_csv(gene_exp_path, sep="\t", index_col=0)
te_exp = pd.read_csv(te_exp_path, sep="\t", index_col=0)
stages = sorted(set(gene_exp.columns).intersection(te_exp.columns))

if not stages:
    sys.exit(1)

for stage in stages:
    for context in ("CG", "CHG", "CHH"):
        path = os.path.join(pre_dir, f"pre3_{stage}_{context}.bed")
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            sys.exit(1)
PYTHON_HELPER
}

#####################################
# Run pipeline interactively
#####################################

CURRENT_STEP="interactive module selection"
echo "Select modules to run (y/n):"
ask_yes_no run0 "Module 0. Preprocessing?" "y"
ask_yes_no run1 "Module 1. TE distribution?" "y"
ask_yes_no run2 "Module 2. Promoter-embedded TE families?" "y"
ask_yes_no run3 "Module 3. TE impact distance?" "y"
ask_yes_no run4 "Module 4. Single-condition correlation?" "y"
ask_yes_no run5 "Module 5. Cross-condition correlation?" "y"

CURRENT_STEP="data directory selection"
read -r -p "Data working directory (default $(pwd)): " data_dir
data_dir=${data_dir:-$(pwd)}
data_dir=$(abs_path "$data_dir")
mkdir -p "$data_dir"
echo "[INFO] Using data directory: $data_dir"

CURRENT_STEP="environment check"
run_environment_check

#####################################
# Collect input files (ask once for shared ones)
#####################################

gene_bed=""
te_bed=""
promoter_bed=""
faidx=""
deg_file=""
dete_file=""
gene_exp=""
te_exp=""
meth_files=()
cds_bed=""
exon_bed=""
utr5_bed=""
utr3_bed=""
te_family=""
unexp="n"
limit=""
tick=""
window=""
wd_num=""
smooth=""
ylim_cg=""
ylim_chg=""
ylim_chh=""
ylim_te_ch=""
ylim_te_cg=""
module0_dir="$data_dir/module_0"
module1_dir="$data_dir/module_1"
module2_dir="$data_dir/module_2"
module3_dir="$data_dir/module_3"
module4_dir="$data_dir/module_4"
module5_dir="$data_dir/module_5"

CURRENT_STEP="collect module inputs"
# BED files shared
if [[ "$run0" == "y" || "$run1" == "y" || "$run3" == "y" ]]; then
  ask_file gene_bed "Gene BED file: " ".bed|.bed.gz" "Modules 0/1/3"
fi

if [[ "$run0" == "y" || "$run1" == "y" || "$run2" == "y" || "$run3" == "y" ]]; then
  ask_file te_bed "TE BED file: " ".bed|.bed.gz" "Modules 0/1/2/3"
fi

# Promoter BED (required for module 0)
if [[ "$run0" == "y" ]]; then
  ask_file promoter_bed "Promoter BED file: " ".bed|.bed.gz" "Module 0"
fi

# Genome FASTA index file (shared across modules 0 and 1)
if [[ "$run0" == "y" || "$run1" == "y" ]]; then
  ask_file faidx "Genome fasta index (.fai): " ".fai" "Modules 0/1"
fi

# DEG + DETE files (run0, run3-2, run4)
if [[ "$run0" == "y" || "$run3" == "y" || "$run4" == "y" ]]; then
  ask_file deg_file "DEG file: " ".txt" "Modules 0/3/4"
  ask_file dete_file "DETE file: " ".txt" "Modules 0/3/4"
fi

# Convert DEG + DETE files to expression.txt (run0, run3-2, run4)
if [[ "$run0" == "y" || "$run3" == "y" || "$run4" == "y" ]]; then
  CURRENT_STEP="prepare expression matrices"
  gene_exp="$data_dir/gene_expression.txt"
  te_exp="$data_dir/TE_expression.txt"
  Rscript - "$deg_file" "$gene_exp" "$dete_file" "$te_exp" <<'R_HELPER'
args <- commandArgs(trailingOnly=TRUE)
extract_expression_only <- function(infile, outfile) {
  df <- read.table(infile, header=TRUE, row.names=1, sep="	", check.names=FALSE)
  keep <- !grepl("^(logFC_|PValue_|FDR_)", colnames(df))
  df_expr <- df[, keep, drop=FALSE]
  write.table(df_expr, file=outfile, sep="	", quote=FALSE, col.names=NA)
}
extract_expression_only(args[1], args[2])
extract_expression_only(args[3], args[4])
R_HELPER
  CURRENT_STEP="collect module inputs"
fi

# Include unexpressed TEs in the analyses (run3, run4, run5)
if [[ "$run3" == "y" || "$run4" == "y" || "$run5" == "y" ]]; then
  ask_yes_no unexp "Include unexpressed TEs?" "n"
fi

# Methylation (module 0 or module 3 when preprocessing is required)
module3_has_pre="n"
if [[ "$run3" == "y" ]] && module3_preprocess_complete "$data_dir" "$gene_exp" "$te_exp"; then
  module3_has_pre="y"
  echo "[INFO] Module 3 | ${SCRIPT_NAME} | input collection: found complete pre_step3 outputs in $data_dir; module 3 will go directly to plotting."
fi

if [[ "$run0" == "y" || ( "$run3" == "y" && "$module3_has_pre" != "y" ) ]]; then
  ask_optional_list meth_files "Methylation CGmap.gz files (space separated): " ".cgmap.gz" "Modules 0/3"
fi

# Module 1 inputs
if [[ "$run1" == "y" ]]; then
  ask_file cds_bed "CDS BED file: " ".bed|.bed.gz" "Module 1"
  ask_file exon_bed "Exon BED file: " ".bed|.bed.gz" "Module 1"
  ask_file utr5_bed "5'UTR BED file: " ".bed|.bed.gz" "Module 1"
  ask_file utr3_bed "3'UTR BED file: " ".bed|.bed.gz" "Module 1"
fi

# Module 2 inputs
if [[ "$run2" == "y" ]]; then
  ask_file te_family "TE family file: " ".txt" "Module 2"
fi

# Module 3 plotting parameters
if [[ "$run3" == "y" ]]; then
  ask_with_default limit "Limit up-/down-stream range (bp, default 15000): " 15000
  ask_with_default tick "Tick size (bp, default 5000): " 5000
  ask_with_default window "Window size (bp, default 200): " 200
fi

# Module 4 parameters
if [[ "$run4" == "y" ]]; then
  ask_with_default smooth "Sliding window smooth multiplier (1-5, default 3): " 3
  ask_with_default ylim_cg "y-axis limit gene exp vs TE/promoter mC (CG, default 50): " 50
  ask_with_default ylim_chg "y-axis limit gene exp vs TE/promoter mC (CHG, default 10): " 10
  ask_with_default ylim_chh "y-axis limit gene exp vs TE/promoter mC (CHH, default 10): " 10
  ask_with_default ylim_te_ch "y-axis limit TE exp vs TE mC (CH, default 15): " 15
  ask_with_default ylim_te_cg "y-axis limit TE exp vs TE mC (CG, default 30): " 30
fi

#####################################
# Execute steps
#####################################

run_step_0() {
  (
    mkdir -p "$module0_dir"
    cd "$module0_dir"
    bash "$SCRIPT_DIR/0_preprocessing.sh" -g "$gene_bed" -t "$te_bed" -p "$promoter_bed" -eg "$gene_exp" -et "$te_exp" -fai "$faidx" -m "${meth_files[@]}"
  )
}

if [[ "$run0" == "y" ]]; then
  run_step_0
fi

# Module 1
if [[ "$run1" == "y" ]]; then
  (
    mkdir -p "$module1_dir"
    cd "$module1_dir"
    bash "$SCRIPT_DIR/1_TE_distribution.sh"       -g "$gene_bed" -c "$cds_bed"       -5 "$utr5_bed" -e "$exon_bed" -3 "$utr3_bed"       -p "$promoter_bed" -t "$te_bed" -fai "$faidx"
  )
fi

# Module 2
if [[ "$run2" == "y" ]]; then
  overlap_path="$module0_dir/TE_overlap_promoter.bed"
  [[ ! -f "$overlap_path" ]] && die "Module 2 requires $overlap_path from module 0, but the file was not found"
  (
    mkdir -p "$module2_dir"
    cd "$module2_dir"
    Rscript "$SCRIPT_DIR/2_TE_families.R" -a "$te_bed" -i "$overlap_path" -T "$te_family" -o "$module2_dir"
  )
fi

# Module 3
if [[ "$run3" == "y" ]]; then
  (
    mkdir -p "$module3_dir"
    cd "$module3_dir"
    if [[ ${#meth_files[@]} -gt 0 ]]; then
      bash "$SCRIPT_DIR/3_TE_impact_distance.sh" -m "${meth_files[@]}" -g "$gene_bed" -t "$te_bed" -i "$module0_dir/TE_overlap_promoter.bed" -eg "$gene_exp" -et "$te_exp" -lim "$limit" -tick "$tick" -WD "$window" -unexp "$unexp"
    else
      bash "$SCRIPT_DIR/3_TE_impact_distance.sh" -g "$gene_bed" -t "$te_bed" -i "$module0_dir/TE_overlap_promoter.bed" -eg "$gene_exp" -et "$te_exp" -lim "$limit" -tick "$tick" -WD "$window" -unexp "$unexp"
    fi
  )
fi

# Module 4
if [[ "$run4" == "y" ]]; then
  (
    mkdir -p "$module4_dir"
    cd "$module4_dir"
    Rscript "$SCRIPT_DIR/4_single_condition_correlation.R" --eg "$gene_exp" --et "$te_exp" --module0-dir "$module0_dir" --outdir "$module4_dir" --unexp "$unexp" --smooth "$smooth" --ylim_CG "$ylim_cg" --ylim_CHG "$ylim_chg" --ylim_CHH "$ylim_chh" --ylim_TEexpTEmC_CH "$ylim_te_ch" --ylim_TEexpTEmC_CG "$ylim_te_cg"
  )
fi

# Module 5 (uses raw DEG/DETE directly)
if [[ "$run5" == "y" ]]; then
  (
    mkdir -p "$module5_dir"
    cd "$module5_dir"
    Rscript "$SCRIPT_DIR/5_cross_condition_correlation.R" --module0-dir "$module0_dir" --outdir "$module5_dir" --unexp "$unexp"
  )
fi

echo "[DONE] Pipeline completed successfully."
