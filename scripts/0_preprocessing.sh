#!/usr/bin/env bash
set -euo pipefail

MODULE_NAME="Module 0"
SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
CURRENT_STEP="initialization"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BIN_DIR="$SCRIPT_DIR/bin"
WIG_TO_BIGWIG="$BIN_DIR/wigToBigWig"
BIGWIG_AVERAGE="$BIN_DIR/bigWigAverageOverBed"
CGMAP2WIGGLE="$SCRIPT_DIR/CGmap2Wiggle.pl"

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

require_executable() {
    local label=$1
    local path=$2
    [[ -x $path ]] || die "required executable not found or not executable (${label}): ${path}"
}

trap 'rc=$?; echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: command failed with exit code ${rc}" >&2; exit ${rc}' ERR

require_executable "wigToBigWig" "$WIG_TO_BIGWIG"
require_executable "bigWigAverageOverBed" "$BIGWIG_AVERAGE"
[[ -f "$CGMAP2WIGGLE" ]] || die "required file not found (CGmap2Wiggle.pl): $CGMAP2WIGGLE"

#####################################
# Usage
#####################################
usage() {
    echo "Usage: bash 0_preprocessing.sh -g gene.bed -t TE.bed -p promoter.bed -eg expression_gene.txt -et expression_TE.txt -fai genome.fa.fai {-m sample1.CGmap.gz [sample2.CGmap.gz ...] | -m1 stage1_rep.CGmap.gz [...] -m2 stage2_rep.CGmap.gz [...] [-m3 ... -m4 ... -m5 ...]}" >&2
    exit 1
}

expression_stage_order() {
    python3 - "$1" "$2" <<'PYTHON_HELPER'
import re
import sys

gene_exp, te_exp = sys.argv[1:3]
de_prefix = re.compile(r"^(logFC_|PValue_|FDR_|dexp_|dTEexp_|PV_|DEG_)")
id_cols = {"", "id", "gene", "gene_id", "te", "te_id", "transcript_id"}

def stages(path):
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header and header[0].strip().lower() in id_cols:
        header = header[1:]
    out = []
    seen = set()
    for col in header:
        if not col or de_prefix.match(col):
            continue
        stage = col.split("_", 1)[0]
        if stage not in seen:
            out.append(stage)
            seen.add(stage)
    return out

gene_stages = stages(gene_exp)
te_stages = stages(te_exp)
if not gene_stages:
    raise SystemExit(f"no expression stage columns found in {gene_exp}")
if gene_stages != te_stages:
    raise SystemExit(
        "gene and TE expression stage columns must have the same order: "
        f"gene={gene_stages}; TE={te_stages}"
    )
for stage in gene_stages:
    print(stage)
PYTHON_HELPER
}

#####################################
# Parse args
#####################################
GENE=""
TE=""
PROMOTER=""
GENE_EXP=""
TE_EXP=""
faidx=""
METH_FILES=()
METH_GROUP1=()
METH_GROUP2=()
METH_GROUP3=()
METH_GROUP4=()
METH_GROUP5=()
SPLIT_METH_USED=0

while [[ $# -gt 0 ]]; do
    case $1 in
        -g) GENE=$2; shift 2;;
        -t) TE=$2; shift 2;;
        -p) PROMOTER=$2; shift 2;;
        -eg) GENE_EXP=$2; shift 2;;
        -et) TE_EXP=$2; shift 2;;
        -fai) faidx=$2; shift 2;;
        -m)
            shift
            while [[ $# -gt 0 ]] && [[ ! $1 =~ ^- ]]; do
                METH_FILES+=("$1")
                shift
            done
            ;;
        -m1)
            SPLIT_METH_USED=1
            shift
            while [[ $# -gt 0 ]] && [[ ! $1 =~ ^- ]]; do
                METH_GROUP1+=("$1")
                shift
            done
            ;;
        -m2)
            SPLIT_METH_USED=1
            shift
            while [[ $# -gt 0 ]] && [[ ! $1 =~ ^- ]]; do
                METH_GROUP2+=("$1")
                shift
            done
            ;;
        -m3)
            SPLIT_METH_USED=1
            shift
            while [[ $# -gt 0 ]] && [[ ! $1 =~ ^- ]]; do
                METH_GROUP3+=("$1")
                shift
            done
            ;;
        -m4)
            SPLIT_METH_USED=1
            shift
            while [[ $# -gt 0 ]] && [[ ! $1 =~ ^- ]]; do
                METH_GROUP4+=("$1")
                shift
            done
            ;;
        -m5)
            SPLIT_METH_USED=1
            shift
            while [[ $# -gt 0 ]] && [[ ! $1 =~ ^- ]]; do
                METH_GROUP5+=("$1")
                shift
            done
            ;;
        *) die "unknown option: $1" ;;
    esac
done

CURRENT_STEP="argument validation"
missing_args=()
[[ -n ${GENE:-} ]] || missing_args+=("-g gene.bed")
[[ -n ${TE:-} ]] || missing_args+=("-t TE.bed")
[[ -n ${PROMOTER:-} ]] || missing_args+=("-p promoter.bed")
[[ -n ${GENE_EXP:-} ]] || missing_args+=("-eg expression_gene.txt")
[[ -n ${TE_EXP:-} ]] || missing_args+=("-et expression_TE.txt")
[[ -n ${faidx:-} ]] || missing_args+=("-fai genome.fa.fai")
if [[ $SPLIT_METH_USED -eq 1 && ${#METH_FILES[@]} -gt 0 ]]; then
    die "use either -m or -m1/-m2/-m3/-m4/-m5, not both"
fi
if [[ $SPLIT_METH_USED -eq 0 && ${#METH_FILES[@]} -eq 0 ]]; then
    missing_args+=("-m sample1.CGmap.gz [sample2.CGmap.gz ...] or -m1/-m2...")
fi
(( ${#missing_args[@]} == 0 )) || usage_missing "${missing_args[*]}"

CURRENT_STEP="input validation"
require_file "-g" "$GENE"
require_file "-t" "$TE"
require_file "-p" "$PROMOTER"
require_file "-eg" "$GENE_EXP"
require_file "-et" "$TE_EXP"
require_file "-fai" "$faidx"

mapfile -t EXPR_STAGES < <(expression_stage_order "$GENE_EXP" "$TE_EXP")
METH_INPUT_MODE="-m"
METH_MAP_LINES=()
if [[ $SPLIT_METH_USED -eq 1 ]]; then
    METH_INPUT_MODE="-m1..-m5"
    (( ${#EXPR_STAGES[@]} <= 5 )) || die "-m1..-m5 supports at most 5 expression stages; found ${#EXPR_STAGES[@]}"
    group_sizes=("${#METH_GROUP1[@]}" "${#METH_GROUP2[@]}" "${#METH_GROUP3[@]}" "${#METH_GROUP4[@]}" "${#METH_GROUP5[@]}")
    for ((i = 0; i < 5; i++)); do
        if (( i < ${#EXPR_STAGES[@]} )); then
            (( group_sizes[$i] > 0 )) || die "missing methylation files for expression stage '${EXPR_STAGES[$i]}' (expected -m$((i+1)))"
        else
            (( group_sizes[$i] == 0 )) || die "-m$((i+1)) was provided but expression table has only ${#EXPR_STAGES[@]} stages"
        fi
    done
    GROUPED_METH_DIR="$PWD/tmp/module_0_methylation_inputs"
    rm -rf "$GROUPED_METH_DIR"
    mkdir -p "$GROUPED_METH_DIR"
    METH_FILES=()
    for ((i = 0; i < ${#EXPR_STAGES[@]}; i++)); do
        stage="${EXPR_STAGES[$i]}"
        [[ $stage =~ ^[A-Za-z0-9.-]+$ ]] || die "stage name '$stage' is not supported for -m1..-m5; use letters/numbers/dot/hyphen only"
        current_group=()
        case $((i+1)) in
            1) current_group=("${METH_GROUP1[@]}") ;;
            2) current_group=("${METH_GROUP2[@]}") ;;
            3) current_group=("${METH_GROUP3[@]}") ;;
            4) current_group=("${METH_GROUP4[@]}") ;;
            5) current_group=("${METH_GROUP5[@]}") ;;
        esac
        rep=1
        for original in "${current_group[@]}"; do
            require_file "Methylation file for ${stage}" "$original"
            prepared="$GROUPED_METH_DIR/${stage}_$(printf '%02d' "$rep").CGmap.gz"
            ln -sf "$(readlink -f "$original")" "$prepared"
            METH_FILES+=("$prepared")
            METH_MAP_LINES+=("${stage}"$'\t'"-m$((i+1))"$'\t'"${rep}"$'\t'"$(readlink -f "$original")"$'\t'"${prepared}")
            rep=$((rep + 1))
        done
    done
else
    for i in "${!METH_FILES[@]}"; do
        require_file "-m" "${METH_FILES[$i]}"
        base=$(basename "${METH_FILES[$i]}" .CGmap.gz)
        stage=${base%%_*}
        METH_MAP_LINES+=("${stage}"$'\t'"-m"$'\t'"$((i+1))"$'\t'"$(readlink -f "${METH_FILES[$i]}")"$'\t'"$(readlink -f "${METH_FILES[$i]}")")
    done
fi

{
    echo "PeTEM Module 0 input mapping"
    echo "created_at_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "gene_expression (-eg): $(readlink -f "$GENE_EXP")"
    echo "TE_expression (-et): $(readlink -f "$TE_EXP")"
    echo "methylation_input_mode: $METH_INPUT_MODE"
    echo
    echo "expression_stage_order:"
    for i in "${!EXPR_STAGES[@]}"; do
        printf '%s\t%s\n' "$((i+1))" "${EXPR_STAGES[$i]}"
    done
    echo
    echo -e "stage\tm_flag\treplicate_index\toriginal_cgmap\tprepared_cgmap"
    printf '%s\n' "${METH_MAP_LINES[@]}"
} > "$PWD/LOG_0_input_mapping.log"

echo "[INFO] Gene BED: $GENE"
echo "[INFO] TE BED:   $TE"
echo "[INFO] Promoter BED: $PROMOTER"
echo "[INFO] Gene expression: $GENE_EXP"
echo "[INFO] TE expression:   $TE_EXP"
echo "[INFO] Genome fasta index: $faidx"
echo "[INFO] Methylation CGmaps: ${METH_FILES[*]}"

#####################################
# Timer start
#####################################
START=$(date +%s)
LOG="$PWD/LOG_0_preprocessing.log"
: > "$LOG"

#####################################
# 1. Annotations preprocessing
#####################################

CURRENT_STEP="step 1 - preprocess annotations"
echo "[STEP 1] Preprocessing annotations..."

# TE overlap with promoters (using provided promoter BED)
bedtools intersect -a "$TE" -b "$PROMOTER" -wa -wb > TE_overlap_promoter.bed

# count TE overlap with promoters
Rscript - "$TE" "$GENE" "$TE_EXP" "$GENE_EXP" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
te_file <- args[1]
gene_file <- args[2]
te_exp_file <- args[3]
gene_exp_file <- args[4]

te_bed <- read.table(te_file, sep="\t", header=FALSE)
gene_bed <- read.table(gene_file, sep="\t", header=FALSE)
if (file.info("TE_overlap_promoter.bed")$size > 0) {
  overlap <- read.table("TE_overlap_promoter.bed", sep="\t", header=FALSE)
} else {
  overlap <- data.frame()
}

te_exp <- read.table(te_exp_file, header=TRUE, row.names=1, check.names=FALSE)
gene_exp <- read.table(gene_exp_file, header=TRUE, row.names=1, check.names=FALSE)

te_exp <- te_exp[rowSums(te_exp) != 0, , drop=FALSE]
gene_exp <- gene_exp[rowSums(gene_exp) != 0, , drop=FALSE]

expressed_te <- rownames(te_exp)
expressed_gene <- rownames(gene_exp)

# TE table
all_te <- te_bed$V4
embedded_te <- if (nrow(overlap) > 0) unique(overlap$V4) else character(0)

a <- sum(all_te %in% intersect(embedded_te, expressed_te))
b <- sum(all_te %in% setdiff(embedded_te, expressed_te))
c <- sum(all_te %in% setdiff(expressed_te, embedded_te))
d <- sum(all_te %in% setdiff(all_te, union(embedded_te, expressed_te)))

te_tab <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE,
                 dimnames=list(Embedded=c("Yes","No"), Expressed=c("Yes","No")))
te_chi <- chisq.test(te_tab)

# Gene table
all_gene <- gene_bed$V4
embedded_gene <- if (nrow(overlap) > 0) unique(overlap$V10) else character(0)

a <- sum(all_gene %in% intersect(embedded_gene, expressed_gene))
b <- sum(all_gene %in% setdiff(embedded_gene, expressed_gene))
c <- sum(all_gene %in% setdiff(expressed_gene, embedded_gene))
d <- sum(all_gene %in% setdiff(all_gene, union(embedded_gene, expressed_gene)))

gene_tab <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE,
                   dimnames=list(Promoter_has_TE=c("Yes","No"), Expressed=c("Yes","No")))
gene_chi <- chisq.test(gene_tab)

# Output txt
sink("OUTPUT_0_embedded_TE_gene_number.txt")
cat("########## TE expression vs promoter embedding ##########\n\n")
print(te_tab)
cat("\nChi-squared p =", te_chi$p.value, "\n\n\n\n")

cat("########## Gene expression vs promoter has TE ##########\n\n")
print(gene_tab)
cat("\nChi-squared p =", gene_chi$p.value, "\n")
sink()
EOF



# Promoter regions with TE insertions
bedtools intersect -a "$PROMOTER" -b "$TE" -wa > overlapped_promoter.bed
bedtools subtract -a overlapped_promoter.bed -b "$TE" > overlapped_promoterselves.bed
sort overlapped_promoterselves.bed | uniq > overlapped_promoterselves_uniq.bed #remove those duplicated regions

# Rename duplicated promoter names
Rscript - <<'EOF'
df <- read.table("overlapped_promoterselves_uniq.bed", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
df$V4 <- paste0(df$V4, "_", ave(seq_along(df$V4), df$V4, FUN = seq_along))
write.table(df, file="overlapped_promoterselves_uniq_rename.bed", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
EOF

normalize_score() {
  local infile=$1
  local outfile=$2
  awk 'BEGIN{OFS="\t"}{$5 = ($5 ~ /^-?[0-9]+$/ ? $5 : 0); print}' "$infile" > "$outfile"
}
GENE_NORM="$PWD/gene_norm.bed"
PROMOTER_NORM="$PWD/promoter_norm.bed"
PROMOTERSELF_NORM="$PWD/overlapped_promoterselves_uniq_rename_norm.bed"
TE_NORM="$PWD/TE_norm.bed"
normalize_score "$GENE" "$GENE_NORM"
normalize_score "$PROMOTER" "$PROMOTER_NORM"
normalize_score "overlapped_promoterselves_uniq_rename.bed" "$PROMOTERSELF_NORM"
normalize_score "$TE" "$TE_NORM"

#####################################
# 2. Methylation preprocessing
#####################################
CURRENT_STEP="step 2 - preprocess methylation data"
echo "[STEP 2] Preprocessing methylation data..."
# Keep the chromosome names exactly as they appear in the FAI after upstream
# normalization so wigToBigWig and downstream tables use the same reference set.
awk 'NF >= 2 {print $1"\t"$2}' "$faidx" | awk '!seen[$1]++' > chrom.size

# clean old wig/bw/tab and set work dir for new ones
WIG_DIR="$PWD/wig"
export WIG_DIR
rm -rf "$WIG_DIR"
mkdir -p "$WIG_DIR"
rm -f ./*.wig ./*.bw ./*.tab

for f in "${METH_FILES[@]}"; do
(
    echo "[INFO] Processing $f"
    base=$(basename "$f" .CGmap.gz)  # e.g. FB_01

    # CGmap -> wig
    perl "$CGMAP2WIGGLE" "$f"

    # clean wig: drop any log lines (e.g. "processing chromosomes.....") before conversion
    for ctx in CG CHG CHH; do
        wigfile=${base}.${ctx}.wig
        [[ -f "$wigfile" ]] || continue
        awk '(/^processing / || /^processing chromosomes/){next}
             /^track /{print; next}
             /^variableStep /{print; next}
             ($1 ~ /^[^ \t]+$/ && $2 ~ /^-?[0-9.eE]+$/){print}' "$wigfile" > "${wigfile}.clean" \
          && mv "${wigfile}.clean" "$wigfile"
        mv "$wigfile" "$WIG_DIR/"
    done

    # wig -> bigWig
    #wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig
    #chmod +x wigToBigWig
    "$WIG_TO_BIGWIG" "$WIG_DIR/${base}.CG.wig"  chrom.size "$WIG_DIR/${base}.CG.bw"
    "$WIG_TO_BIGWIG" "$WIG_DIR/${base}.CHG.wig" chrom.size "$WIG_DIR/${base}.CHG.bw"
    "$WIG_TO_BIGWIG" "$WIG_DIR/${base}.CHH.wig" chrom.size "$WIG_DIR/${base}.CHH.bw"

    # bigWigAverageOverBed for annotations
    # wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigAverageOverBed
    # chmod +x bigWigAverageOverBed
    for ctx in CG CHG CHH; do
        "$BIGWIG_AVERAGE" "$WIG_DIR/${base}.${ctx}.bw" "$TE_NORM" "$WIG_DIR/${base}_TE_${ctx}.tab" 2>>"$LOG"
        "$BIGWIG_AVERAGE" "$WIG_DIR/${base}.${ctx}.bw" "$GENE_NORM" "$WIG_DIR/${base}_gene_${ctx}.tab" 2>>"$LOG"
        "$BIGWIG_AVERAGE" "$WIG_DIR/${base}.${ctx}.bw" "$PROMOTER_NORM" "$WIG_DIR/${base}_promoter_${ctx}.tab" 2>>"$LOG"
        "$BIGWIG_AVERAGE" "$WIG_DIR/${base}.${ctx}.bw" "$PROMOTERSELF_NORM" "$WIG_DIR/${base}_promoterselves_${ctx}.tab" 2>>"$LOG"
    done
) &

done


wait

#####################################
# 3. Merge in R (per-stage averages)
#####################################
CURRENT_STEP="step 3 - merge methylation tables"
echo "[STEP 3] Merging methylation tables (stage averages) in R..."

Rscript - "$@" <<'EOF'
wig_dir <- Sys.getenv("WIG_DIR", ".")

process_group <- function(feature, ctx){
  tabs <- list.files(path=wig_dir, pattern=paste0("_",feature,"_",ctx,".tab$"), full.names=TRUE)
  if (length(tabs) == 0) {
    message("No tables for ", feature, " ", ctx, "; skipping.")
    return(NULL)
  }
  lst <- lapply(tabs, function(f){
    df <- read.table(f, header=FALSE, stringsAsFactors=FALSE)
    if (nrow(df) == 0) return(NULL)
    # bigWigAverageOverBed yields: name size covered sum mean (5 cols) or more with extra stats
    val_col <- ncol(df)   # use last column as methylation value
    depth_col <- if (ncol(df) >= 3) 3 else val_col
    fname <- basename(f)
    stage_rep <- sub(paste0("_",feature,"_",ctx,".tab$"),"",fname)
    parts <- strsplit(stage_rep, "_")[[1]]
    stage <- parts[1]
    replicate <- parts[2]
    
    # especially for promoterselves data, merge the methylation values of same promoter
    if(feature == "promoterselves"){
      df$V1 <- sub("_[0-9]+$", "", df$V1)
      df_grouped <- aggregate(. ~ V1, data = df, FUN = sum)
      df_grouped$V6 <- df_grouped$V4 / df_grouped$V3
      df_grouped[is.na(df_grouped)] <- 0
      df <- df_grouped
    }

    data.frame(ID=df[[1]],
               value=df[[val_col]],
               depth=df[[depth_col]],
               stage=stage,
               replicate = replicate,
               stringsAsFactors=F)
  })

  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) {
    message("No non-empty tables for ", feature, " ", ctx, "; skipping.")
    return(NULL)
  }

  df <- do.call(rbind, lst)

  # remove features with depth < 5
  C3 <- aggregate(depth ~ ID, data=df, min)
  keep <- C3$ID[C3$depth >= 5]
  df <- df[df$ID %in% keep, ]

  # calculate the methylation value of each ID at each stage 
  avg_df <- aggregate(value ~ ID + stage, data=df, FUN=function(x) mean(x, na.rm=TRUE))

  # Wide format: ID -> row, stage -> column
  wide <- reshape(avg_df, idvar="ID", timevar="stage", direction="wide")

  # remove "value." in column names
  names(wide) <- sub("^value\\.", "", names(wide))

  # output
  out <- file.path(getwd(), paste0("Tab_",feature,"_",ctx,".txt"))
  write.table(wide, file=out, sep="\t", quote=F, row.names=F)
  message("Wrote: ", out)
}

for(f in c("TE","gene","promoter","promoterselves")){
  for(ctx in c("CG","CHG","CHH")){
    process_group(f,ctx)
  }
}
EOF


#####################################
# Timer end
#####################################
END=$(date +%s)
RUNTIME=$((END-START))
echo "[DONE] All outputs generated."
echo "[TIME] Total runtime: $RUNTIME seconds ($((RUNTIME/60)) min)"

rm overlapped_promoter.bed overlapped_promoterselves.bed overlapped_promoterselves_uniq.bed overlapped_promoterselves_uniq_rename.bed
# consolidate wig/bw/tab files
mkdir -p wig
shopt -s nullglob
files=( *.tab *.bw *.wig )
shopt -u nullglob
if (( ${#files[@]} )); then
  mv "${files[@]}" wig/
fi
