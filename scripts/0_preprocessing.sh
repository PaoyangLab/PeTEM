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
    echo "Usage: bash 0_preprocessing.sh -g gene.bed -t TE.bed -p promoter.bed -eg expression_gene.txt -et expression_TE.txt -fai genome.fa.fai -m sample1.CGmap.gz [sample2.CGmap.gz ...]" >&2
    exit 1
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
(( ${#METH_FILES[@]} > 0 )) || missing_args+=("-m sample1.CGmap.gz [sample2.CGmap.gz ...]")
(( ${#missing_args[@]} == 0 )) || usage_missing "${missing_args[*]}"

CURRENT_STEP="input validation"
require_file "-g" "$GENE"
require_file "-t" "$TE"
require_file "-p" "$PROMOTER"
require_file "-eg" "$GENE_EXP"
require_file "-et" "$TE_EXP"
require_file "-fai" "$faidx"
for f in "${METH_FILES[@]}"; do
    require_file "-m" "$f"
done

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
