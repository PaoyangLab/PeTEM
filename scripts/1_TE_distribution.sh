#!/usr/bin/env bash
set -euo pipefail

MODULE_NAME="Module 1"
SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
CURRENT_STEP="argument parsing"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

die() {
  echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: $*" >&2
  exit 1
}

usage() {
  echo "Usage: bash 1_TE_distribution.sh -g gene.bed -c CDS.bed -5 5UTR.bed -e exon.bed -3 3UTR.bed -p promoter.bed -t TE.bed -fai genome.fa.fai" >&2
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

while [[ $# -gt 0 ]]; do
  case $1 in
    -g) gene=$2; shift 2;;
    -c) cds=$2; shift 2;;
    -5) utr5=$2; shift 2;;
    -e) exon=$2; shift 2;;
    -3) utr3=$2; shift 2;;
    -p) promoter=$2; shift 2;;
    -t) te=$2; shift 2;;
    -fai) faidx=$2; shift 2;;
    *) die "unknown option: $1";;
  esac
done

missing_args=()
[[ -n ${gene:-} ]] || missing_args+=("-g gene.bed")
[[ -n ${cds:-} ]] || missing_args+=("-c CDS.bed")
[[ -n ${utr5:-} ]] || missing_args+=("-5 5UTR.bed")
[[ -n ${exon:-} ]] || missing_args+=("-e exon.bed")
[[ -n ${utr3:-} ]] || missing_args+=("-3 3UTR.bed")
[[ -n ${promoter:-} ]] || missing_args+=("-p promoter.bed")
[[ -n ${te:-} ]] || missing_args+=("-t TE.bed")
[[ -n ${faidx:-} ]] || missing_args+=("-fai genome.fa.fai")
(( ${#missing_args[@]} == 0 )) || usage_missing "${missing_args[*]}"

CURRENT_STEP="input validation"
require_file "-g" "$gene"
require_file "-c" "$cds"
require_file "-5" "$utr5"
require_file "-e" "$exon"
require_file "-3" "$utr3"
require_file "-p" "$promoter"
require_file "-t" "$te"
require_file "-fai" "$faidx"

#------------------
# Sorting
CURRENT_STEP="step 1 - sorting and interval preparation"
sort -k1,1 -k2,2n "$gene" > gene_sort.bed
sort -k1,1 -k2,2n "$cds" > CDS_sort.bed
sort -k1,1 -k2,2n "$exon" > exon_sort.bed
sort -k1,1 -k2,2n "$utr5" > UTR5_sort.bed
sort -k1,1 -k2,2n "$utr3" > UTR3_sort.bed
sort -k1,1 -k2,2n "$promoter" > promoter_sort.bed
sort -k1,1 -k2,2n "$te" > TE_sort.bed
# Intron
bedtools subtract -a gene_sort.bed -b exon_sort.bed -s > intron.bed
sort -k1,1 -k2,2n -k3,3n intron.bed > intron_sort.bed
# IGR
bedtools merge -i gene_sort.bed > merge_2strands_gene.bed

awk '{print $1"\t0\t"$2}' "$faidx" > genome.bed

bedtools subtract -a genome.bed -b merge_2strands_gene.bed > IGR.bed
sort -k1,1 -k2,2n -k3,3n IGR.bed > IGR_sort.bed

# Merge with no strand info
bedtools merge -i CDS_sort.bed > merge_2strands_CDS.bed
bedtools merge -i exon_sort.bed > merge_2strands_exon.bed
bedtools merge -i UTR5_sort.bed > merge_2strands_5UTR.bed
bedtools merge -i UTR3_sort.bed > merge_2strands_3UTR.bed
bedtools merge -i promoter_sort.bed > merge_2strands_promoter.bed
bedtools merge -i TE_sort.bed > merge_2strands_TE.bed

bedtools merge -i intron_sort.bed > merge_2strands_intron.bed
bedtools merge -i IGR_sort.bed > merge_2strands_IGR.bed

# Intersect TE with regions
for region in gene CDS 5UTR exon intron 3UTR promoter IGR; do
    bedtools intersect -a merge_2strands_TE.bed -b merge_2strands_${region}.bed > TE_on_${region}.bed
done

#------------------
# Call R script
CURRENT_STEP="step 2 - TE distribution plotting"
Rscript "$SCRIPT_DIR/1_TE_distribution.R" \
  TE_on_gene.bed TE_on_CDS.bed TE_on_5UTR.bed TE_on_exon.bed TE_on_intron.bed TE_on_3UTR.bed TE_on_promoter.bed TE_on_IGR.bed \
  merge_2strands_TE.bed merge_2strands_gene.bed merge_2strands_CDS.bed merge_2strands_5UTR.bed merge_2strands_exon.bed \
  merge_2strands_intron.bed merge_2strands_3UTR.bed merge_2strands_promoter.bed merge_2strands_IGR.bed genome.bed 

CURRENT_STEP="cleanup"
rm -f *_sort.bed
rm -f merge_2strands_*bed
rm -f TE_on_*bed
