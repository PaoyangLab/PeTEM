#!/usr/bin/env bash

set -euo pipefail

MODULE_NAME="Module 3"
SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
CURRENT_STEP="argument parsing"

die() {
    echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: $*" >&2
    exit 1
}

usage() {
    echo "Usage: bash 3_1_TE_impact_distance_preprocess.sh -m sample1.CGmap.gz [sample2.CGmap.gz ...]" >&2
    exit 1
}

require_file() {
    local path=$1
    [[ -f $path ]] || die "input file not found: $path"
}

trap 'rc=$?; echo "[ERROR] ${MODULE_NAME} | ${SCRIPT_NAME} | ${CURRENT_STEP}: command failed with exit code ${rc}" >&2; exit ${rc}' ERR

METH_FILES=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -m)
            shift
            while [[ $# -gt 0 && ! $1 =~ ^- ]]; do
                METH_FILES+=("$1")
                shift
            done
            ;;
        *)
            die "unknown option: $1"
            ;;
    esac
done

if [[ ${#METH_FILES[@]} -eq 0 ]]; then
    shopt -s nullglob
    defaults=(*.CGmap.gz)
    shopt -u nullglob
    if [[ ${#defaults[@]} -eq 0 ]]; then
        die "no CGmap.gz files were provided with -m and none were found in the current directory"
    fi
    METH_FILES=("${defaults[@]}")
fi

CURRENT_STEP="input validation"
for f in "${METH_FILES[@]}"; do
    require_file "$f"
done

LOG="LOG_3_1_TE_impact_distance_preprocess.log"
start_step=$(date +%s)

CURRENT_STEP="step 1 - unzip and filter CGmap files"
echo "[`date`] Preprocess methylation files: (1) unzip + filter" | tee -a "$LOG"

rm -rf pre_step3
mkdir -p pre_step3

for f in "${METH_FILES[@]}"; do
(
    start=$(date +%s)
    base=$(basename "$f" .CGmap.gz)
    gunzip -c "$f" | awk '$8>=4 {print $1"	"$3"	"$2"	"$4"	"$6}' > "pre_step3/${base}.txt"
    end=$(date +%s)
    echo "[INFO] Preprocessed $f in $((end-start)) sec" | tee -a "$LOG"
)&
done
wait
echo "[INFO] All replicates done." | tee -a "$LOG"

CURRENT_STEP="step 2 - compute per-stage average methylation"
echo "[`date`] Preprocess methylation files: (2) calculate average mC of each site at each stage" | tee -a "$LOG"

stages=($(printf '%s
' "${METH_FILES[@]}" | xargs -n1 basename | cut -d'_' -f1 | sort -u))

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
    df = pd.read_csv(f, sep="	", header=None,
                     names=["chr", "site", "nt", "CNN", "mC"])
    dfs.append(df)

combined = pd.concat(dfs, axis=0, ignore_index=True)
m = combined.groupby(["chr", "site", "nt", "CNN"], as_index=False)["mC"].mean()
m = m.rename(columns={"mC": f"{stage}_mC"})

m["strand"] = m.nt.replace({"C": "+", "G": "-"})
m["name"] = ["site_" + str(i + 1) for i in range(len(m))]

for type in ["CG", "CHG", "CHH"]:
    sub = m[m.CNN == type][["chr", "site", "site", "name", f"{stage}_mC", "strand"]].dropna()
    sub.to_csv(f"pre_step3/pre3_{stage}_{type}.bed", sep="	", index=False, header=False)
EOF
    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec" | tee -a "$LOG"
)&
done
wait

end_step=$(date +%s)
echo "[INFO] All stages done." | tee -a "$LOG"
echo "[`date`] Preprocessed methylation files finished in $((end_step-start_step)) sec" | tee -a "$LOG"
