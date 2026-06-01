#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

R_PACKAGES=(
  optparse
  dplyr
  tidyr
  zoo
  reshape2
  stringr
  ggplot2
  gplots
  ggalluvial
  ggpointdensity
  RColorBrewer
  viridis
  rlang
  ggbreak
)

PYTHON_PACKAGES=(
  pandas
)

fail() {
  echo "[env-check] ERROR: $*" >&2
  exit 1
}

check_command() {
  local cmd=$1
  command -v "$cmd" >/dev/null 2>&1 || fail "Missing required command: $cmd"
  echo "[env-check] OK command: $cmd -> $(command -v "$cmd")"
}

check_executable() {
  local path=$1
  [[ -x "$path" ]] || fail "Missing or non-executable file: $path"
  echo "[env-check] OK executable: $path"
}

check_python_packages() {
  python3 - "${PYTHON_PACKAGES[@]}" <<'PY'
import importlib
import sys

missing = []
for pkg in sys.argv[1:]:
    try:
        importlib.import_module(pkg)
    except Exception:
        missing.append(pkg)

if missing:
    print("[env-check] Missing Python packages: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)

print("[env-check] OK Python packages: " + ", ".join(sys.argv[1:]))
PY
}

check_r_packages() {
  Rscript - "${R_PACKAGES[@]}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
missing <- args[!vapply(args, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("[env-check] Missing R packages: ", paste(missing, collapse = ", "))
  quit(status = 1)
}
message("[env-check] OK R packages: ", paste(args, collapse = ", "))
RS
}

main() {
  echo "[env-check] Starting PeTEM environment validation"
  check_command bash
  check_command python3
  check_command Rscript
  check_command bedtools
  check_command samtools
  check_command perl
  check_command awk
  check_command sort
  check_command uniq
  check_command gunzip

  check_executable "$SCRIPT_DIR/wigToBigWig"
  check_executable "$SCRIPT_DIR/bigWigAverageOverBed"
  [[ -f "$SCRIPT_DIR/CGmap2Wiggle.pl" ]] || fail "Missing file: $SCRIPT_DIR/CGmap2Wiggle.pl"
  echo "[env-check] OK file: $SCRIPT_DIR/CGmap2Wiggle.pl"

  echo "[env-check] bash version: $(bash --version | head -n 1)"
  echo "[env-check] python version: $(python3 --version 2>&1)"
  echo "[env-check] R version: $(Rscript --version 2>&1)"
  echo "[env-check] bedtools version: $(bedtools --version 2>&1)"
  echo "[env-check] samtools version: $(samtools --version | head -n 1)"

  check_python_packages
  check_r_packages
  echo "[env-check] Environment validation passed"
}

main "$@"
