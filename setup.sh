#!/usr/bin/env bash
# Bootstrap system dependencies, Python packages, and R libraries required by PeTEM.
set -euo pipefail

APT_PACKAGES=(
  bash
  python3
  python3-pip
  python3-venv
  bedtools
  samtools
  perl
  gzip
  gawk
  libxml2-dev
  libcurl4-openssl-dev
  libssl-dev
  libfontconfig1-dev
  libharfbuzz-dev
  libfribidi-dev
  libfreetype6-dev
  libpng-dev
  libtiff5-dev
  libjpeg-dev
  g++
  make
)

PYTHON_PACKAGES=(
  pandas
)

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

install_apt_packages() {
  if command -v apt-get >/dev/null 2>&1; then
    echo "[setup] Installing APT packages: ${APT_PACKAGES[*]}"
    sudo apt-get update
    sudo apt-get install -y --no-install-recommends "${APT_PACKAGES[@]}"
  else
    echo "[setup] Warning: apt-get not found. Please install the following packages manually:" >&2
    printf '  - %s\n' "${APT_PACKAGES[@]}" >&2
  fi
}

install_python_packages() {
  if ! command -v python3 >/dev/null 2>&1; then
    echo "[setup] Error: python3 not found. Install Python 3 before proceeding." >&2
    return 1
  fi
  if command -v pip3 >/dev/null 2>&1; then
    echo "[setup] Installing Python packages: ${PYTHON_PACKAGES[*]}"
    pip3 install --user --no-cache-dir "${PYTHON_PACKAGES[@]}"
  else
    echo "[setup] Warning: pip3 not found. Skipping Python packages." >&2
  fi
}

install_r_packages() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "[setup] Error: Rscript not found. Install R (>=4.2) before proceeding." >&2
    return 1
  fi
  echo "[setup] Installing R packages: ${R_PACKAGES[*]}"
  Rscript - <<'RS'
packages <- c(
  "optparse","dplyr","tidyr","zoo","reshape2","stringr",
  "ggplot2","gplots","ggalluvial","ggpointdensity",
  "RColorBrewer","viridis","rlang","ggbreak"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) == 0) {
  message("All R packages already installed.")
} else {
  repos <- getOption("repos")
  if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  install.packages(missing, dependencies = TRUE)
}
RS
}

main() {
  install_apt_packages
  install_python_packages
  install_r_packages
  if [[ -x "$(dirname "$0")/env_check.sh" ]]; then
    echo "[setup] Running environment validation."
    bash "$(dirname "$0")/env_check.sh"
  fi
  echo "[setup] Completed dependency installation."
}

main "$@"
