#!/usr/bin/env bash
# Download common genomes (FA, GTF, GFF3) for popular species.
# Usage:
#   bash database.sh [all|human|mouse|fly|arabidopsis ...] [-o OUTDIR] [--keep-gz]
# Defaults: download all species into ./genomes and gunzip the files.

set -euo pipefail

OUTDIR="genomes"
KEEP_GZ=0
declare -a TARGETS=()

usage() {
  echo "Usage: $0 [all|human|mouse|fly|arabidopsis ...] [-o OUTDIR] [--keep-gz]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    all)
      TARGETS=(human mouse fly arabidopsis zebrafish maize rice soybean magnaporthe_oryzae truffle celegans)
      shift
      ;;
    human|mouse|fly|arabidopsis|zebrafish|maize|rice|soybean|soybeans|magnaporthe_oryzae|truffle|celegans|c_elegans)
      # normalize aliases
      case "$1" in
        soybeans) TARGETS+=("soybean") ;;
        c_elegans) TARGETS+=("celegans") ;;
        *) TARGETS+=("$1") ;;
      esac
      shift
      ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    --keep-gz) KEEP_GZ=1; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Default to all if none specified
if [[ ${#TARGETS[@]} -eq 0 ]]; then
  TARGETS=(human mouse fly arabidopsis zebrafish maize rice soybean magnaporthe_oryzae truffle celegans)
fi

mkdir -p "$OUTDIR"

require() {
  command -v "$1" >/dev/null 2>&1 || { echo "Error: $1 not found in PATH."; exit 1; }
}
require curl
require gunzip

download_and_unpack() {
  local url="$1"
  local dest="$2"
  echo "[INFO] Downloading ${url}"
  curl -L --fail --retry 3 --retry-delay 5 -o "$dest" "$url"
  if [[ $KEEP_GZ -eq 0 && $dest == *.gz ]]; then
    gunzip -f "$dest"
  fi
}

# Ensembl release numbers are fixed for reproducibility.
ENSEMBL_REL=112
ENSEMBL_PLANTS_REL=60
ENSEMBL_FUNGI_REL=60

declare -A URLS

# Human (GRCh38)
URLS[human_fa]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
URLS[human_gtf]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_REL}.gtf.gz"
URLS[human_gff3]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_REL}.gff3.gz"

# Mouse (GRCm39)
URLS[mouse_fa]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
URLS[mouse_gtf]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_REL}.gtf.gz"
URLS[mouse_gff3]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gff3/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_REL}.gff3.gz"

# Fly (Drosophila melanogaster BDGP6.46)
URLS[fly_fa]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"
URLS[fly_gtf]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.${ENSEMBL_REL}.gtf.gz"
URLS[fly_gff3]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.${ENSEMBL_REL}.gff3.gz"

# Zebrafish (GRCz11)
URLS[zebrafish_fa]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
URLS[zebrafish_gtf]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gtf/danio_rerio/Danio_rerio.GRCz11.${ENSEMBL_REL}.gtf.gz"
URLS[zebrafish_gff3]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gff3/danio_rerio/Danio_rerio.GRCz11.${ENSEMBL_REL}.gff3.gz"

# C. elegans (WBcel235)
URLS[celegans_fa]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
URLS[celegans_gtf]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.${ENSEMBL_REL}.gtf.gz"
URLS[celegans_gff3]="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.${ENSEMBL_REL}.gff3.gz"

# Arabidopsis (TAIR10, EnsemblPlants release 60)
URLS[arabidopsis_fa]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
URLS[arabidopsis_gtf]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${ENSEMBL_PLANTS_REL}.gtf.gz"
URLS[arabidopsis_gff3]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${ENSEMBL_PLANTS_REL}.gff3.gz"

# Maize (Zea mays, B73 RefGen v5)
URLS[maize_fa]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v5.dna.toplevel.fa.gz"
URLS[maize_gtf]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gtf/zea_mays/Zea_mays.B73_RefGen_v5.${ENSEMBL_PLANTS_REL}.gtf.gz"
URLS[maize_gff3]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gff3/zea_mays/Zea_mays.B73_RefGen_v5.${ENSEMBL_PLANTS_REL}.gff3.gz"

# Rice (Oryza sativa Japonica IRGSP-1.0)
URLS[rice_fa]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
URLS[rice_gtf]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.${ENSEMBL_PLANTS_REL}.gtf.gz"
URLS[rice_gff3]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.${ENSEMBL_PLANTS_REL}.gff3.gz"

# Soybean (Glycine max v4.0)
URLS[soybean_fa]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/fasta/glycine_max/dna/Glycine_max.Glycine_max_v4.0.dna.toplevel.fa.gz"
URLS[soybean_gtf]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gtf/glycine_max/Glycine_max.Glycine_max_v4.0.${ENSEMBL_PLANTS_REL}.gtf.gz"
URLS[soybean_gff3]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${ENSEMBL_PLANTS_REL}/gff3/glycine_max/Glycine_max.Glycine_max_v4.0.${ENSEMBL_PLANTS_REL}.gff3.gz"

# Magnaporthe oryzae (MG8, EnsemblFungi release 60)
URLS[magnaporthe_oryzae_fa]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-${ENSEMBL_FUNGI_REL}/fasta/magnaporthe_oryzae/dna/Magnaporthe_oryzae.MG8.dna.toplevel.fa.gz"
URLS[magnaporthe_oryzae_gtf]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-${ENSEMBL_FUNGI_REL}/gtf/magnaporthe_oryzae/Magnaporthe_oryzae.MG8.${ENSEMBL_FUNGI_REL}.gtf.gz"
URLS[magnaporthe_oryzae_gff3]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-${ENSEMBL_FUNGI_REL}/gff3/magnaporthe_oryzae/Magnaporthe_oryzae.MG8.${ENSEMBL_FUNGI_REL}.gff3.gz"

# Truffle (Tuber melanosporum Mel28, EnsemblFungi release 60)
URLS[truffle_fa]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-${ENSEMBL_FUNGI_REL}/fasta/tuber_melanosporum/dna/Tuber_melanosporum.Mel28.dna.toplevel.fa.gz"
URLS[truffle_gtf]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-${ENSEMBL_FUNGI_REL}/gtf/tuber_melanosporum/Tuber_melanosporum.Mel28.${ENSEMBL_FUNGI_REL}.gtf.gz"
URLS[truffle_gff3]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-${ENSEMBL_FUNGI_REL}/gff3/tuber_melanosporum/Tuber_melanosporum.Mel28.${ENSEMBL_FUNGI_REL}.gff3.gz"

for sp in "${TARGETS[@]}"; do
  echo "[INFO] === ${sp} ==="
  sp_dir="${OUTDIR}/${sp}"
  mkdir -p "$sp_dir"
  for type in fa gtf gff3; do
    key="${sp}_${type}"
    url="${URLS[$key]:-}"
    if [[ -z $url ]]; then
      echo "[WARN] No URL configured for ${sp} ${type}, skipping."
      continue
    fi
    fname=$(basename "$url")
    dest="${sp_dir}/${fname}"
    download_and_unpack "$url" "$dest"
  done
done

echo "[DONE] Downloads finished. Files are under ${OUTDIR}/"
