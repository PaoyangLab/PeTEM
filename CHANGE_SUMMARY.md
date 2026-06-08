# PeTEM Change Summary

This document summarizes the major updates currently applied to this PeTEM workspace.

## 1. Reproducible Environment

- Added `Dockerfile`
- Added `environment.yml`
- Added `setup.sh`
- Added `env_check.sh`
- Added an environment preflight check before pipeline execution
- Explicit dependency validation now covers:
  - `bash`
  - `python3`
  - `Rscript`
  - `bedtools`
  - `samtools`
  - `perl`
  - `gzip/gunzip`
  - `awk`, `sort`, `uniq`
  - required R/Python packages

## 2. CLI Entry Point

- Added `petem` CLI wrapper
- Supported commands:
  - `petem --0` / `0_preprocessing`
  - `petem --1` / `1_TE_distribution`
  - `petem --2` / `2_TE_families`
  - `petem --3` / `3_TE_impact_distance`
  - `petem --4` / `4_single_condition_correlation`
  - `petem --5` / `5_cross_condition_correlation`

## 3. Module Naming and Error Handling

- Main user-facing labels changed from `Step 0~5` to `Module 0~5`
- `step` is now only used for internal substeps
- Main shell scripts use `set -euo pipefail`
- Standardized fatal error reporting:
  - module
  - script
  - step
  - missing file / missing argument / exit code
- Missing required arguments now show usage with clearer context

## 4. Path Handling and Output Structure

- Reduced reliance on the current working directory
- Added explicit module-level output handling
- Module 0 derived annotation files are stored under:
  - `OUTPUT_DIR/tmp/module_0_annotation/`
- Module 4 and Module 5 now explicitly read dependencies from `--module0-dir`

## 5. GFF / Annotation Processing

- Added `gff_to_bed.R`
- Supported outputs from GFF:
  - `gene.bed`
  - `CDS.bed`
  - `exon.bed`
  - `UTR5.bed`
  - `UTR3.bed`
- `gene` parsing supports:
  - feature filtering by column 3
  - optional `protein_coding` filtering from column 9
- Gene name parsing supports multiple formats, including:
  - `ID=name;`
  - `gene_id=name;`
  - `gene_id "name";`
- `petem --0` can now start from a GFF file directly and auto-generate:
  - `gene.bed`
  - `promoter.bed`

## 6. TE Input Refactor

- Added `te_to_bed_family.R`
- TE annotation can now be split into:
  - `TE.bed`
  - `TE_family.txt`

## 7. Module 1 Updates

- Confirmed Module 1 reads:
  - `gene.bed`
  - `CDS.bed`
  - `exon.bed`
  - `UTR5.bed`
  - `UTR3.bed`
  - `promoter.bed`
- Genomic feature order changed to:
  - `Promoter`
  - `Gene`
  - `5'UTR`
  - `CDS`
  - `Exon`
  - `Intron`
  - `3'UTR`
  - `IGR`
- Added short CLI usage:
  - `petem --1 -g <genome_bed_dir> -t <TE.bed> -f <fai> -o <output_dir>`

## 8. Module 2 Updates

- Enrichment is now consistently `log2`
- `TE Types` follow `TE_family.txt` directly
- Plot styling adjusted
- Added `petem --2`

## 9. Module 3 Updates

- Added `petem --3`
- Module 3 can run from:
  - `gene.bed`
  - or a directory containing `gene.bed`
- Auto skip logic:
  - if `pre_step3` is complete, preprocessing is skipped
- Added optional control plotting:
  - `-c y`
  - normal plots use `wTE`
  - control plots use `woTE`
- Added new options:
  - `-d` distance limit, default `15000`
  - `-p` top/bottom expression percentile, default `10`
  - `-w` window size, default `100`
  - `-l poly` polynomial smoothing with CI
- Current default line behavior:
  - moving average
  - `100 bp` window
  - `100 bp` stride
- Dot transparency reduced

## 10. Module 4 Updates

- Added `petem --4`
- Module 4 now reads Module 0 dependencies from `--module0-dir`
- Correlation bar plots updated:
  - custom CG / CHG / CHH colors
  - TE shown in gray
  - grouped y-axis labels
  - no y-axis ticks
  - TE separator lines spanning the panel
  - uppercase x-axis labels:
    - `Pearson correlation coefficient`
    - `Spearman correlation coefficient`
- `TEexp_TEmC_line` plots updated:
  - legend moved above the panel
  - left y-axis label:
    - `TE CG methylation level (%)`
  - right y-axis label:
    - `TE CH methylation level (%)`
- `geneexp_TEm*` plots updated:
  - captions moved outside the plotting panel
  - `Promoter with TEs` / `Promoter without TEs` wording standardized
  - per-line `rho` labels added near the line starts

## 11. Module 5 Updates

- Added `petem --5`
- Module 5 requires:
  - `--DEG`
  - `--DETE`
  - `--module0-dir`
- Scatter plots changed from point-density style to regular points

## 12. Input Constraints

- Module 4 requires exactly two shared conditions/stages
- Module 5 requires exactly one cross-condition comparison between the same two stages
- DEG / DETE / methylation inputs must be stage-consistent

## 13. Language and Script Cleanup

- Removed Chinese comments/messages from code paths that were updated
- Standardized plotting defaults through `plot_defaults.R`
- Folder creation uses `mkdir -p`

## 14. Generated Example Outputs

Example test outputs were produced during development under:

- `test/module_1`
- `test/module_2`
- `test/module_3`
- `test/module_3_control`
- `test/module_4`
- `test/module_5`

These are development/test outputs and are not required as source inputs.
