# PeTEM Module 0, 2, and 5 CLI

Run these commands from the PeTEM repository root. Replace example filenames with your own paths.

## Module 0: prepare two conditions

```bash
./petem --0 \
  -g annotation.gff3 -t TE_with_family.txt \
  -eg gene_expression.txt -et TE_expression.txt \
  -f genome.fa.fai \
  -m1 stage1_rep1.CGmap.gz stage1_rep2.CGmap.gz \
  -m2 stage2_rep1.CGmap.gz stage2_rep2.CGmap.gz \
  -o module0
```

- `-m1` and `-m2` correspond to the first and second stage columns in the gene and TE expression tables. For more stages, use `-m3` through `-m5` as needed. Each stage needs at least one CGmap file.
- The original `-m file1.CGmap.gz ...` mode is still supported. Do not mix `-m` with `-m1` through `-m5`.
- Optional `-DEG DEG.txt -DETE DETE.txt` supplies both differential-expression tables instead of generating them; provide the two options together.
- Module 0 writes `PETEM_MODULE0_MANIFEST.json` and `LOG_0_input_mapping.log` in its output directory. The manifest carries prepared TE and gene BED paths for modules 2 and 5.

## Module 2: TE features

```bash
./petem --2 \
  --manifest module0/PETEM_MODULE0_MANIFEST.json \
  --feature y --fasta genome.fa -o module2
```

- `--feature` defaults to `n`. With `y`, module 2 makes separate GC content, TE length, and nearest TE-edge-to-TSS distance bar plots for all TEs and promoter-embedded TEs. Error bars show SEM.
- `--fasta` must have a matching `genome.fa.fai`; `bedtools` is required. The gene BED comes from the manifest, or can be set with `--gene-bed`.
- Without a manifest, provide `-t TE.bed -T TE_family.txt -i TE_overlap_promoter.bed` and, for `--feature y`, `--gene-bed gene.bed`.

## Module 5: change filters and TE features

```bash
./petem --5 \
  --manifest module0/PETEM_MODULE0_MANIFEST.json \
  -expchange pvalue 0.05 -mcchange 0 \
  --feature y --fasta genome.fa -o module5
```

- `-expchange pvalue 0.05` is the default expression-significance filter. Use `-expchange fdr 0.05` to filter by FDR instead; the cutoff must be between 0 and 1.
- `-mcchange 0` is the default minimum absolute TE methylation change in percentage points. Use a larger non-negative value to tighten the filter.
- `--feature` defaults to `n`. With `y`, module 5 makes separate GC content, TE length, and TE-TSS distance bar plots, plus a four-bar TE-family composition plot. Each plot compares all TEs, promoter-embedded TEs, Positive TEs, and Negative TEs. Feature error bars show SEM.
- Positive means Q1/Q3 and Negative means Q2/Q4 in the gene-expression versus TE-methylation change plots, after the selected filters. TEs are counted once per group; the two groups can overlap across methylation contexts or gene pairs. `--positive y` additionally enables the existing Q1/Q3 boxplots, but is not required for the feature plots.
- `--fasta` and its `.fai` index are required with `--feature y`. The manifest supplies TE BED, gene BED, and the original TE-family table. Override them with `--te-bed`, `--gene-bed`, and `--te-family` when needed. The family table must have either two columns (TE ID, family) or at least seven columns (TE ID in column 1, family in column 7).
- Without `--manifest`, use `--module0-dir module0`. Module 5 looks for the prepared BEDs and local manifest there; if no family table is available from the local manifest, provide `--te-family` explicitly.
