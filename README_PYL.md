# Promoter-embedded TE Methylation (PeTEM) analyzer
PeTEM is designed to analyse the impact of TE methylation on neighbouring genes. It integrates genome annotations with methylome and transcriptome data, enabling users to study TE distribution and measure the correlation of promoter-embedded TE methylation with gene expression in different conditions.   


## Pipeline
<img width="962" height="751" alt="image" src="https://github.com/user-attachments/assets/445236f0-74e0-43f2-b58d-4004144b601a" />

### Tutorial
Please follow the [tutorial](https://github.com/yc811/PeTEM/blob/main/Tutorial.md) of example use case.

## Installation

Clone the repository:

```bash
git clone https://github.com/yc811/PeTEM.git
cd PeTEM
```

### Docker image

We provide a container image that bundles the complete runtime used by the pipeline:

```bash
# from repo root
docker build -t petem -f tools/PeTEM/Dockerfile --build-arg PETEM_SRC=tools/PeTEM .

# or from inside tools/PeTEM
docker build -t petem -f Dockerfile .

# validate the image environment
docker run --rm --entrypoint bash petem /opt/petem/scripts/env_check.sh

# launch the interactive runner
docker run --rm -it -v /absolute/path/to/your/data:/data petem
```

Inside the container the runner starts in `/data`. Mount any directory containing your BED, CGmap, and expression files to that path before launching the container.

### Conda setup

Use the checked-in environment definition:

```bash
cd /work4/home/peiyu/tools/PeTEM
conda env create -f environment.yml
conda activate petem
bash scripts/env_check.sh
```

### Local setup

Run the helper script to install required system tools, Python packages, and R libraries (Ubuntu/Debian):

```bash
bash scripts/setup.sh
```

> The script uses `apt-get`, `pip3 --user`, and `Rscript` to install dependencies, then runs `bash scripts/env_check.sh`. If `apt-get` is not available the script prints the package list to install manually.

## System requirements
### Runtime dependencies
* Bash
* Perl
* gzip / gunzip
* awk / sort / uniq

### R environment
* R version ≥ 4.2 (tested on 4.3.2)
* Required R packages:
    * optparse
    * dplyr
    * tidyr
    * zoo
    * reshape2
    * stringr
    * ggplot2
    * gplots
    * ggalluvial
    * ggpointdensity
    * RColorBrewer
    * viridis
    * rlang
    * ggbreak

### Python environment
* Python ≥ 3.8 (tested on 3.8.10+)
* Required Python packages:
    * pandas (≥ 1.2.4)

### Bioinformatics tools and bundled binaries
* samtools
* bedtools
* `wigToBigWig` (bundled in this repository)
* `bigWigAverageOverBed` (bundled in this repository)

### Environment check before running

Before starting the pipeline, validate the runtime:

```bash
cd /work4/home/peiyu/tools/PeTEM
bash scripts/env_check.sh
```

`env_check.sh` verifies:
* `bash`, `python3`, `Rscript`, `bedtools`, `samtools`, `perl`, `awk`, `sort`, `uniq`, `gunzip`
* bundled executables `wigToBigWig`, `bigWigAverageOverBed`, and script `CGmap2Wiggle.pl`
* Python package `pandas`
* all required R packages listed above

## Input Files
Depending on which steps you choose to run, you need some or all of the following files:

### Annotation conversion helpers
If your annotation source is still in GFF3/GTF-like format, generate the BED inputs first:

```bash
Rscript scripts/gff_to_bed.R --gff annotation.gff3 --outdir annotation_bed
Rscript scripts/te_to_bed_family.R --input te_annotation.gff3 --outdir annotation_bed
```

`gff_to_bed.R` writes `gene.bed`, `CDS.bed`, `exon.bed`, `UTR5.bed`, and `UTR3.bed`. Gene rows are filtered from the feature name in column 3 and can optionally be restricted to `protein_coding` entries detected in column 9. Gene names are parsed from column 9 and support formats such as `ID=name;`, `gene_id=name;`, and `gene_id "name";`.

`te_to_bed_family.R` splits one TE annotation source into `TE.bed` and `TE_family.txt`.

### Genome / Annotation
* gene.bed – Gene coordinates (BED)
* TE.bed – Transposable element coordinates (BED)
* CDS.bed – Coding sequence coordinates (Step 1 only)
* UTR5.bed – 5′ UTR coordinates (Step 1 only)
* exon.bed – Exon coordinates (Step 1 only)
* UTR3.bed – 3′ UTR coordinates (Step 1 only)
> BED file format includes 6 columns: chromosome, start, end, name, score, strand
```
Chr1    3631    5899    AT1G01010       0       +
Chr1    6788    9130    AT1G01020       0       -
Chr1    11649   13714   AT1G01030       0       -
Chr1    23121   31227   AT1G01040       0       +
```

* genome.fa.fai – FASTA index
> The fai index file is generated from genome.fasta file, including 5 columns: name, length, offset, linebases, linewidth
```
Chr1    30427671    74              79      80
Chr2    19698289    30812981        79      80
Chr3    23459830    50760691        79      80
Chr4    18585056    74517556        79      80
Chr5    26975502    93337941        79      80
ChrC    154478      120654981       79      80
ChrM    367808      120811562       70      71
```

* TE_family.txt – TE family annotation (Step 2 only)
> TE family annotation includes 2 columns: names of each TE and their family
```
AT1TE52125      LTR/Gypsy
AT1TE42735      LTR/Copia
AT1TE36140      LTR/Copia
AT1TE21850      RC/Helitron
AT1TE95105      RC/Helitron
```

### Expression Data
* DEG.txt – Differentially expressed genes (Step 0, 3, 4, 5)
* DETE.txt – Differentially expressed TEs (Step 0, 3, 4, 5)
> In the DEG/DETE files, the row names are the gene and TE names, followed by columns showing average expression level (RPKM) of each conditions. The rest of columns shows log2 fold change, p value, and FDR comparing each two conditions.
> The column names should be: conditions names, "logFC_condition1_condition2", "PValue_condition1_condition2", "FDR_condition1_condition2", "logFC_condition2_condition3", ... etc.
```
             WT      drdd    logFC_drdd_WT   PValue_drdd_WT  FDR_drdd_WT
AT1G01010    1.58    2.39    0.61            0.21            1
AT1G01020    5.60    5.43    -0.04           0.87            1
AT1G01030    1.49    3.39    1.17            0.01            0.12
```  
These files will automatically be converted into expression matrices (gene_expression.txt, TE_expression.txt) for steps 0, 3, and 4.

### Methylation Data
* *.CGmap.gz files – Per-sample methylation CGmap files
> CGmap files includes 8 columns: chromosome, C or G (forward or reverse strand), position, context (CG/CHG/CHH), dinucleotide, methylation level (0-1), # of reads supporting methylation, depth
```
Chr3    C    556    CG     CG    0.877551    43    49
Chr3    G    557    CG     CG    0.787879    26    33
Chr3    G    558    CHG    CC    0.405405    15    37
Chr3    G    560    CHH    CA    0.102564    4     39
```  


## Pipeline Steps
Upon running run_pipeline.sh, you will be asked which steps to execute (y/n). You will also provide all required files and parameters upfront. Steps are modular:

### Step 0. Preprocessing
* Generate promoter regions (promoter.bed)
* Integrate methylation and expression data
* __Inputs:__
    * gene.bed, TE.bed, genome.fa.fai, DEG.txt, DETE.txt, *.CGmap.gz
    * Promoter upstream/downstream length (default: 1500 / 500)

### Step 1. TE Distribution
* Analyze TE distribution across genomic features
* __Inputs:__
    * gene.bed, CDS.bed, UTR5.bed, exon.bed, UTR3.bed, TE.bed, genome.fa.fai
    * promoter.bed: generated automatically in Step 0

### Step 2. Promoter-embedded TE Families
* Identify enriched TE families overlapping with promoters
* __Inputs:__
    * TE.bed, promoter.bed (from Step 0), TE_family.txt

### Step 3. TE Impact Distance
* Module 3 now runs as a single entry point.
* If `pre_step3` is incomplete or missing required stage/context files, PeTEM runs preprocessing first.
* If `pre_step3` already contains the required `pre3_<stage>_{CG,CHG,CHH}.bed` files, PeTEM skips preprocessing and runs plotting directly.
* __Inputs:__
    * gene.bed, TE.bed, DEG.txt and DETE.txt (will be converted to gene_expression.txt and TE_expression.txt)
    * Methylation `*.CGmap.gz` files are only required when preprocessing must be run
    * Parameters: limit range, tick size, window size
    * Optionally include unexpressed TEs (y/n, default n)

### Step 4. Correlation (Single Condition)
* Correlate gene expression with TE/promoter methylation and TE expression with TE methylation
* __Inputs:__
    * DEG.txt and DETE.txt (will be converted to gene_expression.txt and TE_expression.txt)
    * Exactly two stages / conditions shared across expression and methylation inputs
    * Parameters: window number, y-axis limits (gene exp vs TE mC:CG, CHG, CHH; TE exp vs TE mC: CG, CH)

### Step 5. Correlation (Across Conditions)
* Examine the correlations between changes in TE methylation, TE expression, and gene expression across different conditions
* __Inputs:__
    * DEG.txt, DETE.txt
    * Exactly two stages / conditions and exactly one pairwise comparison shared across DEG, DETE, and methylation inputs
    * Optionally include unexpressed TEs (y/n, default n)

## Usage
Run the interactive pipeline:
```
bash scripts/run_PeTEM.sh
```

`scripts/run_PeTEM.sh` now performs `scripts/env_check.sh` automatically before any pipeline step starts.

```
Select steps to run (y/n):
0. Preprocessing? (Y/n): y
1. TE distribution? (Y/n): y
2. Promoter-embedded TE families? (Y/n): y
3. TE impact distance? (Y/n): y
4. Correlation single condition? (Y/n): y
5. Correlation across conditions? (Y/n): y
Data working directory (default current directory): /analysis/petem_run
```
> According to the selected steps, you will be prompted for the required files and parameters. Press **Enter** to accept the default shown in parentheses, or type a different value. Invalid paths will be rejected and you will be asked again.
```
Gene BED file: /analysis/petem_run/gene.bed
TE BED file: /analysis/petem_run/TE.bed
Genome fasta index (.fai): /analysis/petem_run/TAIR10.fa.fai
DEG file: /analysis/petem_run/DEG.txt
DETE file: /analysis/petem_run/DETE.txt
Methylation CGmap.gz files (space separated): /analysis/petem_run/WT_01.CGmap.gz /analysis/petem_run/WT_02.CGmap.gz
Include unexpressed TEs? (y/N): n
Promoter upstream length from TSS (bp, default 1500): 1500
Promoter downstream length from TSS (bp, default 500): 500
TE family file: /analysis/petem_run/TE_family.txt
Limit up-/down-stream range (bp, default 15000): 15000
Tick size (bp, default 5000): 5000
Window size (bp, default 200): 200
Window number (default 156): 156
y-axis limit gene exp vs TE/promoter mC (CG, default 50): 50
y-axis limit gene exp vs TE/promoter mC (CHG, default 10): 10
y-axis limit gene exp vs TE/promoter mC (CHH, default 10): 10
y-axis limit TE exp vs TE mC (CH, default 15): 15
y-axis limit TE exp vs TE mC (CG, default 30): 30
```

All intermediate files and final results are written inside the selected data working directory.
