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

## System requirements
### R environment:
*	R version ≥ 4.2 (tested on 4.3.2)
*	Required R packages:
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

### Python environment:
*	Python ≥ 3.8 (tested on 3.8.10)
*	Required Python packages:
    * pandas (≥ 1.2.4)
    * (uses built-in: glob, os, time)

### Bioinformatics tools:
*	samtools (tested on 1.10)
*	bedtools (tested on v2.27.1)
*	wigToBigWig
*	bigWigAverageOverBed

## Input Files
Depending on which steps you choose to run, you need some or all of the following files:

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
* Some genome annotations are commonly used and can be downloaded via the [provided links](https://paoyang.ipmb.sinica.edu.tw/Termitomyces/reference_genome/). The available species are:
    * Animals: Human, mouse, zebrafish, fly
    * Plants: Arabidopsis, rice, maize, soybean
    * Fungi: Truffle, Botrytis cinerea, Magnaporthe oryzae; Basidiomycetes: Pleurotus ostreatus, Puccinia graminis f. sp. Tritici

### Expression Data
* DEG.txt – Differentially expressed genes (Step 0, 3-2, 4, 5)
* DETE.txt – Differentially expressed TEs (Step 0, 3-2, 4, 5)
> In the DEG/DETE files, the row names are the gene and TE names, followed by columns showing average expression level (RPKM) of each conditions. The rest of columns shows log2 fold change, p value, and FDR comparing each two conditions.
> The column names should be: conditions names, "logFC_condition1_condition2", "PValue_condition1_condition2", "FDR_condition1_condition2", "logFC_condition2_condition3", ... etc.
```
             WT      drdd    logFC_drdd_WT   PValue_drdd_WT  FDR_drdd_WT
AT1G01010    1.58    2.39    0.61            0.21            1
AT1G01020    5.60    5.43    -0.04           0.87            1
AT1G01030    1.49    3.39    1.17            0.01            0.12
```  
These files will automatically be converted into expression matrices (gene_expression.txt, TE_expression.txt) for steps 0, 3-2, and 4.

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

### Parameter: Optionally include unexpressed TEs
* Default: `n` (do not include unexpressed TEs)
* This parameter determines whether to include TEs with zero expression across all samples in the analysis.
* Setting this to `y` includes all TEs regardless of expression level, which may increase the number of analyzed TEs but also add background noise.
* Recommended to keep `n` unless specifically investigating silent or lowly expressed TEs.

### Step 0. Preprocessing
* Generate promoter regions (promoter.bed)
* Integrate methylation and expression data
* __Inputs:__
    * `gene.bed`, `TE.bed`, `genome.fa.fai`, `DEG.txt`, `DETE.txt`, `*.CGmap.gz`
* __Parameters:__ 
    * __Promoter region:__ The default promoter is defined as `1500` bp upstream to `500` bp downstream from the transcription start site (TSS). Users can customize this range by entering other upstream/downstream length from TSS. 

### Step 1. TE Distribution
* Analyze TE distribution across genomic features
* __Inputs:__
    * `gene.bed`, `CDS.bed`, `UTR5.bed`, `exon.bed`, `UTR3.bed`, `TE.bed`, `genome.fa.fai`
    * `promoter.bed`: generated automatically in Step 0

### Step 2. Promoter-embedded TE Families
* Identify enriched TE families overlapping with promoters
* __Inputs:__
    * `TE.bed`, `promoter.bed` (from Step 0), `TE_family.txt`

### Step 3. TE Impact Distance
* 3-1 Preprocessing: Prepare methylation files required in Step 3-2
* 3-2 Plotting: Visualize distance impact of TE methylation on gene expression
* __Inputs:__
    * `gene.bed`, `TE.bed`, `DEG.txt` and `DETE.txt` (will be converted to `gene_expression.txt` and `TE_expression.txt`)
* __Parameters:__ 
    * __Limit range:__ The total up- and downstream distance (in bp) to consider for TE–gene impact analysis. (Default `15000` means ±15 kb around genes will be analyzed)
    * __Window size:__ Sliding window size (bp) used to smooth the TE methylation level curve. (Default: `200` )
    * __Tick size:__ The spacing (bp) between x-axis ticks in the resulting plot. Recommended: approximately 1/3 to 1/4 of the limit range. (Default: `5000` )

### Step 4. Correlation (Single Condition)
* Correlate gene expression with TE/promoter methylation and TE expression with TE methylation
* __Inputs:__
    * `DEG.txt` and `DETE.txt` (will be converted to `gene_expression.txt` and `TE_expression.txt`)
* __Parameters:__ 
    * __Window number:__ Number of sliding windows used to smooth the correlation curves. (Default: `156`).
    * __Y-axis limits:__ Controls the maximum value shown in the y-axis of each correlation plot:
        * ylim_CG: gene expression vs TE/promoter CG methylation (Default: `50`)
        * ylim_CHG: gene expression vs TE/promoter CHG methylation (Default: `10`)
        * ylim_CHH: gene expression vs TE/promoter CHH methylation (Default: `10`)
        * ylim_TEexpTEmC_CH: TE expression vs TE CHG/CHH methylation (Default: `15`)
        * ylim_TEexpTEmC_CG: TE expression vs TE CG methylation (Default: `30`)

### Step 5. Correlation (Across Conditions)
* Examine the correlations between changes in TE methylation, TE expression, and gene expression across different conditions
* __Inputs:__
    * `DEG.txt`, `DETE.txt`

## Usage
Run the interactive pipeline:
```
bash run_PeTEM.sh
```

```
Select steps to run (y/n):
0. Preprocessing? (y/n): y
1. TE distribution? (y/n): y
2. Promoter-embedded TE families? (y/n): y
3-1. TE impact distance: preprocessing? (y/n): y
3-2. TE impact distance: plot? (y/n): y
4. Correlation single condition? (y/n): y
5. Correlation across conditions? (y/n): y
```
> According to the selected steps, users need to give the input names or parameters.
```
Gene BED file: gene.bed
TE BED file: TE.bed
Genome fasta index: genome.fa.fai
DEG file: DEG.txt
DETE file: DETE.txt
Methylation CGmap.gz files (space separated): WT_01.CGmap.gz WT_02.CGmap.gz...
Include unexpressed TEs? (y/n, default n): n
Promoter upstream length from TSS (default 1500): 1500
Promoter downstream length from TSS (default 500): 500
TE family file: TE_family.txt
Limit up-/down-stream range (bp)(e.g. 15000): 15000
Tick size (bp)(e.g. 5000): 5000
Window size (bp)(e.g. 200): 200
Window number (default 156):  156
y-axis limit for gene expression vs TE/promoter mC plot (CG, default 50):  50
y-axis limit for gene expression vs TE/promoter mC plot (CHG, default 10):  10
y-axis limit for gene expression vs TE/promoter mC plot (CHH, default 10):  10
y-axis limit for TE expression vs TE mC plot (CH, default 15):  15
y-axis limit for TE expression vs TE mC plot (CG, default 30):  30
```

