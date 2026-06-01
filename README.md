# Promoter-embedded TE Methylation (PeTEM) analyzer
PeTEM is designed to analyse the association between promoter-embedded TE methylation and neighbouring gene expression. It integrates genome annotations with methylome and transcriptome data, enabling users to evaluate the genome-wide correlations between TE methylation and gene expression, and identify specific TE and gene pairs with associated methylation and expression changes across different conditions (such as stages, treatments, or phenotypes).

<img width="865" height="439" alt="image" src="https://github.com/user-attachments/assets/56516953-0937-4c9c-951f-f474a6b0ee67" />

## Table of Contents

## Table of Contents

- [Tutorial](#tutorial)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Input Files](#input-files)
  - [Genome Annotation](#genome-annotation)
  - [Methylation Data](#methylation-data)
  - [Expression Data](#expression-data)
- [Pipeline Modules](#pipeline-modules)
  - [Module 0. Preprocessing](#module-0-preprocessing)
  - [Module 1. TE Distribution](#module-1-te-distribution)
  - [Module 2. Enriched promoter-embedded TE Families](#module-2-enriched-promoter-embedded-te-families)
  - [Module 3. TE Methylation Near Gene](#module-3-te-methylation-near-gene)
  - [Module 4. Correlation Between Methylation and Expression](#module-4-correlation-between-methylation-and-expression)
  - [Module 5. Associated TE and Gene Pairs](#module-5-associated-te-and-gene-pairs)

## Tutorial
Please see the [tutorial](https://github.com/PaoyangLab/PeTEM/blob/main/Tutorial.md) for an example workflow.

## System requirements

### Runtime dependencies 
* Bash
* Perl
* gzip / gunzip
* awk / sort / uniq

### Environment setup
> We provide three alternative methods for setting up the environment [here](#set-up-environment).

<details>
<summary><b><ins>R version ≥ 4.2 (tested on 4.3.2)</ins></b></summary>
  
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
  * ggbreak
  * RColorBrewer
  * viridis
  * rlang

</details>

<details>
<summary><b><ins>Python version ≥ 3.8 (tested on 3.8.10)</ins></b></summary>
  
  * pandas (≥ 1.2.4)

</details>

<details>
<summary><b><ins>Bioinformatics tools</ins></b></summary>
  
  *	samtools (tested on 1.10)
  *	bedtools (tested on v2.27.1)
  *	wigToBigWig (bundled in this repository)
  *	bigWigAverageOverBed (bundled in this repository)

</details>


## Installation

### Clone repository:
```bash
git clone https://github.com/PaoyangLab/PeTEM.git
cd PeTEM
```

### Download example data:
```bash
wget https://github.com/PaoyangLab/PeTEM/releases/download/v0.0.3/PeTEM_data.tar.gz
tar -xzvf PeTEM_data.tar.gz
```

### Set up environment 
#### Docker image
  
> Inside the container the runner starts in `/data`. Mount any directory containing your BED, CGmap, and expression files to that path before launching the container.

```bash
docker build -t petem:local .
docker run --rm petem:local --help
```

#### Conda setup

> Use the checked-in environment definition:

```bash
conda env create -f environment.yml
conda activate petem
bash env_check.sh ##optional
```

#### Local setup

> The script uses apt-get, pip3 --user, and Rscript to install dependencies, then runs bash env_check.sh. If apt-get is not available the script prints the package list to install manually

```bash
bash setup.sh
bash env_check.sh ##optional
```

## Input Files
PeTEM integrates inputs data including [genome annotations](#genome-annotation), [DNA methylation data](#methylation-data), and [expression data](#expression-data). 
In PeTEM, running the module 1 and 2 rely solely on annotation data, while running the remaining modules additionally require methylation and expression data.

### Genome Annotation
#### General features annotations file (GFF3 format) 
> The genome annotation file contains genomic feature coordinates and hierarchical annotations, including genes, transcripts, CDS regions, exons, and untranslated regions (5′UTR and 3′UTR).



`genomic.gff`
| sequence ID | source | feature type | feature start | feature end | score | strand | phase | attributes |
|---|---|---|---|---|---|---|---|---|
| Chr1 | Araport11 | gene | 3631 | 5899 | . | + | . | ID=AT1G01010;Name=AT1G01010;full_name=NAC domain containing protein 1 |
| Chr1 | Araport11 | mRNA | 3631 | 5899 | . | + | . | ID=AT1G01010.1;Name=AT1G01010.1;Parent=AT1G01010 |
| Chr1 | Araport11 | CDS | 3760 | 3913 | . | + | 0 | ID=AT1G01010:CDS:1;Parent=AT1G01010.1 |
| Chr1 | Araport11 | exon | 3631 | 3913 | . | + | . | ID=AT1G01010:exon:1;Parent=AT1G01010.1 |
| Chr1 | Araport11 | five_prime_UTR | 3631 | 3759 | . | + | . | ID=AT1G01010:five_prime_UTR:1 |
| Chr1 | Araport11 | three_prime_UTR | 5631 | 5899 | . | + | . | ID=AT1G01010:three_prime_UTR:1 |

<details>
<summary><b><ins>Sources of commonly used genome annotations</ins></b></summary>
   
* Animals:
   * [Human (GRCh38/hg38)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
   * [Mouse (GRCm39/mm39)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/)
   * [Zebrafish (GRCz11/danRer11)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002035.6/)
   * [Fruit fly (dm6)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4/)
* Plants:
   * [Arabidopsis (Araport11)](https://www.arabidopsis.org/download/list?dir=Genes%2FAraport11_genome_release)
   * [Rice (IRGSP-1.0)](https://rice.uga.edu/download_osa1r7.shtml)
   * [Maize (Zea mays cv. B73, RefGen_v5)](https://www.maizegdb.org/download)
   * [Soybean (Glycine max cv. Williams 82, Glycine_max_v4.0)](http://ncbi.nlm.nih.gov/datasets/genome/GCF_000004515.6/)
* Fungi
   * [12 species](https://urgi.versailles.inra.fr/download/fungi/TEs/)

</details>

#### Gene annotation file (FASTA index) 
> The genome FASTA index file is generated from genome fasta file (usage: `samtools faidx genome.fa`), providing the names and lengths of each chromosome.

`genome.fa.fai`
| name | length | offset | linebases | linewidth |
|---|---|---|---|---|
| Chr1 | 30427671 | 74 | 79 | 80 |
| Chr2 | 19698289 | 30812981 | 79 | 80 |
| Chr3 | 23459830 | 50760691 | 79 | 80 |
| Chr4 | 18585056 | 74517556 | 79 | 80 |
| Chr5 | 26975502 | 93337941 | 79 | 80 |
| ChrC | 154478 | 120654981 | 79 | 80 |
| ChrM | 367808 | 120811562 | 70 | 71 |

#### Transposable element coordinates 
> The transposable element annotation file contains genomic coordinates and classification information for annotated transposable elements (TEs), including TE family and strand orientation.

`TE.txt`
| TE name | chromosome | start | end | score | strand | TE family |
|---|---|---|---|---|---|---|
| AT1TE00010 | Chr1 | 11897 | 11976 | 0 | + | LTR/Copia |
| AT1TE00020 | Chr1 | 16883 | 17009 | 0 | - | RC/Helitron |
| AT1TE00025 | Chr1 | 17024 | 18924 | 0 | + | RC/Helitron |
| AT1TE00030 | Chr1 | 18331 | 18642 | 0 | - | DNA/HAT |


### Methylation Data
> CGmap files includes 8 columns: chromosome, C or G (forward or reverse strand), position, context (CG/CHG/CHH), dinucleotide, methylation level (0-1), # of reads supporting methylation, depth

`*.CGmap.gz files`
| chromosome | nucleotide | site | context | dinucleotide | methylation level | C site | C+T site |
|---|---|---|---|---|---|---|---|
| Chr3 | C | 556 | CG | CG | 0.877551 | 43 | 49 |
| Chr3 | G | 557 | CG | CG | 0.787879 | 26 | 33 |
| Chr3 | G | 558 | CHG | CC | 0.405405 | 15 | 37 |
| Chr3 | G | 560 | CHH | CA | 0.102564 | 4 | 39 |


### Expression Data
> The expression data includes differentially expressed genes `gene_expression.txt` and differentially expressed TEs `TE_expression.txt`.

> In these files, the row names are the gene and TE names, followed by columns showing average expression level (RPKM) of each conditions. The rest of columns shows log2 fold change, p value, and FDR comparing each two conditions.

> The column names should be arranged by: **conditions names** ("root", "leaf"), **logFC_condition1_condition2** ("logFC_root_leaf"), **PValue_condition1_condition2** ("PValue_root_leaf"), **FDR_condition1_condition2** ("FDR_root_leaf"), logFC_condition2_condition3, ... etc.

`gene_expression.txt`
|  | root | leaf | logFC_root_leaf | PValue_root_leaf | FDR_root_leaf |
|---|---|---|---|---|---|
| AT1G01010 | 11.64 | 10.93 | 0.09 | 0.60 | 1 |
| AT1G01020 | 6.39 | 5.95 | 0.10 | 0.62 | 1 |
| AT1G01030 | 0.65 | 1.01 | -0.63 | 0.19 | 0.84 |

`TE_expression.txt`
|  | root | leaf | logFC_root_leaf | PValue_root_leaf | FDR_root_leaf |
|---|---|---|---|---|---|
| AT1G01010 | 386.06 | 240.27 | 0.63 | 0.73 | 1 |
| AT1G01020 | 0 | 0 | 0 | 1 | 1 |
| AT1G01030 | 0 | 0 | 0 | 1 | 1 |

## Pipeline Modules
Users must run module 0 at the first time to preprocess the input files before running module 1 to 5. Module 1 to 5 are independent and not sequential. 

<img width="865" height="638" alt="image" src="https://github.com/user-attachments/assets/6cafee51-1022-412e-b2e1-ab7e3cb1a213" />


### Module 0. Preprocessing
> Generate promoter regions (`promoter.bed`) and integrate methylation and expression data
* __Inputs:__
    * `genomic.gff`, `TE.txt`, `genome.fa.fai`, `gene_expression.txt`, `TE_expression.txt`, `*.CGmap.gz`
* __Outputs:__
    * `OUTPUT_0_embedded_TE_gene_number.txt`, `promoter.bed`, `TE_overlap_promoter.bed`, `Tab_*.txt`
* __Usage:__
  ```bash
  ./petem --0 \
    -g <ANNOTATION_GFF3> \
    -t <TE_BED> \
    -eg <GENE_EXPRESSION> \
    -et <TE_EXPRESSION> \
    -f <GENOME_FAI> \
    -m <CGMAP_FILES> \
    -o <OUTPUT_DIR>
  ```
* __Parameters:__
  
  | Parameter | Description |
  |---|---|
  | `--0` | Run module 0. |
  | `-g` | Path to the genome annotation file in GFF3 format. |
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-eg` | Path to the gene expression table. |
  | `-et` | Path to the transposable element expression table. |
  | `-f` | Path to the genome FASTA index file (`.fai`). This is used to automatically generate gene, promoter, intergenic region, and other BED files. |
  | `-m` | One or more CGmap files, separated by spaces. Compressed files such as `.CGmap.gz` are supported. |
  | `-o` | Output directory for module 0 results. |




### Module 1. TE Distribution
> Analyze TE distribution across genomic features
* __Inputs:__
    * `genomic.gff`, `TE.txt`, `genome.fa.fai`
    * `promoter.bed`: generated automatically in Step 0
* __Outputs:__
    * `OUTPUT_1_TE_distribution_enrichment.png`, `OUTPUT_1_TE_distribution_enrichment.png`
* __Usage:__
  ```bash
  ./petem --1 \
    -g <MODULE0_ANNOTATION_DIR> \
    -t <TE_BED> \
    -f <GENOME_FAI> \
    -o <OUTPUT_DIR>
  ```
* __Parameters:__

  | Parameter | Description |
  |---|---|
  | `--1` | Run module 1. Equivalent to `petem 0_preprocessing`. |
  | `-g` | Path to the annotation directory generated by module 0 (`module_0_annotation`). |
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-f` | Path to the genome FASTA index file (`.fai`). |
  | `-o` | Output directory for module 1 results. |
  
### Module 2. Enriched promoter-embedded TE Families
> Identify enriched TE families overlapping with promoters
* __Inputs:__
    * `TE.txt`, `promoter.bed` (from Step 0)
* __Outputs:__
    * `OUTPUT_2_Promoter_embedded_TE_family_enrichment.png`, `OUTPUT_2_Promoter_embedded_TE_family.txt`
* __Usage:__
  ```bash
  ./petem --2 \
    -t <TE_BED> \
    -T <TE_FAMILY_TABLE> \
    -i <TE_OVERLAP_PROMOTER_BED> \
    -o <OUTPUT_DIR>
  ```
* __Parameters:__

  | Parameter | Description |
  |---|---|
  | `--2` | Run module 2. |
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-T` | Path to the transposable element family annotation table. |
  | `-i` | Path to the `TE_overlap_promoter.bed` file generated by module 0. |
  | `-o` | Output directory for module 2 results. |

### Module 3. TE methylation near gene
> Visualize distance impact of TE methylation on gene expression
* __Inputs:__
    * `gene.bed`(from Step 0), `TE.txt`, `gene_expression.txt` and `TE_expression.txt`
* __Outputs:__
    * `OUTPUT_3_gene_TE_number.txt`, `OUTPUT_3_gene_proximal_TE_*.png`
    
* __Usage:__
  ```bash
  ./petem --3 \
    -g <MODULE0_ANNOTATION_DIR> \
    -t <TE_BED> \
    -eg <GENE_EXPRESSION> \
    -et <TE_EXPRESSION> \
    -m <CGMAP_FILES> \
    -o <OUTPUT_DIR> \
    -d <PROMOTER_DISTANCE> \
    -p <MIN_COVERAGE> \
    -w <WINDOW_SIZE> \
    -c <CONTROL_MODE> \
    -l <CONTROL_TE_CLASS>
  ```

* __Parameters:__

  | Parameter | Description |
  |---|---|
  | `--3` | Run module 3. |
  | `-g` | Path to the annotation directory generated by module 0 (`module_0_annotation`). |
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-eg` | Path to the gene expression table. |
  | `-et` | Path to the transposable element expression table. |
  | `-m` | One or more CGmap files, separated by spaces. Compressed files such as `.CGmap.gz` are supported. |
  | `-o` | Output directory for module 3 results. |
  | `-d` | The total up- and downstream distance (in bp) to consider gene-proximal TEs. (Default `4000` means ±4 kb around genes will be analyzed) |
  | `-p` | The top/bottom X% of genes will be considered as the highly/lowly expressed genes. (Default: `10`, means top/bottom 10% of genes will be the highly/lowly expressed genes) |
  | `-w` | For choosing `average`. Sliding window size (bp) used to smooth the TE methylation level curve. (Default: `100`) |
  | `-c` | Enable (`y`) or disable (`n`) the control analysis. |
  | `-l` | There are several options to show the pattern of gene-proximal TE methylation, including (1) average TE methylation within each window, (2) linear regression line, (3) second-degree polynomial regression line, and (4) local regression line (`average`, `linear`, `poly2`, and `poly`) |

### Module 4. Correlation between methylation and expression
> Correlate gene expression with TE/promoter methylation and TE expression with TE methylation
* __Inputs:__
    * `gene_expression.txt` and `TE_expression.txt`
* __Outputs:__
    * `OUTPUT_4_geneexp/TEexp_*.png`, `OUTPUT_4_*_correlation_*.png`
* __Usage:__
  ```bash
  ./petem --4 \
    -eg <GENE_EXPRESSION> \
    -et <TE_EXPRESSION> \
    --module0-dir <MODULE0_OUTPUT_DIR> \
    -o <OUTPUT_DIR>
  ```

* __Parameters:__

  | Parameter | Description |
  |---|---|
  | `--4` | Run module 4. |
  | `-eg` | Path to the gene expression table. |
  | `-et` | Path to the transposable element expression table. |
  | `--module0-dir` | Path to the module 0 output directory. Required for loading methylation and annotation results generated in module 0. |
  | `-o` | Output directory for module 4 results. |
  | `--window-number` | Number of sliding windows used to smooth the correlation curves. Default: `100`. |
  | `--ylim_CG` | Maximum y-axis value for gene expression vs. TE/promoter CG methylation correlation plots. Default: `80`. |
  | `--ylim_CHG` | Maximum y-axis value for gene expression vs. TE/promoter CHG methylation correlation plots. Default: `40`. |
  | `--ylim_CHH` | Maximum y-axis value for gene expression vs. TE/promoter CHH methylation correlation plots. Default: `40`. |
  | `--ylim_TEexpTEmC_CH` | Maximum y-axis value for TE expression vs. TE CHG/CHH methylation correlation plots. Default: `40`. |
  | `--ylim_TEexpTEmC_CG` | Maximum y-axis value for TE expression vs. TE CG methylation correlation plots. Default: `80`. |


### Module 5. Associated TE and gene pairs
> Examine the correlations between changes in TE methylation, TE expression, and gene expression across different conditions.
* __Inputs:__
    * `gene_expression.txt`, `TE_expression.txt`
* __Outputs:__
    * `OUTPUT_5_*_scatter.png & .txt`, `OUTPUT_5_Q2/4_boxplot_*.png`, `OUTPUT_5_Q2/4_*.txt`
* __Usage:__
  ```bash
  ./petem --5 \
    --DEG <DEG_TABLE> \
    --DETE <DETE_TABLE> \
    --module0-dir <MODULE0_OUTPUT_DIR> \
    -o <OUTPUT_DIR>
  ```

* __Parameters:__

  | Parameter | Description |
  |---|---|
  | `--5` | Run module 5. |
  | `--DEG` | Path to the differentially expressed gene (DEG) table. |
  | `--DETE` | Path to the differentially expressed transposable element (DETE) table. |
  | `--module0-dir` | Path to the module 0 output directory. Required for loading methylation and annotation results generated in module 0. |
  | `-o` | Output directory for module 5 results. |
