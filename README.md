# Promoter-embedded TE Methylation (PeTEM) analyzer
PeTEM is designed to analyse the association between promoter-embedded TE methylation and neighbouring gene expression. It integrates genome annotations with methylome and transcriptome data, enabling users to evaluate the genome-wide correlations between TE methylation and gene expression, and identify specific TE and gene pairs with associated methylation and expression changes across different conditions (such as stages, treatments, or phenotypes).

<img width="865" height="439" alt="image" src="https://github.com/user-attachments/assets/56516953-0937-4c9c-951f-f474a6b0ee67" />


## Tutorial
Please follow the [tutorial](https://github.com/PaoyangLab/PeTEM/blob/main/Tutorial.md) of example use case.

## Installation

#### Clone the repository from Github:
```bash
git clone https://github.com/yc811/PeTEM.git
cd PeTEM
```

## System requirements

#### Set up environment through one of following methods 
<details>
<summary>Docker</summary>

```bash
docker build -t petem:local .
docker run --rm petem:local --help
```

</details>

<details>
<summary>Anaconda</summary>

```bash
conda env create -f environment.yml
conda activate petem
bash env_check.sh ##optional
```

</details>

<details>
<summary>Local</summary>

```bash
bash setup.sh
bash env_check.sh ##optional
```

</details>

#### Details of required environment and packages

<details>
<summary>R version ≥ 4.2 </summary>
  
(packages tested on 4.3.2)
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

</details>

<details>
<summary>Python version ≥ 3.8 </summary>
  
(packages tested on 3.8.10)
* pandas (≥ 1.2.4)
* built-in packages:
  * glob
  * os
  * time

</details>

<details>
<summary>Bioinformatics tools</summary>
   
*	samtools (tested on 1.10)
*	bedtools (tested on v2.27.1)
*	wigToBigWig
*	bigWigAverageOverBed

</details>

## Input Files
PeTEM integrates inputs data including [genome annotations](#genome-annotation), [genome-wide DNA methylation](#methylation-data), and [expression data](#expression-data). 
In PeTEM, running the first two modules rely solely on annotation data, while running the remaining modules additionally require methylation and expression data.

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

#### Gene annotation file (FASTA index) 
> The genome FASTA index file is generated from genome fasta file (Usage: `samtools faidx genome.fa`), providing the names and lengths of each chromosome.

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

#### Transposable element coordinates (BED format) 
> The transposable element annotation file contains genomic coordinates and classification information for annotated transposable elements (TEs), including TE family and strand orientation.

`TE.txt` 
| TE name | chromosome | start | end | score | strand | TE family |
|---|---|---|---|---|---|---|
| AT1TE00010 | Chr1 | 11897 | 11976 | 0 | + | LTR/Copia |
| AT1TE00020 | Chr1 | 16883 | 17009 | 0 | - | RC/Helitron |
| AT1TE00025 | Chr1 | 17024 | 18924 | 0 | + | RC/Helitron |
| AT1TE00030 | Chr1 | 18331 | 18642 | 0 | - | DNA/HAT |

<details>
<summary>Sources of commonly used genome annotations</summary>
   
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

   

### Methylation Data
* *.CGmap.gz files – Per-sample methylation CGmap files
   > CGmap files includes 8 columns: chromosome, C or G (forward or reverse strand), position, context (CG/CHG/CHH), dinucleotide, methylation level (0-1), # of reads supporting methylation, depth


| chromosome | nucleotide | site | context | dinucleotide | methylation level | C site | C+T site |
|---|---|---|---|---|---|---|---|
| Chr3 | C | 556 | CG | CG | 0.877551 | 43 | 49 |
| Chr3 | G | 557 | CG | CG | 0.787879 | 26 | 33 |
| Chr3 | G | 558 | CHG | CC | 0.405405 | 15 | 37 |
| Chr3 | G | 560 | CHH | CA | 0.102564 | 4 | 39 |


### Expression Data
* gene_expression.txt – Differentially expressed genes (Step 0, 3-2, 4, 5)
* TE_expression.txt – Differentially expressed TEs (Step 0, 3-2, 4, 5)
   > In the expression files, the row names are the gene and TE names, followed by columns showing average expression level (RPKM) of each conditions. The rest of columns shows log2 fold change, p value, and FDR comparing each two conditions.

   > The column names should be: conditions names, "logFC_condition1_condition2", "PValue_condition1_condition2", "FDR_condition1_condition2", "logFC_condition2_condition3", ... etc.

|  | WT | drdd | logFC_drdd_WT | PValue_drdd_WT | FDR_drdd_WT |
|---|---|---|---|---|---|
| AT1G01010 | 1.58 | 2.39 | 0.61 | 0.21 | 1 |
| AT1G01020 | 5.60 | 5.43 | -0.04 | 0.87 | 1 |
| AT1G01030 | 1.49 | 3.39 | 1.17 | 0.01 | 0.12 |

   These files will automatically be converted into expression matrices (gene_expression.txt, TE_expression.txt) for steps 0, 3-2, and 4.



## Pipeline Modules
Upon running run_pipeline.sh, users will be asked which modules to execute (y/n), and need to provide all required files and parameters upfront. Users must run module 0 at the first time to preprocess the input files, before running module 1 to 5. Module 1 to 5 are dependent and not sequential. 

### Module 0. Preprocessing
* Generate promoter regions (promoter.bed)
* Integrate methylation and expression data
* __Inputs:__
    * `genomic.gff`, `TE.txt`, `genome.fa.fai`, `gene_expression.txt`, `TE_expression.txt`, `*.CGmap.gz`
* __Outputs:__
    * `OUTPUT_0_embedded_TE_gene_number.txt`, `promoter.bed`, `TE_overlap_promoter.bed`, `Tab_*.txt`
* __Parameters:__ 
    * __Promoter region:__ The default promoter is defined as `1500` bp upstream to `500` bp downstream from the transcription start site (TSS). Users can customize this range by entering other upstream/downstream length from TSS. 

### Module 1. TE Distribution
* Analyze TE distribution across genomic features
* __Inputs:__
    * `genomic.gff`, `TE.txt`, `genome.fa.fai`
    * `promoter.bed`: generated automatically in Step 0
* __Outputs:__
    * `OUTPUT_1_TE_distribution_enrichment.png`, `OUTPUT_1_TE_distribution_enrichment.png`
  
### Module 2. Promoter-embedded TE Families
* Identify enriched TE families overlapping with promoters
* __Inputs:__
    * `TE.txt`, `promoter.bed` (from Step 0)
* __Outputs:__
    * `OUTPUT_2_Promoter_embedded_TE_family_enrichment.png`, `OUTPUT_2_Promoter_embedded_TE_family.txt`
  
### Module 3. Gene-proximal TE methylation
* 3-1 Preprocessing: Prepare methylation files required in Step 3-2
* 3-2 Plotting: Visualize distance impact of TE methylation on gene expression[
* __Inputs:__
    * `gene.bed`(from Step 0), `TE.txt`, `gene_expression.txt` and `TE_expression.txt`
* __Outputs:__
    * `OUTPUT_3_gene_TE_number.txt`, `OUTPUT_3_gene_proximal_TE_*.png`
* __Parameters:__ 
    * __Limit range:__ The total up- and downstream distance (in bp) to consider gene-proximal TEs. (Default `4000` means ±4 kb around genes will be analyzed)
    * __Tick size:__ The spacing (bp) between x-axis ticks in the resulting plot. Recommended: approximately 1/3 to 1/4 of the limit range. (Default: `1000` )
    *	__Percent of genes__: The top/bottom X %of genes will be considered as the highly/lowly expressed genes. (Default: `10`, means top/bottom 10% of genes will be the highly/lowly expressed genes)
    *	__Point__: To show the points representing the methylation level of each TE site on the plot or not. (Default: `yes`)
    *	__95% CI__: To show the 95% CI on the plot or not. (Default: `no`)
    *	__Pattern of TE methylation__: There are several options to show the pattern of gene-proximal TE methylation, including (1) average TE methylation within each window, (2) linear regression line, (3) second-degree polynomial regression line, and (4) local regression line (`average`, `linear`, `poly2`, and `loess`) For choosing average TE methylation, the users need to setup the window, while for choosing local regression, the users need to setup the span. (Default: `average`)
    *	__Window size__: For choosing `average`. Sliding window size (bp) used to smooth the TE methylation level curve. (Default: `100`)
    *	__Window step__: For choosing `average`. The length of each step of window (bp). If the step = window size means that it just calculates the average TE methylation of each non–overlapped bin. If the step < window size means that it is a sliding window. (Default: `100`)
    *	__Span__: For choosing `loess`. This parameter which represents the proportion of total data points used to compute each smoothed value, controls the degree of smoothing. (Default: `0.02`)


### Module 4. Correlation (Single Condition)
* Correlate gene expression with TE/promoter methylation and TE expression with TE methylation
* __Inputs:__
    * `gene_expression.txt` and `TE_expression.txt`
* __Outputs:__
    * `OUTPUT_4_geneexp/TEexp_*.png`, `OUTPUT_4_*_correlation_*.png`
* __Parameters:__ 
    * __Window number:__ Number of sliding windows used to smooth the correlation curves. (Default: `100`).
    * __Y-axis limits:__ Controls the maximum value shown in the y-axis of each correlation plot:
        * ylim_CG: gene expression vs TE/promoter CG methylation (Default: `80`)
        * ylim_CHG: gene expression vs TE/promoter CHG methylation (Default: `40`)
        * ylim_CHH: gene expression vs TE/promoter CHH methylation (Default: `40`)
        * ylim_TEexpTEmC_CH: TE expression vs TE CHG/CHH methylation (Default: `40`)
        * ylim_TEexpTEmC_CG: TE expression vs TE CG methylation (Default: `80`)

### Module 5. Correlation (Across Conditions)
* Examine the correlations between changes in TE methylation, TE expression, and gene expression across different conditions.
* __Inputs:__
    * `gene_expression.txt`, `TE_expression.txt`
* __Outputs:__
    * `OUTPUT_5_*_scatter.png & .txt`, `OUTPUT_5_Q2/4_boxplot_*.png`, `OUTPUT_5_Q2/4_*.txt`
### Additional parameter for optionally include unexpressed TEs
* This parameter determines whether to include TEs with zero expression across all samples in the analysis (all plots of Step 3; correlation between gene expression and TE methylation in Step 4 and Step 5). The default is “no (n)”. Setting this to “yes (y)” includes all TEs regardless of expression level, which may increase the number of analyzed TEs but also add background noise. Recommended to keep “no (n)” unless specifically investigating silent or lowly expressed TEs.

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
Genomic features annotaion file: genomic.gff
TE BED file: TE.bed
Genome fasta index: genome.fa.fai
Gene expression file: gene_expression.txt
TE expression file: TE_expression.txt
Methylation CGmap.gz files (space separated): WT_01.CGmap.gz WT_02.CGmap.gz...

Include unexpressed TEs? (y/n, default n): n
Promoter upstream length from TSS (default 1500): 1500
Promoter downstream length from TSS (default 500): 500
Limit up-/down-stream range (bp)(default 4000): 4000
Tick size (bp)(default 1000): 1000
Window size (bp)(default 200): 200
Window number (default 100): 100
y-axis limit for gene expression vs TE/promoter mC plot (CG, default 80): 80
y-axis limit for gene expression vs TE/promoter mC plot (CHG, default 40): 40
y-axis limit for gene expression vs TE/promoter mC plot (CHH, default 40): 40
y-axis limit for TE expression vs TE mC plot (CH, default 40): 40
y-axis limit for TE expression vs TE mC plot (CG, default 80): 80
```

