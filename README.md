# Promoter-embedded TE Methylation (PeTEM) analyzer
PeTEM is designed to analyze the association between **promoter-embedded TE methylation** and **neighboring gene expression**. It integrates genome annotation, methylome, and transcriptome data to evaluate genome-wide correlations between TE methylation and gene expression, and identify TE–gene pairs showing coordinated methylation and expression changes across conditions.

<img width="2759" height="1432" alt="image" src="https://github.com/user-attachments/assets/30d94baf-d631-4789-9318-ae78b7744540" />


## Tutorial
See the [tutorial](https://github.com/PaoyangLab/PeTEM/blob/main/Tutorial.md) for an example workflow.


## Table of Contents

- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Input Files](#input-files)
- [Pipeline Modules](#pipeline-modules)
  - [Module 0. Preprocessing](#module-0-preprocessing)
  - [Module 1. Profile TE genomic distribution](#module-1-profile-te-genomic-distribution)
  - [Module 2. Identify enriched promoter-embedded TE families](#module-2-identify-enriched-promoter-embedded-te-families)
  - [Module 3. Visualize TE methylation near gene](#module-3-visualize-te-methylation-near-gene)
  - [Module 4. Calculate correlation coefficients](#module-4-calculate-correlation-coefficients)
  - [Module 5. Identify associated TE and gene pairs](#module-5-identify-associated-te-and-gene-pairs)


## System Requirements

  - **Runtime dependencies** 
    - Bash
    - Perl (v5.22.0+)
    - gzip / gunzip
    - awk / sort / uniq

  - **Environment setup**
  
    > The following sections show the required environments and packages. <br>
    > In [installation](#installation) section, we provide three alternative methods for setting up the environment.
  
    <details>
    <summary> 👉 <b>R version ≥ 4.2 (tested on 4.3.2)</b></summary>
    <br> 
      
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
    <summary> 👉 <b>Python version ≥ 3.8 (tested on 3.8.10)</b></summary>
    <br> 
      
    * pandas (≥ 1.2.4)
  
    </details>
  
    <details>
    <summary> 👉 <b>Bioinformatics tools</b></summary>
    <br> 
      
    *	samtools (tested on 1.10)
    *	bedtools (tested on v2.27.1)
    *	wigToBigWig (bundled in this repository)
    *	bigWigAverageOverBed (bundled in this repository)
  
    </details>


## Installation

  - **Clone repository**
     
    ```bash
    git clone https://github.com/PaoyangLab/PeTEM.git
    cd PeTEM
    ```
  
  - **Download example data**
     
    ```bash
    wget https://github.com/PaoyangLab/PeTEM/releases/download/v0.0.3/PeTEM_data.tar.gz
    tar -xzvf PeTEM_data.tar.gz
    ```

  - **Set up environment**
    > Choose **one** of the following installation methods to set up the PeTEM environment.
  
    - **Conda setup (👍recommended)**
      > Use the provided environment definition:
      ```bash
      conda env create -f environment.yml
      conda activate petem
      bash env_check.sh ##optional
      ```

    - **Docker image**
      ```bash
      docker build -t petem:local .
      docker run --rm petem:local --help
      ```
  
    - **Local setup**
      > The script installs dependencies using `apt-get`, `pip3 --user`, and `Rscript`, then runs `bash env_check.sh`. <br>
      > If `apt-get` is not available the script prints the package list to install manually.
      ```bash
      bash setup.sh
      bash env_check.sh ##optional
      ```

## Input Files
PeTEM integrates **genome annotations**, **DNA methylation data**, and **expression data** as input files. 
In PeTEM, modules 1 and 2 rely solely on annotation data, while running the remaining modules additionally require methylation and expression data.

- **Genome Annotation**

  - **General features annotations file (GFF format)** 
    > The genome annotation file contains genomic feature coordinates and hierarchical annotations.
  
    Format: (`genomic.gff`)
    
    | Column | Description |
    |----------|----------|
    | seqid | Sequence ID, chromosome, or scaffold name |
    | source | Annotation source |
    | type | Genomic feature type |
    | start | Start coordinate (1-based) |
    | end | End coordinate |
    | score | Annotation score |
    | strand | Strand information (`+`, `-`, or `.`) |
    | phase | CDS phase information |
    | attributes | Feature attributes and hierarchical annotation information |
  
    Example:
  
    ```text
    Chr1  Araport11  gene             3631  5899  .  +  .  ID=AT1G01010;Name=AT1G01010;full_name=NAC domain containing protein 1
    Chr1  Araport11  mRNA             3631  5899  .  +  .  ID=AT1G01010.1;Name=AT1G01010.1;Parent=AT1G01010
    Chr1  Araport11  CDS              3760  3913  .  +  0  ID=AT1G01010:CDS:1;Parent=AT1G01010.1
    Chr1  Araport11  exon             3631  3913  .  +  .  ID=AT1G01010:exon:1;Parent=AT1G01010.1
    Chr1  Araport11  five_prime_UTR   3631  3759  .  +  .  ID=AT1G01010:five_prime_UTR:1
    Chr1  Araport11  three_prime_UTR  5631  5899  .  +  .  ID=AT1G01010:three_prime_UTR:1
    ```

    <details>
    <summary> ℹ️ Click for file sources </summary>
      
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
    * Fungi: [10 species](https://urgi.versailles.inra.fr/download/fungi/TEs/)
    
        | File ID | Species |
        |--------------|---------|
        | BCINB05_10_18Chr / Bcin0510 / BcinT4 | *Botrytis cinerea* |
        | Bgra | *Blumeria graminis* |
        | Colletotrichum_higginsianum | *Colletotrichum higginsianum* |
        | Lmac | *Leptosphaeria maculans* |
        | Mlar | *Melampsora larici-populina* |
        | Mory | *Magnaporthe oryzae* |
        | Mvio | *Microbotryum violaceum* |
        | Pgra | *Puccinia graminis sp. Triticeae* |
        | Sscl | *Sclerotinia sclerotiorum* |
        | Tmel | *Tuber melanosporum* |
    
    </details>

  - **Genome index file (FASTA index)** 
  
    > The genome FASTA index file provides the names and lengths of each chromosome.
  
    Format: (`genome.fa.fai`)
  
    | Column | Description |
    |----------|----------|
    | name | Chromosome or sequence name |
    | length | Sequence length (bp) |
    | offset | Byte offset of the sequence in the FASTA file |
    | linebases | Number of bases per sequence line |
    | linewidth | Number of bytes per sequence line, including newline characters |
  
    Example:
  
    ```text
    Chr1  30427671  74         79  80
    Chr2  19698289  30812981   79  80
    Chr3  23459830  50760691   79  80
    Chr4  18585056  74517556   79  80
    Chr5  26975502  93337941   79  80
    ChrC  154478    120654981  79  80
    ChrM  367808    120811562  70  71
    ```

    
    <details>
    <summary> ℹ️ Click for file sources </summary>
    
    > FASTA index file is generated from genome fasta file using samtools:
    ```
    samtools faidx genome.fa
    ```
    
    </details>
    
  
  - **Transposable element coordinates**  
  
    > The transposable element annotation file contains genomic coordinates and classification information for annotated transposable elements (TEs), including TE family and strand orientation.
    
    Format: (`TE.txt`)
  
    | Column | Description |
    |----------|----------|
    | TE name | Transposable element identifier |
    | chromosome | Chromosome or scaffold name |
    | start | Start coordinate (0-based) |
    | end | End coordinate |
    | score | Annotation score |
    | strand | Strand information (`+`, `-`, or `.`) |
    | TE family | Transposable element family classification |
  
    Example:
  
    ```text
    AT1TE00010  Chr1  11897  11976  0  +  LTR/Copia
    AT1TE00020  Chr1  16883  17009  0  -  RC/Helitron
    AT1TE00025  Chr1  17024  18924  0  +  RC/Helitron
    AT1TE00030  Chr1  18331  18642  0  -  DNA/HAT
    ```

    <details>
    <summary> ℹ️ Click for file sources </summary>
    
    > Processed TE annotation files from 18 species are available [here](https://github.com/PaoyangLab/PeTEM/releases/tag/TE_annotation)
 
    **Download files**
    ```bash
    wget https://github.com/PaoyangLab/PeTEM/releases/download/TE_annotation/TE_files.tar.gz
    tar -xzvf TE_files.tar.gz
    ```
    
    **Included species**
    | Kingdom | Species | File |
    |----------|----------|----------|
    | Animal | *Homo sapiens* | `human_TE.txt` |
    | Animal | *Mus musculus* | `mouse_TE.txt` |
    | Animal | *Danio rerio* | `zebrafish_TE.txt` |
    | Animal | *Drosophila melanogaster* | `fruit_fly_TE.txt` |
    | Plant | *Arabidopsis thaliana* | `Arabidopsis_TE.txt` |
    | Plant | *Oryza sativa* | `rice_TE.txt` |
    | Plant | *Zea mays* | `maize_TE.txt` |
    | Plant | *Glycine max* | `soybean_TE.txt` |
    | Fungi | *Botrytis cinerea* | `Botrytis_cinerea_TE.txt` |
    | Fungi | *Blumeria graminis* | `Blumeria_graminis_TE.txt` |
    | Fungi | *Colletotrichum higginsianum* | `Colletotrichum_higginsianum_TE.txt` |
    | Fungi | *Leptosphaeria maculans* | `Leptosphaeria_maculans_TE.txt` |
    | Fungi | *Melampsora larici-populina* | `Melampsora_larici-populina_TE.txt` |
    | Fungi | *Magnaporthe oryzae* | `Magnaporthe_oryzae_TE.txt` |
    | Fungi | *Microbotryum violaceum* | `Microbotryum_violaceum_TE.txt` |
    | Fungi | *Puccinia graminis* f. sp. *tritici* | `Puccinia_graminis_TE.txt` |
    | Fungi | *Sclerotinia sclerotiorum* | `Sclerotinia_sclerotiorum_TE.txt` |
    | Fungi | *Tuber melanosporum* | `Tuber_melanosporum_TE.txt` |

    </details>

- **Methylation Data**  
  > CGmap is a file format for storing single-base resolution DNA methylation data.


  Format: (`*.CGmap.gz files`)

  | Column | Description |
  |----------|----------|
  | chromosome | Chromosome or scaffold name |
  | nucleotide | Cytosine (`C`) or guanine (`G`) on the forward or reverse strand |
  | position | Genomic position of the methylation site (1-based) |
  | context | Methylation context (`CG`, `CHG`, or `CHH`) |
  | dinucleotide | Dinucleotide sequence surrounding the methylation site |
  | methylation level | Methylation level ranging from 0 to 1 |
  | #C site | Number of reads supporting methylated cytosine |
  | #C+T site | Total sequencing depth at the site |

  Example:

  ```text
  Chr3  C  556  CG   CG  0.877551  43  49
  Chr3  G  557  CG   CG  0.787879  26  33
  Chr3  G  558  CHG  CC  0.405405  15  37
  Chr3  G  560  CHH  CA  0.102564   4  39
  ```

  <details>
    
    <summary> ℹ️ Click for file sources </summary>
    
    > CGmap files can be generated from WGBS data using:
    * [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)
    * [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) + [CGmapTools](https://cgmaptools.github.io/)

  </details>

- **Expression Data**

  > Expression data consists of gene expression (`gene_expression.txt`) and transposable element expression (`TE_expression.txt`) tables.

  > Each row represents a gene or TE. The first columns contain the average expression level (RPKM) for each condition, followed by differential expression statistics, including log2 fold change (logFC), p-value, and false discovery rate (FDR) for pairwise condition comparisons.

  Column naming convention

  | Column pattern | Description | Example |
  |----------|----------|----------|
  | condition_name | Average expression level (RPKM) of a condition | `root`, `leaf` |
  | logFC_condition1_condition2 | Log2 fold change between two conditions | `logFC_root_leaf` |
  | PValue_condition1_condition2 | Statistical significance of differential expression | `PValue_root_leaf` |
  | FDR_condition1_condition2 | Multiple-testing adjusted p-value | `FDR_root_leaf` |

  Format: 

  | Column | Description |
  |----------|----------|
  | Row name | Gene/Transposon identifier |
  | condition_name | Average expression level (RPKM) of each condition |
  | logFC_condition1_condition2 | Log2 fold change between two conditions |
  | PValue_condition1_condition2 | Differential expression p-value |
  | FDR_condition1_condition2 | Adjusted p-value (FDR) |

  Example: (`gene_expression.txt`)

  ```text
  AT1G01010  11.64  10.93   0.09   0.60  1
  AT1G01020   6.39   5.95   0.10   0.62  1
  AT1G01030   0.65   1.01  -0.63   0.19  0.84
  ```

  Example: (`TE_expression.txt`)

  ```text
  AT1TE00010  386.06  240.27  0.63  0.73  1
  AT1TE00020    0.00    0.00  0.00  1.00  1
  AT1TE00025    0.00    0.00  0.00  1.00  1
  ```

  <details>
    <summary> ℹ️ Click for file sources </summary>
    
    > Expression files can be generated from expression/count matrices using `edgeR_for_expression.R` bundled in this repository
    ```
    Rscript differential_expression.R \
      --eg gene_expression.txt \
      --et TE_expression.txt \
      -o results
    ```

  </details>

## Pipeline Modules
Users must run module 0 at the first time to preprocess the input files before running module 1 to 5. Module 1 to 5 are independent and not sequential. 

<img width="1197" height="837" alt="image" src="https://github.com/user-attachments/assets/ab60c2e9-87fe-44c6-a000-1b3cbe0a410d" />


### Module 0. Preprocessing
> To prepare the inputs required for all downstream modules, module 0 generates promoter regions (`promoter.bed`) and integrate methylation and expression data.

* __Inputs:__
    * `genomic.gff`, `TE.txt`, `genome.fa.fai`, `gene_expression.txt`, `TE_expression.txt`, `*.CGmap.gz`
* __Outputs:__
    * `OUTPUT_0_embedded_TE_gene_number.txt`, `promoter.bed`, `TE_overlap_promoter.bed`, `Tab_*.txt`
    * `PETEM_MODULE0_MANIFEST.json`: Module 0 generates a json file for following module to simplify the command line of input path. 

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
  
  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Required&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `--0` | Run module 0. |
  | `-g` | Path to the genome annotation file in GFF3 format. |
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-eg` | Path to the gene expression table. |
  | `-et` | Path to the transposable element expression table. |
  | `-f` | Path to the genome FASTA index file (`.fai`). This is used to automatically generate gene, promoter, intergenic region, and other BED files. |
  | `-m` | One or more CGmap files, separated by spaces. Compressed files such as `.CGmap.gz` are supported. |

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Optional&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `-up` | Promoter region at the upstream of TSS. (Default: `1500`) |
  | `-dn` | Promoter region at the downstream of TSS. (Default: `500`) |
  | `-o` | Output directory for module 0 results. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|


### Module 1. Profile TE genomic distribution
> Analyze TE distribution across genomic features.

> This reveals whether TEs preferentially accumulate in specific genomic regions.

* __Inputs:__
    * `genomic.gff`, `TE.txt`, `genome.fa.fai`
    * `promoter.bed` (from Module 0)
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

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Required&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `--1` | Run module 1. Equivalent to `petem 0_preprocessing`. |
  | `--manifest`  | Use manifest from module 0 to simplify input path. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Optional&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `-g` | Path to the annotation directory generated by module 0 (`module_0_annotation`). |
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-f` | Path to the genome FASTA index file (`.fai`). |
  | `-o` | Output directory for module 1 results. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|


  
### Module 2. Identify enriched promoter-embedded TE families
> Identify enriched TE families overlapping with promoters.

> This identifies TE families enriched in promoters and highlights candidates that may influence host gene expression.

* __Inputs:__
    * `TE.txt`
    * `promoter.bed` (from Module 0)
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

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Required&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `--2` | Run module 2. |
  | `--manifest`  | Use manifest from module 0 to simplify input path. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Optional&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-T` | Path to the transposable element family annotation table. |
  | `-i` | Path to the `TE_overlap_promoter.bed` file generated by module 0. |
  | `-o` | Output directory for module 2 results. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|


### Module 3. Visualize TE methylation near gene
> Visualize distance impact of TE methylation on gene expression.
 
> This provides a genome-wide view of how nearby methylated TEs influence neighboring gene expression.

* __Inputs:__
    * `TE.txt`, `gene_expression.txt`, `TE_expression.txt`
    * `gene.bed`(from Module 0)
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

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Required&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `--3` | Run module 3. |
  | `--manifest`  | Use manifest from module 0 to simplify input path |
  | `-m` | One or more CGmap files, separated by spaces. Compressed files such as `.CGmap.gz` are supported. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Optional&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `-g` | Path to the annotation directory generated by module 0 (`module_0_annotation`). |
  | `-t` | Path to the transposable element annotation file in BED format. |
  | `-eg` | Path to the gene expression table. |
  | `-et` | Path to the transposable element expression table. |
  | `-o` | Output directory for module 3 results. |
  | `-d` | The total up- and downstream distance (in bp) to consider gene-proximal TEs. (Default `4000` means ±4 kb around genes will be analyzed) |
  | `-p` | The top/bottom X% of genes will be considered as the highly/lowly expressed genes. (Default: `10`, means top/bottom 10% of genes will be the highly/lowly expressed genes) |
  | `-w` | For choosing `average`. Sliding window size (bp) used to smooth the TE methylation level curve. (Default: `100`) |
  | `-l` | There are several options to show the pattern of gene-proximal TE methylation, including (1) average TE methylation within each window, (2) linear regression line, (3) second-degree polynomial regression line, and (4) local regression line (`average`, `linear`, `poly2`, and `poly`). Default: `average` |
  | `-CI` | To show the 95% CI on the plot or not. Default: `no` |
  | `-border` | To show the border that TE methylation near lowly expressed genes is significantly higher than that near highly expressed genes on the plot or not. Default: `no` |
  | `-unexp` | To include TEs with zero expression across all samples in the analysis. Default: `no` |
  | `-nTE` | To show the methylation level of non transposon sites around genes. Default: `no` |

### Module 4. Calculate correlation coefficients
> Correlate gene expression with TE/promoter methylation, and TE expression with TE methylation.

> This module quantifies the genome-wide regulatory effects of TE methylation.

* __Inputs:__
    * `gene_expression.txt`, `TE_expression.txt`
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

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Required&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `--4` | Run module 4. |
  | `--manifest`  | Use manifest from module 0 to simplify input path |


  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Optional&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `-eg` | Path to the gene expression table. |
  | `-et` | Path to the transposable element expression table. |
  | `--module0-dir` | Path to the module 0 output directory. Required for loading methylation and annotation results generated in module 0. |
  | `-o` | Output directory for module 4 results. |
  | `--smooth` | Adjust the degree of curve smoothing. Values range from `1` to `5`; `1` preserves the original curve shape (no smoothing), while `5` produces the smoothest curve. Default: `3`.|

### Module 5. Identify associated TE and gene pairs
> Examine the correlations between changes in TE methylation, TE expression, and gene expression across different conditions.

> This identifies TE–gene pairs whose condition-specific expression changes are associated with TE methylation dynamics.

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

  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Required&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `--5` | Run module 5. |
  | `--manifest`  | Use manifest from module 0 to simplify input path &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|


  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Optional&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
  |---|---|
  | `--DEG` | Path to the differentially expressed gene (DEG) table. |
  | `--DETE` | Path to the differentially expressed transposable element (DETE) table. |
  | `--module0-dir` | Path to the module 0 output directory. Required for loading methylation and annotation results generated in module 0. |
  | `-o` | Output directory for module 5 results. |
  | `--positive` | Specify the correlation direction to analyze. Set to `yes` for positive correlation analysis or `no` for negative correlation analysis. Default: `yes`.|

