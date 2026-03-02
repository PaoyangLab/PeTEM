# Tutorial
This tutorial explains how to run the PeTEM pipeline using the example files provided in the PeTEM_data folder.

## Installation
```
# Commands to be updated. Add commands for installation and environment check.
```
#### Download from github
```bash
git clone https://github.com/yc811/PeTEM.git
cd PeTEM
```
This will download all scripts and the `run_PeTEM.sh` interactive pipeline.

#### Download the required tools: wigToBigWig, bigWigAverageOverBed
```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig
chmod +x wigToBigWig

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigAverageOverBed
chmod +x bigWigAverageOverBed
```

##  Example Data
We provide example input files under [`PeTEM_data/`](PeTEM_data) to test the pipeline.
| File                                                                           | Description                     |
| ------------------------------------------------------------------------------ | ------------------------------- |
| `genome.fa.fai`                                                                | Genome index file               |
| `TE.txt`                                                                       | TE annotation file                   |
| `gene_expression.txt`                                                          | Gene expression file            |
| `TE_expression.txt`                                                            | TE expression file              |
| `leaf_01.CGmap.gz`, `leaf_02.CGmap.gz`, `root_01.CGmap.gz`, `root_02.CGmap.gz` | Methylation CGmap files         |


## Running the pipeline
PeTEM provides an interactive bash script `run_PeTEM.sh` to select steps and input files.

### Step 0: Link the example data to the current path
```
ln -sf PeTEM_data/* .
```

### Step 1: Launch the interactive pipeline
```
bash run_PeTEM.sh
```
You will be prompted to select which modules to run (enter `y` or `n`):
```
0. Preprocessing? (y/n):
1. TE distribution? (y/n):
2. Promoter-embedded TE families? (y/n):
3-1. Gene-proximal TE methylation: preprocessing? (y/n):
3-2. Gene-proximal TE methylation: plot? (y/n):
4. Correlation single condition? (y/n):
5. Correlation across conditions? (y/n):
```

## Step 2: Input files and parameters
After selecting modules, you will be prompted to provide paths to necessary files.
Use the example files in PeTEM_data:
```
Genome index file: genome.fa.fai
TE annotation file: TE.txt
Gene expression file: gene_expression.txt
TE expresion file: TE_expression.txt
Methylation CGmap files: leaf_01.CGmap.gz, leaf_02.CGmap.gz, root_01.CGmap.gz, root_02.CGmap.gz
```
For parameters like including unexpressed TEs, setting promoter regions, or plotting options, you can press `Enter` to use default values.
```
Include unexpressed TEs? (y/n, default n): n
Promoter upstream length from TSS (default 1500): 1500
Promoter downstream length from TSS (default 500): 500
Limit up-/down-stream range (bp)(default 15000): 15000
Tick size (bp)(default 5000): 5000
Window size (bp)(default 200): 400
Window number (default 156):  100
y-axis limit for gene expression vs TE/promoter mC plot (CG, default 50):  80
y-axis limit for gene expression vs TE/promoter mC plot (CHG, default 10):  40
y-axis limit for gene expression vs TE/promoter mC plot (CHH, default 10):  10
y-axis limit for TE expression vs TE mC plot (CH, default 15):  40
y-axis limit for TE expression vs TE mC plot (CG, default 30):  60
```

## Step 3: Pipeline modules description
| Module | Description                                            |
| ---- | ------------------------------------------------------ |
| 0    | Preprocessing of genome, genes, TEs, methylation files |
| 1    | TE distribution across genomic features                |
| 2    | Promoter-embedded TE families analysis                 |
| 3-1  | Gene-proximal TE methylation preprocessing             |
| 3-2  | Gene-proximal TE methylation plotting                  |
| 4    | Correlation between TEs and genes (single condition)   |
| 5    | Correlation between TEs and genes across conditions    |

## Step 4: Output files
After successful completion, outputs will include:
* Module 0:
    * `OUTPUT_0_embedded_TE_gene_number.txt`: Table showing numbers of promoter-embedded TEs and genes with TE embedded in promoters
    * `promoter.bed`
    * `TE_overlap_promoter.bed`: Names and location of promoter-embedded TEs and their neighbouring genes
    * `Tab_*.txt`: Tables showing CG/CHG/CHH methylation level of each TE, gene, promoter region under every condition
* Module 1:
    * `OUTPUT_1_TE_distribution_enrichment.png`
    * `OUTPUT_1_TE_distribution_percentage.png`
    
      <img width="602" height="261" alt="image" src="https://github.com/user-attachments/assets/261594e5-2f99-4ce6-a73a-f6b0dd55c24d" />

* Module 2:
    * `OUTPUT_2_Promoter_embedded_TE_family_enrichment.png`
    * `OUTPUT_2_Promoter_embedded_TE_family.txt`
      <img width="1557" height="420" alt="image" src="https://github.com/user-attachments/assets/0b1f2ee7-6b26-41c8-86fe-ab359850d527" />

* Module 3:
    * `OUTPUT_3_gene_TE_number.txt`: In every figure, it shows the upstream and downstream border (bp) of TE methylation impact, and the number of highly/lowly expressed genes
    * `OUTPUT_3_gene_proximal_TE_*.png`: Figures of TE CG/CHG/CHH methylation around highly/lowly expressed genes under every condition
       * Following are the examples using different methods to show the TE methylation pattern within 4 kb around top/bottom 10 % expressed genes:

         - (A) Using “average TE methylation within each window” (`average`) with following parameters: 100 bp window size, 100 bp window step, showing 95% CI, not showing the points of each TE site.
         
         - (B) Using “second-degree polynomial regression line” (`poly2`) with following parameters: showing 95% CI, showing the points of each TE site.
         
         - (C) Using “local regression line” (`loess `) with following parameters: span 0.02, showing 95% CI, showing the points of each TE site.
           <img width="1930" height="2243" alt="figure" src="https://github.com/user-attachments/assets/ee2d3a51-049e-474b-bdf3-7948a4ed5ccd" />



* Module 4:
    * `OUTPUT_4_geneexp/TEexp_*.png`: Line plots showing correlations between (1) TE methylation vs TE expression, (2) TE methylation, methylation of promoters with TEs/with no TEs vs gene expression, (3) TE expression vs gene expression under each stage. The TE-gene pairs number and the smoothing window size also showed in each figure.
    * `OUTPUT_4_correlation_*.png`: Bar plots showing pearson's or spearman's correlation coefficient between (1) TE CG/CHG/CHH methylation vs TE expression, (2) TE CG/CHG/CHH methylation vs gene expression, (3) TE expression vs gene expression under each stage
  <img width="1022" height="934" alt="image" src="https://github.com/user-attachments/assets/44ce1c5e-7d89-48b0-bc22-cfef639a2fef" />


* Module 5:
    * `OUTPUT_5_*_scatter.png`: Scatter plots showing correlations between changes of (1) TE CG/CHG/CHH methylation vs TE expression, (2) TE CG/CHG/CHH methylation vs gene expression, (3) TE expression vs gene expression comparing each two stages. The TE-gene pairs number and the regression line also showed in each figure.
    * `OUTPUT_5_Q2/4_boxplot_*.png`: Box plots showing the expression and methylation level changes of those negatively correlated gene-TE pairs under each two stages.
    * `OUTPUT_5_Q2/4_*.txt`: The name lists of those negatively correlated gene-TE pairs are also in the output.
  <img width="1138" height="722" alt="image" src="https://github.com/user-attachments/assets/68b8bce7-a32d-415e-b066-81f75e03aa21" />

