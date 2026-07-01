# Tutorial
> This tutorial demonstrates how to run the PeTEM pipeline using the provided example dataset. (⏳Runtime: 10 minutes) <br>
> For detailed information on installation, input formats, parameters, and outputs, please refer to the [README](https://github.com/PaoyangLab/PeTEM/tree/main#promoter-embedded-te-methylation-petem-analyzer).


**Tested Environment**
| Component | Specification |
|------------|------------|
| Operating System | Ubuntu 20.04.2 LTS (Focal Fossa) |
| Linux Kernel | 5.15.0-138-generic |
| CPU | Intel Xeon Gold 6238 @ 2.10 GHz |
| Memory | 1.4 TB RAM |
| R | 4.2.2 |
| Python | 3.12.3 |
| Conda | 24.11.3 |

## Table of Contents

**Setup**
1. [Clone repository](#clone-repository) 
2. [Download example data](#download-example-data) 
3. [Set up environment](#set-up-environment) 

**Pipeline Workflow** (Modules 1–5 can be run independently after Module 0)
- [Module 0: Preprocessing](#module-0-preprocessing) 
- [Module 1: TE Distribution](#module-1-te-distribution) 
- [Module 2: Enriched promoter-embedded TE Families](#module-2-enriched-promoter-embedded-te-families) 
- [Module 3: TE methylation near gene](#module-3-te-methylation-near-gene) 
- [Module 4: Correlation between methylation and expression](#module-4-correlation-between-methylation-and-expression) 
- [Module 5: Associated TE and gene pairs](#module-5-associated-te-and-gene-pairs) 

## Clone repository
```bash
git clone https://github.com/PaoyangLab/PeTEM.git
cd PeTEM
```

## Download example data
```bash
wget https://github.com/PaoyangLab/PeTEM/releases/download/Example_data/PeTEM_data.tar.gz
tar -xzvf PeTEM_data.tar.gz
```

> This command downloads and extracts the `PeTEM_data` directory, which contains 9 files.

| File                                                                           | Description                 |
| ------------------------------------------------------------------------------ | --------------------------- |
| `genome.fa.fai`                                                                | Genome index                |
| `genomic.gff`                                                                  | Genome annotation           |
| `TE.txt`                                                                       | TE annotation               |
| `gene_expression.txt`                                                          | Gene expression             |
| `TE_expression.txt`                                                            | TE expression               |
| `leaf_01.CGmap.gz`, `leaf_02.CGmap.gz`, `root_01.CGmap.gz`, `root_02.CGmap.gz` | Methylation CGmap           |

## Set up environment 
> This tutorial uses the Conda setup as an example. Alternative setup methods, including **Docker** and **Local setup**, are described in the [README](https://github.com/PaoyangLab/PeTEM#set-up-environment).

```bash
conda env create -f environment.yml
conda activate petem
bash ./scripts/env_check.sh ##optional
```

## Module 0: Preprocessing
* __Command:__

  ```
  ## Run module 0
  ./petem --0 \
    -g ./PeTEM_data/genomic.gff \
    -t ./PeTEM_data/TE.txt \
    -eg ./PeTEM_data/gene_expression.txt \
    -et ./PeTEM_data/TE_expression.txt \
    -f ./PeTEM_data/genome.fa.fai \
    -m ./PeTEM_data/leaf_01.CGmap.gz \
       ./PeTEM_data/leaf_02.CGmap.gz \
       ./PeTEM_data/root_01.CGmap.gz \
       ./PeTEM_data/root_02.CGmap.gz \
    -up 1500 \
    -dn 500 \
    -o ./example/module0 
  ```

* __Outputs:__

  | File | Description |
  |---|---|
  | `OUTPUT_0_embedded_TE_gene_number.txt` | Summary table showing the numbers of promoter-embedded TEs and genes containing promoter-embedded TEs. |
  | `promoter.bed` | BED file containing promoter coordinates generated from the genome annotation. |
  | `TE_overlap_promoter.bed` | Genomic locations of promoter-embedded TEs and their associated neighbouring genes. |
  | `Tab_*.txt` | Tables containing CG, CHG, and CHH methylation levels for each TE, gene, and promoter region across all conditions. |

## Module 1: TE Distribution
* __Command:__

  ```
  ## Run module 1
  ./petem --1 \
    --manifest ./example/module0/PETEM_MODULE0_MANIFEST.json \
    -o ./example/module1
  ```

* __Outputs:__

  | File | Description |
  |---|---|
  | `OUTPUT_1_TE_distribution_enrichment.png` | Bar plot showing the enrichment of transposable elements across different genomic features. |
  | `OUTPUT_1_TE_distribution_percentage.png` | Bar plot showing the percentage distribution of transposable elements across different genomic features. |

* __Figures:__
  
  <img width="1127" height="572" alt="image" src="https://github.com/user-attachments/assets/0dc997ca-8baa-4961-b7c5-5a146c10a8ec" />



## Module 2: Enriched promoter-embedded TE Families
* __Command:__

  ```
  ## Run module 2
  ./petem --2 \
    --manifest ./example/module0/PETEM_MODULE0_MANIFEST.json \
    -o ./example/module2
  ```

* __Outputs:__

  | File | Description |
  |---|---|
  | `OUTPUT_2_Promoter_embedded_TE_family_enrichment.png` | Bar plot showing the enrichment of TE families among promoter-embedded transposable elements. |
  | `OUTPUT_2_Promoter_embedded_TE_family.txt` | Table containing enrichment statistics and proportions of promoter-embedded TE families. |

* __Figure:__

  <img width="1127" height="525" alt="image" src="https://github.com/user-attachments/assets/e876d0fd-93c4-4715-9d09-88f441f5f433" />


## Module 3: TE methylation near gene
* __Command:__

  ```
  ## Run module 3
  ./petem --3 \
    --manifest ./example/module0/PETEM_MODULE0_MANIFEST.json \
    -o ./example/module3 
  ```

* __Outputs:__

  | File | Description |
  |---|---|
  | `OUTPUT_3_genes_number.txt` | Summary table reporting the numbers of highly and lowly expressed genes included in each analysis. |
  | `OUTPUT_3_TE_methylation_near_genes_*.png` | Plots showing CG, CHG, and CHH methylation patterns of gene-proximal TEs around highly and lowly expressed genes under each condition. Different visualization methods can be selected to display TE methylation trends across genomic regions. |
  | `OUTPUT_3_nonTE_methylation_near_genes_*.png` | Plots showing CG, CHG, and CHH methylation patterns of non-TE regions around highly and lowly expressed genes under each condition. |
   
* __Figures:__

  <img width="766" height="917" alt="image" src="https://github.com/user-attachments/assets/7184c2a6-4b07-4980-b55c-221fd9af020b" />







## Module 4. Correlation between methylation and expression
* __Command:__

  ```
  ## Run module 4
  ./petem --4 \
    --manifest ./example/module0/PETEM_MODULE0_MANIFEST.json \
    -o ./example/module4
  ```

* __Outputs:__

  | File | Description |
  |---|---|
  | `OUTPUT_4_geneexp/TEexp_*.png` | Line plots showing the relationships between (1) TE methylation and TE expression, (2) TE methylation, promoter methylation (with or without embedded TEs), and gene expression, and (3) TE expression and gene expression across conditions. Each figure reports the number of TE–gene pairs analyzed and the smoothing window size used. |
  | `OUTPUT_4_correlation_*.png` | Bar plots summarizing Pearson or Spearman correlation coefficients for (1) TE CG/CHG/CHH methylation vs. TE expression, (2) TE CG/CHG/CHH methylation vs. gene expression, and (3) TE expression vs. gene expression across conditions. |

* __Figures:__

  <img width="1127" height="1650" alt="image" src="https://github.com/user-attachments/assets/b770d81d-e07b-462b-9030-6e2532f1078b" />


## Module 5. Associated TE and gene pairs
* __Command:__

  ```
  ## Run module 5
  ./petem --5 \
    --manifest ./example/module0/PETEM_MODULE0_MANIFEST.json \
    -stage root leaf \
    -o ./example/module5
  
  ```

* __Outputs:__

  | File | Description |
  |---|---|
  | `OUTPUT_5_*_scatter.png` | Scatter plots showing the relationships between changes in (1) TE CG/CHG/CHH methylation and TE expression, (2) TE CG/CHG/CHH methylation and gene expression, and (3) TE expression and gene expression between two conditions. Each figure includes the number of TE–gene pairs analyzed and the fitted regression line. |
  | `OUTPUT_5_Q2/4_boxplot_*.png` | Box plots showing expression and methylation changes of negatively correlated TE–gene pairs identified between two conditions. |
  | `OUTPUT_5_Q2/4_*.txt` | Tables containing the lists of negatively correlated TE–gene pairs identified in the analysis. |

* __Figures:__

  <img width="1127" height="1407" alt="image" src="https://github.com/user-attachments/assets/60694bdd-73cc-4e0e-9e8a-23a29ad2acbf" />


* __Tables:__

  `OUTPUT_5_Q2_gene_TE_pairs_CG_root_leaf.txt`
  | gene_id | TE_id |
  |----------|----------|
  | AT1G63055 | AT1TE77010 |
  | AT1G63057 | AT1TE77010 |
  | AT3G47965 | AT3TE71740 |
  | AT5G28660 | AT5TE38885 |
  
  `OUTPUT_5_Q4_gene_TE_pairs_CG_root_leaf.txt`
  | gene_id | TE_id |
  |----------|----------|
  | AT1G07135 | AT1TE07165 |
  | AT1G32880 | AT1TE38580 |
  | AT1G61340 | AT1TE74615 |
  | AT3G25820 | AT3TE39395 |
  | AT3G27710 | AT3TE42655 |
  | AT3G44260 | AT3TE64560 |
  | AT3G50800 | AT3TE76570 |
  | AT3G58860 | AT3TE88640 |
  | AT4G04840 | AT4TE11400 |
  | AT5G47440 | AT5TE69240 |
