# PeTEM 指令速查

## 互動式全流程
- `bash run_PeTEM.sh`  
  依提示選步驟 (0–5) 並輸入檔案路徑/參數。

## 安裝
- `bash setup.sh`  
  安裝所需系統工具、Python 套件、R 套件。

- `bash database.sh`  
  安裝常見gemo,e

## 個別步驟

python3 ../annotation_to_bed.py \
  --gtf ../genomes/arabidopsis/Arabidopsis_thaliana.TAIR10.60.gtf \
  --fai ../genomes/arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai \
  --promoter-up 1500 --promoter-down 500


### Step 0 — 前處理 (promoter + 甲基化表)
```bash
bash ../0_preprocessing.sh \
  -g annotation_bed/gene.bed -t ../data/TE.bed \
  -p annotation_bed/promoter.bed \
  -eg ../data/gene_expression.txt -et ../data/TE_expression.txt \
  -fai ../genomes/arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai \
  -m ../data/WT_01.CGmap.gz ../data/WT_02.CGmap.gz \
     ../data/drdd_01.CGmap.gz ../data/drdd_02.CGmap.gz 
```

<!-- ## 刪掉promoter計算 -->

<!-- ## TE.bed drdd_02.CGmap.gz 有修改過
sed 's/^Chr//' TE.bed > /work4/home/peiyu/tools/PeTEM/data/TE.bed
zcat drdd_01.CGmap.gz | sed 's/^Chr//' | gzip > /work4/home/peiyu/tools/PeTEM/data/drdd_01.CGmap.gz
zcat drdd_02.CGmap.gz | sed 's/^Chr//' | gzip > /work4/home/peiyu/tools/PeTEM/data/drdd_02.CGmap.gz
zcat WT_01.CGmap.gz | sed 's/^Chr//' | gzip > /work4/home/peiyu/tools/PeTEM/data/WT_01.CGmap.gz
zcat WT_02.CGmap.gz | sed 's/^Chr//' | gzip > /work4/home/peiyu/tools/PeTEM/data/WT_02.CGmap.gz -->


### Step 1 — TE 分布
```bash
bash ../1_TE_distribution.sh \
  -g annotation_bed/gene.bed -c annotation_bed/CDS.bed \
  -5 annotation_bed/UTR5.bed -e annotation_bed/exon.bed -3 annotation_bed/UTR3.bed \
  -p annotation_bed/promoter.bed -t ../data/TE.bed -fai ../genomes/arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai
```

## y不會超過100%

### Step 2 — 啟動子內 TE family 富集
```bash
Rscript ../2_TE_families.R \
  -a ../data/TE.bed \
  -i TE_overlap_promoter.bed \
  -T ../data/TE_family.txt
```
## 顏色之後要換掉

<!-- ### Step 3-1 — TE 影響距離前處理
```bash
bash ../3_1_TE_impact_distance_preprocess.sh \
  -m ../data/WT_01.CGmap.gz ../data/WT_02.CGmap.gz \
     ../data/drdd_01.CGmap.gz ../data/drdd_02.CGmap.gz
```

### Step 3-2 — TE 影響距離繪圖
```bash
bash 3_2_TE_impact_distance_plot.sh \
  -g PeTEM_data/gene.bed -t PeTEM_data/TE.bed \
  -eg PeTEM_data/gene_expression.txt -et PeTEM_data/TE_expression.txt \
  -lim 15000 -tick 5000 -WD 200 \
  -unexp n
``` -->

bash ../3_TE_impact_distance.sh \
  -m ../data/WT_01.CGmap.gz ../data/WT_02.CGmap.gz \
     ../data/drdd_01.CGmap.gz ../data/drdd_02.CGmap.gz  -g annotation_bed/gene.bed -t ../data/TE.bed \
  -eg ../data/gene_expression.txt -et ../data/TE_expression.txt \

## 合併兩個 ifpre3有了也不需要重跑
## 修改legned位置

### Step 4 — 單條件相關
```bash
Rscript ../4_single_condition_correlation.R \
  --eg ../data/gene_expression.txt --et ../data/TE_expression.txt 
```

# 改圖改成ggplot legend在外面

### Step 5 — 跨條件相關
```bash
Rscript ../5_cross_condition_correlation.R \
  --DEG ../data/DEG.txt --DETE ../data/DETE.txt 
```
