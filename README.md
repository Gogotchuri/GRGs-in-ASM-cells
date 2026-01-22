# Glucocorticoid-Responsive Genes in Airway Smooth Muscle Cells: A Comparison of Modern and Legacy RNA-seq Pipelines

*This is a replication study of the 2014 original paper [RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells (by Blanca E. Himes and Xiaofeng Jiang )](https://doi.org/10.1371/journal.pone.0099625)*

## Abstract
RNA-seq analysis pipelines have evolved significantly since 2014, with improvements in alignment algorithms, quantification methods and analysis tools. Using a newer pipeline and reference genome,
We reanalyzed the dataset [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) - RNA-seq data from eight airway smooth muscle cells (ASM) from four donors.
The dataset was used in the original paper to characterize transcriptomic changes in human ASM cells treated with dexamethasone (potent glucocorticoid).
The original paper used the Tuxedo tools, TopHat2 for alignment and Cufflinks+Cuffdiff for quantification and differential expression analysis with reference genome **hg19**.<br/>
In this study, we modernized the pipeline using STAR, featureCounts and DESeq2, and also compared the reference genome effect, by running the pipeline on both the original reference genome **hg19** 
and the latest reference genome available **hg38**. The STAR alignment achieved a 97.6% mapping rate compared to 83.4% in the original paper. Junction-spanning reads were also improved from 26% to 45%. 
Despite the technical metric improvements, we have re-discovered and validated 96.2% (304/316) (with **hg18**) of the genes identified by the original paper with **0.998** Pearson correlation between the log2 fold changes.
The modern pipeline with the new reference genome also identified a total of 910 significant differentially expressed genes (padj < 0.05 and |log2FoldChange| > 1), reflecting higher sensitivity.
Although the reference genome effect was not significant, it might still be interesting for biological interpretation.<br/>
We conclude that the original well-designed RNA-seq analysis was robust to pipeline variations and would likely still benefit from the modern tooling.

## 1. Introduction

### 1.1 Pipeline Background
RNA-seq is a widely used method to study transcriptomes in cells. The tools for alignment evolved significantly. One of the most effective tools for alignment we have today - STAR,
was first introduced in a paper published at the end of 2012 ([PMC3530905](https://doi.org/10.1093/bioinformatics/bts635)).<br/>
In 2013 a paper was published (["Systematic evaluation of spliced alignment programs for RNA-seq data"](https://www.nature.com/articles/nmeth.2722)) that compared the performance of several alignment tools.<br>
GSNAP, GSTRUCT, MapSplice, and STAR emerged as the top performers. But the interesting thing about STAR was its performance, it performed about 180x faster than GSNAP and MapSplice while maintaining comparable accuracy and in some cases improved junction-spanning with two passes.
This allowed STAR to be run on a single modern computer and take a few hours to perform alignment instead of days, making discovery and iteration faster, cheaper and overall more accessible. Running the scripts in this repository took around 2 hours in total. 
For assigning sequence reads to genomic features, the pipeline uses [featureCounts (published in 2013 paper)](https://doi.org/10.1093/bioinformatics/btt656) which is another highly efficient tool, order of magnitude faster than comparable tools like [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/).
For fold change and dispersion analysis, the pipeline uses a python version of [DESeq2 (published in 2014 paper)](https://doi.org/10.1186/s13059-014-0550-8).

### 1.2 Dataset
The dataset used in the paper is available at [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)
The dataset actually has 16 samples, combinations of untreated, treated with Albuterol and treated with Dexamethasone.
Since the study mostly focuses on the untreated and treated with dexamethasone (12 hours) samples, we will only use these two groups. 
The samples are listed in config/samples_table.tsv, four donors, two conditions.
//embed: config/samples_table.tsv

### 1.3 Original Study Summary
The original paper used the Tuxedo tools, TopHat2 for alignment and Cufflinks+Cuffdiff for quantification and differential expression analysis with reference genome **hg19**.</br>
Analysis in the original paper found 316 genes that were significantly differentially expressed between the two conditions with p < 0.05.</br>
The **CRISPLD2** gene was identified as the novel finding by the paper. Here are the top gene findings from the original paper (padj < 1E-16 is marked as 0):

| **Gene** | **log2 FC** | **P adj (Q)** |
|----------|-------------|---------------|
| C7       | -3.35       | 0             |
| CCDC69   | -2.92       | 0             |
| DUSP1    | -2.99       | 0             |
| FKBP5    | -3.95       | 0             |
| GPX3     | -3.76       | 0             |
| KLF15    | -4.58       | 0             |
| MAOA     | -3.29       | 0             |
| SAMHD1   | -3.83       | 0             |
| SERPINA3 | -3.34       | 0             |
| SPARCL1  | -4.70       | 0             |
| C13orf15 | -3.27       | 2.5E-13       |
| TSC22D3  | -3.27       | 2.5E-13       |
| CRISPLD2 | -2.70       | 6.9E-13       |


### 1.4 Study Design and Objectives
*You can find the original analysis description in the original paper, on page 9, paragraph "RNA-Seq Data Analysis."*

This study tries to separate **pipeline effects** from **reference genome effects** by performing two part analyses:

| Analysis     | Reference | Pipeline                    | Purpose                      |
|--------------|-----------|-----------------------------|------------------------------|
| **Original** | hg19      | TopHat/Cufflinks + Cuffdiff | Baseline (2014)              |
| **Part 1**   | hg19      | STAR/featureCount + DESeq2  | Isolate pipeline improvement |
| **Part 2**   | hg38      | STAR/featureCount + DESeq2  | Full modernization           |

We will try to validate the findings from the original paper and compare them to the findings from the modern pipeline.

#### 1.4.1 Expected Insights
- **Part 1 vs Original**: Effect of the updated pipeline
- **Part 2 vs Part 1**: Effect of updated genome assembly and annotation
- **Part 2 vs Original**: Combined modernization effect

## 2. Results
### 2.1 Alignment Quality
Modern STAR alignment achieved a 97.6% mapping rate compared to 83.4% in the original paper using TopHat2. Junction-spanning reads were also improved from 26% to 45%.
Indicating improvement in splice site detection by STAR. Although the read counts are lower


| Metric                 | hg19 (TopHat2) | hg19 (STAR+featureCount) | hg38 (STAR+featureCount) |
|------------------------|----------------|--------------------------|--------------------------|
| Reads (avg)            | 58,9M          | 48,9M                    | 48,9M                    |
| Reads (min)            | 44,2M          | 33,7M                    | 33,7M                    |
| Reads (max)            | 71,3M          | 68,6M                    | 68,6M                    |
| Mapped (avg)           | 83.36%         | 97.59%                   | 97.64%                   |
| Mapped (min)           | 81.94%         | 96.26%                   | 96.33%                   |
| Mapped (max)           | 84.34%         | 98.10%                   | 98.13%                   |
| Junctions Mapped (avg) | 26.43%         | 45.57%                   | 44.52%                   |
| Junctions Mapped (min) | NA             | 43.06%                   | 42.07%                   |
| Junctions Mapped (max) | NA             | 47.10%                   | 45.94%                   |

## Scripts
The scripts will be organized and enumerated under `scripts/`. Those should be run in the order they are listed, which will run the complete analysis used in this study.
and will also produce the figures and tables used in the analysis.

### Environment setup
I am using the `environment.yml` file for dependency declaration environment setup with conda.
```bash
conda env create -f environment.yml
conda activate bioinf-grg
```

### 01 - Data Gathering automation
Data gathering is automated with `scripts/01_download_data.py` using. Simply running the script will download all the data required for the analysis.
Configuration for the data to be used is in `config/samples_table.tsv` and `config/configuration.py`.<br>
Reference genomes are downloaded using `wget` from the datasets of [EMBL's European Bioinformatics Institute](https://ebi.ac.uk) and placed under `data/reference/{hg19|hg38}`.
Fasta files are downloaded using `fasterq-dump` from SRA and for compression using `pigz`. Fasta files are placed under `data/raw/`.<br/>
Git simply doesn't allow files of such a size to be uploaded, so they are ignored in `.gitignore`.</br>
Additionally, I have manually downloaded and formatted tables from the original paper supplementary materials. Also, I have created a new table for comparison, with stats described in the results section.
The 316 genes identified by the original paper are from the supplementary [table S3](https://doi.org/10.1371/journal.pone.0099625.s014)
I modified the headers to correspond to the correct columns produced by the newer pipeline. All of the aforementioned files are placed under `config/`.

### 02 - Alignment
Part A: STAR index is built using `scripts/02a_build_star_index.py`. This script takes around an hour to run for the two genomes and needs >= 32GB of RAM.
Part B: Alignment is performed using `scripts/02b_align_reads.py`. This uses STAR index built in Part A and takes around 5 minutes per sample, producing BAM files for each sample.
After this part, every next step is significantly faster.

### 03 - Gene Quantification
*featureCounts* uses the produced BAM files to quantify gene expression. featureCounts is run using `scripts/03_quantify_genes.py`. 
The script contains post-processing and cleaning steps to produce the final table under `results/hg{19,38}/tables/gene_counts.tsv` containing gene counts for each sample.

### 04 - Gene Differential Expression Analysis
Differential expression analysis is performed using `scripts/04_analysis_deseq2.py`. This script uses the gene counts table produced in Step 3 and performs DESeq2 analysis.
The complete list of results with padj < 0.05 is saved under `results/hg{19,38}/tables/deseq2_results.csv`. The significant genes with |log2FoldChange| > 1 are saved under `results/hg{19,38}/tables/top_results_from_deseq2.csv`.
At this point, the results are ready to be processed and visualized for presentation. We are technically done.

### 05 - Comparisons
The results from Part 1 and Part 2 are compared to the original findings and to each other. Respectively by scripts `scripts/05a_comparisons_to_original.py` and `scripts/05b_comparisons_between_genomes`. 
The first script produces tables of validated and discarded genes for each genome from the original findings. And the second script produces tables to evaluate the genome effect on the results. Produced files are placed under `results/general/tables/`.

### 06 - General Stats
This step produces `general_stats_summary.tsv` under each genome tables folder. The stats should be comparable to the original paper stats. The script used is `scripts/06_collect_general_stats.py`.

### 07 - Figures
The figures are produced using `scripts/07_plots.py`. This script uses the tables produced in the previous steps and plots the figures under `results/figures/`.
