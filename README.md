# Glucocorticoid Response Gene That Modulates Cytokine Function in ASM Analysis

## Background
** This is a replication study of the 2014 paper [RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells (by Blanca E. Himes and Xiaofeng Jiang )](https://doi.org/10.1371/journal.pone.0099625) **

I have tried to modernize the pipeline of the paper and try to reproduce the results with modern tools.
I hope to verify the results of the study, and in the second part of the study, I will try to run the analysis on the most recent reference Human Genome HG38 (Instead of HG19)

## Dataset
The dataset used in the paper is available at [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=52778)
The dataset actually has 16 samples, combinations of untreated, treated with Albuterol and treated with Dexamethasone.
Since the study mostly focuses on the untreated and treated with dexamethasone samples, we will only use these two groups. 
The samples are listed in config/samples_table.tsv, four donors, two conditions.

## Two-Part Analysis Design
This study separates **pipeline effects** from **reference genome effects**:

| Analysis     | Reference | Pipeline         | Purpose                      |
|--------------|-----------|------------------|------------------------------|
| **Original** | hg19      | TopHat/Cufflinks | Baseline (2014)              |
| **Part 1**   | hg19      | STAR/DESeq2      | Isolate pipeline improvement |
| **Part 2**   | hg38      | STAR/DESeq2      | Full modernization           |

### Expected Insights
- **Part 1 vs Original**: Effect of modern aligner (STAR) and DE method (DESeq2)
- **Part 2 vs Part 1**: Effect of updated genome assembly and annotation
- **Part 2 vs Original**: Combined modernization effect

## Scripts
The scripts will be organized and enumerated under `scripts/`

### Environment setup
I am using the `environment.yml` file for dependency declaration environment setup with conda.
```bash
conda env create -f environment.yml
conda activate bioinf-grg
```

### Data Gathering automation
Data gathering is automated with `scripts/01_download_data.py` using. Simply running the script will download all the data required for the analysis.
Configuration for the data to be used is in `config/samples_table.tsv` and `config/configuration.py`.
Reference genomes are downloaded using `wget` from the datasets of [EMBL's European Bioinformatics Institute](https://ebi.ac.uk) and placed under `data/reference/{hg19|hg38}`.
Fasta files are downloaded using `fasterq-dump` from SRA and for compression using `pigz`. Fasta files are placed under `data/raw/`.
Git simply doesn't allow files of such a size to be uploaded, so they are ignored in `.gitignore`.