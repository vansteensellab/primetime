<p align="center"><img src="misc/logo.png" alt="primetime" width="80%"></p>

# Table of contents

- [**Introduction**](#introduction)
- [**Installation**](#installation)
- [**Setting up your configuration file**](#setting-up-your-configuration-file)
  - [**Setting up the input files**](#setting-up-the-input-files)
  - [**Setting up the conditions for the comparative analysis**](#setting-up-the-conditions-for-the-comparative-analysis)
  - [**Setting up the output directory**](#setting-up-the-output-directory)
  - [**Optional: Changing the p-value threshold**](#optional-changing-the-p-value-threshold)
  - [**Optional: Changing the library information**](#optional-changing-the-library-information)
- [**Output Files**](#output-files)

# Introduction

**primetime** is a pipeline designed for the analysis of TF prime reporter data. It processes fastq files, counts barcodes, clusters them, annotates them, and performs a differential TF activity analysis. The pipeline compares different conditions of the samples based on the 'condition' field in the configuration file.

# Installation

To install and run **primetime**, follow these steps:

1. Clone the repository:
    ```sh
    git clone https://github.com/vansteensellab/primetime.git
    cd primetime
    ```

2. Make sure you have snakemake installed. If you don't have it, you can install it with conda:
    ```sh
    conda install -c bioconda snakemake
    ```

We recommend trying to run primetime with our test data to check if everything was installed correctly.
This should run without any errors.

```sh
snakemake --configfile config.yaml --use-conda --cores 10 --printshellcmds
```

# Setting up your configuration file

Before running the **primetime**, you need to change the parameters on the `config.yaml` according to your data. Below, we will use our test data (that is meant to access the changes on TF activity on U2OS cells uppon calcitriol treatment) to guide you through the different sections of the configuration file and how to set them up:

> Note: The configuration file is sensitive to the number of spaces before each field. Make sure the format of your file is the same as the examples we provide.

## Setting up the input files

The `INPUT_DATA` parameter specifies the input data for the pipeline. Each sample should have a unique identifier and include information about whether it is pDNA, the condition, and the paths to the fastq files.

Example:
```yaml
INPUT_DATA:
  U2OS_DMSO_12W:
    is_pdna: False
    condition: DMSO
    fastq: 
      - test_data/U2OS_DMSO_1.fastq.gz
      - test_data/U2OS_DMSO_2.fastq.gz
      - test_data/U2OS_DMSO_3.fastq.gz

  U2OS_Calcitriol_12W:
    is_pdna: False
    condition: Calcitriol
    fastq: 
      - test_data/U2OS_calcitriol_1.fastq.gz
      - test_data/U2OS_calcitriol_2.fastq.gz
      - test_data/U2OS_calcitriol_3.fastq.gz

  pDNA:
    is_pdna: True
    condition: None
    fastq: 
      - test_data/pDNA_1.fastq.gz
```

## Setting up the conditions for the comparative analysis

The `COMPARATIVE_ANALYSIS` parameter defines the reference and contrast conditions for the comparative analysis. For our test data, we are comparing the cells treated with calcitriol against the control cells, treated with DMSO.

Example:
```yaml
COMPARATIVE_ANALYSIS:
  REFERENCE_CONDITION: DMSO
  CONTRAST_CONDITION: Calcitriol
```

## Setting up the output directory

The `OUTPUT_DIRECTORY` parameter specifies the directory where the output files will be saved.

Example:
```yaml
OUTPUT_DIRECTORY: test_data_output
```

## Optional: Changing the p-value threshold

The `PVALUE_THRESHOLD` parameter sets the p-value threshold for the differential analysis.

Example:
```yaml
PVALUE_THRESHOLD: 0.05
```

## Optional: Changing the library information

There are some additional parameters on the config file related to the barcodes used in the analysis

> Be aware that changing this parameters will impact the counting of the barcodes in the input reads. Only change this parameters if you know what you are doing.

Example:
```yaml
BARCODE_LENGTH: 12
BARCODE_DOWNSTREAM_SEQUENCE: CATCGTCGCATCCAAGAGGCTAGCTAACTA 
MAX_MISMATCH_DOWNSTREAM_SEQ: 3
BARCODE_ANNOTATION_FILE: misc/bc_annotation_prime.csv
```

# Output Files

**primetime** generates several output files during its execution:

1. **Quality Check Reports**:
    - `reads_QC/multiqc_report.html`: MultiQC report summarizing the quality of the reads.

2. **Barcode Counts**:
    - `barcode_counts/per_sample/{sample}_{replicate}.txt`: Raw barcode counts.
    - `barcode_counts/per_sample/{sample}_{replicate}.stats`: Statistics of barcode counts.
    - `barcode_counts/per_sample/{sample}_{replicate}.cluster.txt`: Clustered barcode counts.
    - `barcode_counts/per_sample/{sample}_{replicate}.cluster.annotated.txt`: Annotated barcode counts.

3. **Design File**:
    - `design.txt`: Design file containing sample, replicate, pDNA, and treatment information.

4. **Activity and QC Plots**:
    - `qc_plots/barcode_correlations.pdf`: QC plots showing barcode correlations.
    - `MPRAnalyze/cDNA_counts.txt`: cDNA counts.
    - `MPRAnalyze/pDNA_counts.txt`: pDNA counts.
    - `MPRAnalyze/activity.txt`: Activity counts.
    - `MPRAnalyze/corrected_activity.txt`: Corrected activity counts.

5. **MPRAnalyze Results**:
    - `results/primetime_results.txt`: Differential TF activity results.
    - `results/primetime_volcano.pdf`: Volcano plot of the differential activity results.

6. **Final Results**:
    - `results/primetime_lollipop.pdf`: Lollipop plot of the final results.

These files provide a comprehensive overview of the TF activity analysis and can be used for further downstream analysis.