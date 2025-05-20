---
title: Setting up your configuration file
nav_order: 3
---

# Setting up your configuration file
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

Before running the **primetime**, you need to change the parameters on the `config.yaml` according to your data. Below, we will use our test data (that is meant to access the changes on TF activity on U2OS cells uppon calcitriol treatment) to guide you through the different sections of the configuration file and how to set them up:


{: .note }
The configuration file is sensitive to the number of spaces before each field. Make sure the format of your file is the same as the examples we provide.

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

{: .note }
Be aware that changing this parameters will impact the counting of the barcodes in the input reads. Only change this parameters if you know what you are doing.

Example:
```yaml
BARCODE_LENGTH: 12
BARCODE_DOWNSTREAM_SEQUENCE: CATCGTCGCATCCAAGAGGCTAGCTAACTA 
MAX_MISMATCH_DOWNSTREAM_SEQ: 3
BARCODE_ANNOTATION_FILE: misc/bc_annotation_prime.csv
EXPECTED_PDNA_COUNTS: misc/expected_pDNA_counts.txt
```