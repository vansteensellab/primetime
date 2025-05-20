---
title: Output files
nav_order: 4
---

# Output Files

**primetime** generates several output files during its execution. All the plots for QC are stored in the `primetime_QC` folder. Then the main result of **primetime** containing the activity and fold-change values is stored in the `primetime_results` folder, together with some additional plots.

Bellow we list all the QC and results files with a brief description.

## 1. Quality Check outputs

Inside the `primetime_QC` folder, several QC plots will be placed.

| File Name| Description|
|----------|------------|
| `barcode_correlations.pdf`              | Correlation of Log2(cDNA/pDNA) for different barcodes of each replicate.                                      |
| `replicate_correlations.pdf`            | Correlation of Log2(cDNA/pDNA) — after averaging the different barcodes — for each sample.                    |
| `bleedthrough_estimation.pdf`           | Estimation of the amount of pDNA bleedthrough (percentage of cDNA counts coming from pDNA) for each replicate.|
| `distribution_of_BC_counts.pdf`         | Distribution of the counts from all the barcodes of each replicate.                                           |
| `expected_vs_observed_pDNA_counts.pdf`  | Correlation of your pDNA counts (observed) with the ones of our lab (expected).                               |
| `read_counts.pdf`                       | Total amount of reads coming from each replicate.                                                             |

## 2. Main Results

| File Name| Description|
|----------|------------|
| `primetime_results.txt`         | Main result file containing adjusted p-values, fold-change values for each TF, and TF activity for each condition.        |
| `primetime_volcano.pdf`         | Volcano plot of the differential activity results.                                                                        |
| `primetime_lollipop.pdf`        | Lollipop plot showing the activity of each TF for both conditions, highlighting the differentially active ones.           |

{: .note }
Click here to read a deeper description on the expecter outcomes of **primetime**.

### Additional files
**primetime** also saves some additional files in the `tmp_primetime` folder, such as the barcode counts, and the results of the barcode clustering.

These files provide a comprehensive overview of the TF activity analysis and can be used for further downstream analysis.