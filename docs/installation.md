---
title: Installation
nav_order: 2
---

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