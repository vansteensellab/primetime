################################################################################
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.2
# =============================================================================
# Description:
#
# This pipeline performs the analysis of the TF prime reporter data. It
# processes the fastq files, counts the barcodes, clusters them, annotates them,
# and finally performs a differential TF activity analysis.
# The diff. TF activity analysis compares the different conditions of the samples
# (based on the 'condition' field of the config file)
#
# Functionality:
#
# The pipeline is divided into the following steps:
# 1) Get barcode counts
# Knowing a constant sequence downstream of the barcode, the pipeline counts the
# number of times each barcode appears in the fastq files. Here, a filter of a
# maximum number of mismatches is applied to the downstream sequence.
# 2) Cluster barcodes
# The pipeline clusters the barcodes based on a distance of 1. This is done to
# deal with possible sequencing errors.
# 3) Annotate the barcodes
# The pipeline annotates the barcodes based on a reference file. This file
# is provided by us, and it contains the information about the barcodes.
# 4) Get activity and QC plots
# After annotating the barcodes, the pipeline calculates the activity of each
# TF, and generates some QC plots. Two important considerations about how the
# activity is calculated:
#   - First, both cDNA and pDNA are normalized by read depth (RPM units). Then, 
#     the cDNA RPM is corrected for DNA bleedthrough using the data of 
#     the random TFs. In the end, we calculate a corrected activity that is:
#     corrected_activity_RPM = (cDNA_RPM - bleedthrough_RPM) / pDNA_RPM.
#   - Second, the pipeline normalizes the corrected activity by the median of the
#     corrected activity of the random TFs, defining the 'normalized_activity'.
#     In other words, the final normalized activity tells how many times the 
#     activity of a TF is higher than the median of the random TFs.
# 6) Run MPRAnalyze
# The pipeline runs the MPRAnalyze package to perform a differential TF activity
# analysis on raw counts (not using the calculated activity. See the MPRAnalyze 
# paper for more details). Also, it important to know that MPRAnalyze will 
# compare the different conditions of the samples (based on the 'condition' field 
# of the config file).
# 7) Save final results
# Main final result is the primetime_results.txt file, which contains the list 
# of TFs with their corrected activity and p-values for the differential activity.
#
# =============================================================================
# Versions:
# 0.1 - Initial version
# 0.2 - Model MRPAnalyze with promoter information
# =============================================================================
__version__ = "0.2"

# Importing libraries
import sys
import inspect
import os

# Defining important output/input variables 
filename = inspect.getframeinfo(inspect.currentframe()).filename
base_dir = os.path.dirname(os.path.abspath(filename))
scripts_dir = os.path.join(base_dir, "scripts")
conda_envs_dir = os.path.join(base_dir, "conda_envs")
output_dir = config["OUTPUT_DIRECTORY"]
# List of paths to the dataframe of annotated barcodes (be used in the rule get_barcode_correlations)
path_unexpanded = os.path.join(output_dir, "barcode_counts/per_sample/{sample}_{replicate}.cluster.annotated.txt")
fastqc_unexpanded = os.path.join(output_dir, "reads_QC/{sample}_{replicate}.fastq.gz")
path_to_fastqc = []
path_to_annotated_dfs = []
for sample_name, info in config["INPUT_DATA"].items():
    for replicate_number in range(len(info['fastq'])):
        path_to_annotated_dfs.append(path_unexpanded.format(sample=sample_name, replicate=replicate_number + 1))
        path_to_fastqc.append(fastqc_unexpanded.format(sample=sample_name, replicate=replicate_number + 1))

################################################################################
### RULE ALL ###################################################################
################################################################################
rule all:
    input:
        #fastqc = os.path.join(output_dir, "reads_QC/multiqc_report.html"),
        mpra_analyze = os.path.join(output_dir, "results/primetime_results.txt"),
        # txt = os.path.join(output_dir, "MPRAnalyze/primetime_results.txt"),
        plots = os.path.join(output_dir, "results/primetime_lollipop.pdf"),
        qc = os.path.join(output_dir, "qc_plots/barcode_correlations.pdf"),

################################################################################
### HELPER FUNCTIONS ###########################################################
################################################################################

# Function to return the sample name
def get_samples_data(config, wildcards):
    """
    This function returns the sample name based on the wildcards.
    defines the wildcards:
    - sample: sample name
    - replicate: replicate number
    ---
    Used in the rule get_barcode_counts
    """
    # Subtract 1 because list indexing starts at 0, but replicate numbers start at 1.
    # this makes the 'replicate' as wildcard 1-based
    index = int(wildcards.replicate) - 1
    return config["INPUT_DATA"][wildcards.sample]["fastq"][index]


def print_config(key, value, space=0):
    """
    Just a helper function to print the config file in a nice way :)
    """
    # Check if value is a dictionary
    if isinstance(value, dict):
        print(f"\033[1m{key}:\033[0m")
        for k, v in value.items():
            if isinstance(v, dict):
                print()
                print_config(k, v, space + 2)
            else:
                print(f'{" " * space} - {k}: ', end="")
                print(v)
    else:
        print(f"\033[1m{key}:\033[0m {value}")

# print("=" * 80)
# print("PrimeTime: TF reporter pipeline")
# print("Version", __version__)
# print("-" * 80)
# for key, value in config.items():
#     print_config(key, value)
# print("=" * 80)

################################################################################
### RULES ######################################################################
################################################################################
# 0) Get quality check of the reads
rule prepare_to_fastqc:
    input:
        lambda wildcards: get_samples_data(config, wildcards),
    output:
        os.path.join(output_dir, "reads_QC/{sample}_{replicate}.fastq.gz")
    conda:
        os.path.join(conda_envs_dir, "fastqc.yaml")
    shell:
        """
        # Create symbolic link to the fastq file with the new name
        ln -s {input} {output}
        """

rule fastqc:
    input: 
        files=path_to_fastqc,
        multiqc_config = os.path.join(output_dir, "reads_QC/multiqc_config.yaml")

    output: os.path.join(output_dir, "reads_QC/multiqc_report.html")
    conda:
        os.path.join(conda_envs_dir, "fastqc.yaml")
    shell:
        """
        fastqc -o {output_dir}/reads_QC {input.files}
        multiqc {output_dir}/reads_QC -o {output_dir}/reads_QC -c {input.multiqc_config}
        """

# 1) Get barcode counts
rule get_barcode_counts:
    input:
        lambda wildcards: get_samples_data(config, wildcards),
    output:
        txt=os.path.join(output_dir, "barcode_counts/per_sample/{sample}_{replicate}.txt"),
        stats=os.path.join(output_dir, "barcode_counts/per_sample/{sample}_{replicate}.stats"),
    params:
        script=os.path.join(scripts_dir, "get_barcode_counts.py"),
        bc_length=config["BARCODE_LENGTH"],
        bc_downstream_seq=config["BARCODE_DOWNSTREAM_SEQUENCE"],
        max_mismatch=config["MAX_MISMATCH_DOWNSTREAM_SEQ"],
    threads: 20
    conda:
        os.path.join(conda_envs_dir, "test_bc_counts.yaml")
    shell:
        """
        python {params.script} \
--fastq {input} \
--bc_length {params.bc_length} \
--bc_downstream_seq {params.bc_downstream_seq} \
--max_mismatch {params.max_mismatch} \
--threads {threads} \
> {output.txt} 2> {output.stats}
        """

# 2) Cluster barcodes
rule cluster_barcodes:
    input:
        os.path.join(output_dir, "barcode_counts/per_sample/{sample}_{replicate}.txt"),
    output:
        os.path.join(output_dir, "barcode_counts/per_sample/{sample}_{replicate}.cluster.txt"),
    conda:
        os.path.join(conda_envs_dir, "starcode.yaml")
    threads: 10
    shell:
        """
        starcode --threads {threads} --print-clusters -i {input} --dist 1 2> /dev/null | \
sort -k1,1 > {output}
        """

# 3) Annotate the barcodes
rule annotate_barcodes:
    input:
        os.path.join(output_dir, "barcode_counts/per_sample/{sample}_{replicate}.cluster.txt"),
    output:
        os.path.join(output_dir, "barcode_counts/per_sample/{sample}_{replicate}.cluster.annotated.txt"),
    params:
        script=os.path.join(scripts_dir, "annotate_barcodes.R"),
        bc_annotation=config['BARCODE_ANNOTATION_FILE'],
        is_pdna=lambda wildcards: config["INPUT_DATA"][wildcards.sample]["is_pdna"],
        treatment=lambda wildcards: config["INPUT_DATA"][wildcards.sample]["condition"],
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        Rscript {params.script} \
--bc_counts {input} \
--bc_annotation {params.bc_annotation} \
--sample {wildcards.sample} \
--replicate_number {wildcards.replicate} \
--is_pdna {params.is_pdna} \
--treatment {params.treatment} \
--output {output}
        """

# Aux rule to get the design dataframe: sample, replicate, pDNA?, treatment
rule get_design_df:
    output:
        os.path.join(output_dir, "design.txt")
    run:
        samples_dict = config["INPUT_DATA"]
        reference_condition = config["COMPARATIVE_ANALYSIS"]["REFERENCE_CONDITION"]
        contrast_condition = config["COMPARATIVE_ANALYSIS"]["CONTRAST_CONDITION"]
        with open(output[0], "w") as f:
            f.write("sample\treplicate\tpDNA\ttreatment\n")
            for sample_name, info in config["INPUT_DATA"].items():
                if info['condition'] == reference_condition:
                    is_pdna = info["is_pdna"]
                    for replicate_number in range(len(info['fastq'])):
                        replicate = f"{sample_name}_{replicate_number + 1}"
                        treatment = info["condition"]
                        f.write(f"{sample_name}\t{replicate}\t{is_pdna}\t{treatment}\n")
            for sample_name, info in config["INPUT_DATA"].items():
                if info['condition'] == contrast_condition:
                    is_pdna = info["is_pdna"]
                    for replicate_number in range(len(info['fastq'])):
                        replicate = f"{sample_name}_{replicate_number + 1}"
                        treatment = info["condition"]
                        f.write(f"{sample_name}\t{replicate}\t{is_pdna}\t{treatment}\n")
            for sample_name, info in config["INPUT_DATA"].items():
                if info['condition'] != reference_condition and info['condition'] != contrast_condition:
                    is_pdna = info["is_pdna"]
                    for replicate_number in range(len(info['fastq'])):
                        replicate = f"{sample_name}_{replicate_number + 1}"
                        treatment = info["condition"]
                        f.write(f"{sample_name}\t{replicate}\t{is_pdna}\t{treatment}\n")

# 4) Get activity and QC plots
rule get_activity_and_qc_plots:
    input:
        dfs = path_to_annotated_dfs,
        design=os.path.join(output_dir, "design.txt")
    output:
        plots = os.path.join(output_dir, "qc_plots/barcode_correlations.pdf"),
        cDNA = os.path.join(output_dir, "MPRAnalyze/cDNA_counts.txt"),
        pDNA = os.path.join(output_dir, "MPRAnalyze/pDNA_counts.txt"),
        activity = os.path.join(output_dir, "MPRAnalyze/activity.txt"),
        corrected_activity = os.path.join(output_dir, "MPRAnalyze/corrected_activity.txt")
    params:
        script=os.path.join(scripts_dir, "get_activity_and_qc_plots.R"),
        df_basedir=os.path.join(output_dir, "barcode_counts/per_sample"),
        plots_basedir=os.path.join(output_dir, "qc_plots/"),
        MPRAnalyze_basedir = os.path.join(output_dir, "MPRAnalyze")
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        mkdir -p {params.plots_basedir}
        mkdir -p {params.MPRAnalyze_basedir}

        Rscript {params.script} \
        --df_basedir {params.df_basedir} \
        --plots_basedir {params.plots_basedir} \
        --mpranalyze_basedir {params.MPRAnalyze_basedir} \
        --design {input.design}
        """

# 5) Run MPRAnalyze
rule run_MPRAnalyze:
    input:
        cDNA = os.path.join(output_dir, "MPRAnalyze/cDNA_counts.txt"),
        pDNA = os.path.join(output_dir, "MPRAnalyze/pDNA_counts.txt"),
        design = os.path.join(output_dir, "design.txt"),
        corrected_activity = os.path.join(output_dir, "MPRAnalyze/corrected_activity.txt")
    output:
        txt = os.path.join(output_dir, "results/primetime_results.txt"),
        plots = os.path.join(output_dir, "results/primetime_volcano.pdf")
    params:
        script=os.path.join(scripts_dir, "run_comparative_analysis.R"),
        out_basedir=os.path.join(output_dir, "MPRAnalyze"),
        pval_threshold=config["PVALUE_THRESHOLD"],
        plot_output_dir=os.path.join(output_dir, "results"),
        contrast_condition=config["COMPARATIVE_ANALYSIS"]["CONTRAST_CONDITION"],
        reference_condition=config["COMPARATIVE_ANALYSIS"]["REFERENCE_CONDITION"]
    conda:
        os.path.join(conda_envs_dir, "mpranalyze.yaml")
    threads: 100
    shell:
        """
        mkdir -p {params.plot_output_dir}
        # MPRAnalyze with fit.se
        Rscript {params.script} \
        --cdna {input.cDNA} \
        --pdna {input.pDNA} \
        --design {input.design} \
        --output {params.out_basedir} \
        --pval_threshold {params.pval_threshold} \
        --contrast_condition {params.contrast_condition} \
        --reference_condition {params.reference_condition} \
        --corrected_activity {input.corrected_activity} \
        --plot_output {params.plot_output_dir}
        """

# 6) Merge activity and MPRAnalyze results, and save final results
rule get_final_results:
    input: 
        corrected_activity = os.path.join(output_dir, "MPRAnalyze/activity.txt"), ###### Note that i changed to the corrected
        mpra_results = os.path.join(output_dir, "results/primetime_results.txt")
    output:
        # txt = os.path.join(output_dir, "MPRAnalyze/primetime_results.txt"),
        plots = os.path.join(output_dir, "results/primetime_lollipop.pdf")
    params:
        script=os.path.join(scripts_dir, "get_final_results_mpranalyze.R"),
        output_dir = os.path.join(output_dir, "MPRAnalyze"),
        plot_output_dir = os.path.join(output_dir, "results"),
        design = os.path.join(output_dir, "design.txt"),
        contrast_condition=config["COMPARATIVE_ANALYSIS"]["CONTRAST_CONDITION"],
        reference_condition=config["COMPARATIVE_ANALYSIS"]["REFERENCE_CONDITION"]
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        Rscript {params.script} \
        --corrected_activity {input.corrected_activity} \
        --mpra_results {input.mpra_results} \
        --output {params.output_dir} \
        --design {params.design} \
        --contrast_condition {params.contrast_condition} \
        --reference_condition {params.reference_condition} \
        --plot_output {params.plot_output_dir}
        """