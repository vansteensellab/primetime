################################################################################
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 1.0.0
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
# TF, and generates some QC plots.
# 6) Run Comparative analysis
# The pipeline runs the BCalm package to perform a differential TF activity
# analysis on raw counts. See the BCalm paper for more details). Also, it is
# important to know that BCalm will compare the different conditions of the
# samples (based on the 'condition' field of the config file).
# 7) Save final results
# Main final result is the primetime_results.txt file, which contains the list
# of TFs with their corrected activity and p-values for the differential activity.
# =============================================================================
__version__ = "1.0.0"

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
path_unexpanded = os.path.join(
    output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.annotated.txt"
)
fastqc_unexpanded = os.path.join(output_dir, "reads_QC/{sample}_{replicate}.fastq.gz")
path_to_fastqc = []
path_to_annotated_dfs = []
for sample_name, info in config["INPUT_DATA"].items():
    for replicate_number in range(len(info["fastq"])):
        path_to_annotated_dfs.append(
            path_unexpanded.format(sample=sample_name, replicate=replicate_number + 1)
        )
        path_to_fastqc.append(
            fastqc_unexpanded.format(sample=sample_name, replicate=replicate_number + 1)
        )


################################################################################
### RULE ALL ###################################################################
################################################################################
rule all:
    input:
        results=os.path.join(output_dir, "primetime_results/primetime_results.txt"),
        loli=os.path.join(output_dir, "primetime_results/primetime_lollipop.pdf"),
        volcano=os.path.join(output_dir, "primetime_results/primetime_volcano.pdf"),
        qc=os.path.join(output_dir, "primetime_QC/bleedthrough_estimation.pdf"),
        replicate_corr=os.path.join(
            output_dir, "primetime_QC/replicate_correlations.pdf"
        ),
        expected_vs_obs=os.path.join(
            output_dir, "primetime_QC/expected_vs_observed_pDNA_counts.pdf"
        ),
        bc_counts_for_top_TFs=os.path.join(
            output_dir, "primetime_results/BC_counts_for_top_TFs.pdf"
        ),


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
        print(f"\033[1m{' ' * space}{key}:\033[0m")
        for k, v in value.items():
            if isinstance(v, dict):
                print_config(k, v, space + 2)
            else:
                if isinstance(v, list):
                    print(f"{' ' * space} - \033[4m{k}\033[0m:")
                    for i in v:
                        print(f"{'   ' * space}{i}")
                else:
                    print(f"{' ' * space} - \033[4m{k}\033[0m: ", end="")
                    print(v)
    else:
        print(f"\033[1m{key}:\033[0m {value}")


print("=" * 80)
print("Primetime: TF reporter pipeline")
print("Version", __version__)
print("-" * 80)
for key, value in config.items():
    print_config(key, value)
print("=" * 80)


################################################################################
### RULES ######################################################################
################################################################################
# 1) Get barcode counts
rule get_barcode_counts:
    input:
        lambda wildcards: get_samples_data(config, wildcards),
    output:
        txt=os.path.join(output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.txt"),
        stats=os.path.join(
            output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.stats"
        ),
    params:
        script=os.path.join(scripts_dir, "primetime_get_barcode_counts.py"),
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
        os.path.join(output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.txt"),
    output:
        os.path.join(
            output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.txt"
        ),
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
        os.path.join(
            output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.txt"
        ),
    output:
        os.path.join(
            output_dir,
            "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.annotated.txt",
        ),
    params:
        script=os.path.join(scripts_dir, "primetime_annotate_barcodes.R"),
        bc_annotation=config["BARCODE_ANNOTATION_FILE"],
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
        os.path.join(output_dir, "tmp_primetime/design.txt"),
    run:
        samples_dict = config["INPUT_DATA"]
        reference_condition = config["COMPARATIVE_ANALYSIS"]["REFERENCE_CONDITION"]
        contrast_condition = config["COMPARATIVE_ANALYSIS"]["CONTRAST_CONDITION"]
        with open(output[0], "w") as f:
            f.write("sample\treplicate\tpDNA\ttreatment\n")
            for sample_name, info in config["INPUT_DATA"].items():
                if info["condition"] == reference_condition:
                    is_pdna = info["is_pdna"]
                    for replicate_number in range(len(info["fastq"])):
                        replicate = f"{sample_name}_{replicate_number+1}"
                        treatment = info["condition"]
                        f.write(f"{sample_name}\t{replicate}\t{is_pdna}\t{treatment}\n")
            for sample_name, info in config["INPUT_DATA"].items():
                if info["condition"] == contrast_condition:
                    is_pdna = info["is_pdna"]
                    for replicate_number in range(len(info["fastq"])):
                        replicate = f"{sample_name}_{replicate_number+1}"
                        treatment = info["condition"]
                        f.write(f"{sample_name}\t{replicate}\t{is_pdna}\t{treatment}\n")
            for sample_name, info in config["INPUT_DATA"].items():
                if (
                    info["condition"] != reference_condition
                    and info["condition"] != contrast_condition
                ):
                    is_pdna = info["is_pdna"]
                    for replicate_number in range(len(info["fastq"])):
                        replicate = f"{sample_name}_{replicate_number+1}"
                        treatment = info["condition"]
                        f.write(f"{sample_name}\t{replicate}\t{is_pdna}\t{treatment}\n")


# 4) Get activity and QC plots
rule get_activity_and_qc_plots:
    input:
        dfs=path_to_annotated_dfs,
        design=os.path.join(output_dir, "tmp_primetime/design.txt"),
    output:
        bc_corr=os.path.join(output_dir, "primetime_QC/barcode_correlations.pdf"),
        replicate_corr=os.path.join(
            output_dir, "primetime_QC/replicate_correlations.pdf"
        ),
        expected_vs_obs=os.path.join(
            output_dir, "primetime_QC/expected_vs_observed_pDNA_counts.pdf"
        ),
        read_counts=os.path.join(output_dir, "primetime_QC/read_counts.pdf"),
        bc_counts=os.path.join(output_dir, "primetime_QC/distribution_of_BC_counts.pdf"),
        bleedthrough=os.path.join(
            output_dir, "primetime_QC/bleedthrough_estimation.pdf"
        ),
        cDNA=os.path.join(output_dir, "tmp_primetime/activity/cDNA_counts.txt"),
        pDNA=os.path.join(output_dir, "tmp_primetime/activity/pDNA_counts.txt"),
    params:
        script=os.path.join(scripts_dir, "primetime_get_activity_and_qc_plots.R"),
        df_basedir=os.path.join(output_dir, "tmp_primetime/bc_counts"),
        expected_pdna_counts=config["EXPECTED_PDNA_COUNTS"],
        plots_basedir=os.path.join(output_dir, "primetime_QC/"),
        activity_basedir=os.path.join(output_dir, "tmp_primetime/activity"),
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        mkdir -p {params.plots_basedir}
        mkdir -p {params.activity_basedir}

        Rscript {params.script} \
        --df_basedir {params.df_basedir} \
        --plots_basedir {params.plots_basedir} \
        --activity_basedir {params.activity_basedir} \
        --design {input.design} \
        --expected_pdna {params.expected_pdna_counts}
        """


# 5) Run comparative analysis
rule run_comparative_analysis:
    input:
        cDNA=os.path.join(output_dir, "tmp_primetime/activity/cDNA_counts.txt"),
        pDNA=os.path.join(output_dir, "tmp_primetime/activity/pDNA_counts.txt"),
        design=os.path.join(output_dir, "tmp_primetime/design.txt"),
    output:
        txt=os.path.join(output_dir, "primetime_results/primetime_results.txt"),
        plots=os.path.join(output_dir, "primetime_results/primetime_volcano.pdf"),
        loli=os.path.join(output_dir, "primetime_results/primetime_lollipop.pdf"),
    params:
        script=os.path.join(scripts_dir, "primetime_run_comparative_analysis_bcalm.R"),
        out_basedir=os.path.join(output_dir, "tmp_primetime/activity"),
        pval_threshold=config["PVALUE_THRESHOLD"],
        plot_output_dir=os.path.join(output_dir, "primetime_results"),
        contrast_condition=config["COMPARATIVE_ANALYSIS"]["CONTRAST_CONDITION"],
        reference_condition=config["COMPARATIVE_ANALYSIS"]["REFERENCE_CONDITION"],
    conda:
        os.path.join(conda_envs_dir, "comparative_analysis.yaml")
    threads: 100
    shell:
        """
        mkdir -p {params.plot_output_dir}
        
        Rscript {params.script} \
        --cdna {input.cDNA} \
        --pdna {input.pDNA} \
        --design {input.design} \
        --output {params.out_basedir} \
        --pval_threshold {params.pval_threshold} \
        --contrast_condition {params.contrast_condition} \
        --reference_condition {params.reference_condition} \
        --plot_output {params.plot_output_dir}
        """


# 6) Get barcode counts for the top-responding TFs
rule get_bc_counts_for_top_tfs:
    input:
        txt=os.path.join(output_dir, "primetime_results/primetime_results.txt"),
        design=os.path.join(output_dir, "tmp_primetime/design.txt"),
        cdna=os.path.join(output_dir, "tmp_primetime/activity/cDNA_counts.txt"),
    output:
        os.path.join(output_dir, "primetime_results/BC_counts_for_top_TFs.pdf"),
    params:
        script=os.path.join(scripts_dir, "primetime_get_bc_counts_for_top_tfs.R"),
        reference_condition=config["COMPARATIVE_ANALYSIS"]["REFERENCE_CONDITION"],
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        Rscript {params.script} \
        --result {input.txt} \
        --design {input.design} \
        --cdna {input.cdna} \
        --output {output} \
        --reference_condition {params.reference_condition}
        """
