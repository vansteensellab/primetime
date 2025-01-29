# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script runs the MPRAnalyze analysis for the TF reporter data.
# It reads the barcode counts, processes the data, and performs differential
# TF activity analysis comparing different conditions of the samples.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================

suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(MPRAnalyze)
        library(reshape2)
        library(BiocParallel)
        library(optparse)
        library(tidyr)
        library(ggrepel)
})

options(dplyr.width = Inf)



##########################################################################################
# Read arguments #########################################################################
##########################################################################################


option_list <- list(
        make_option(c("-p", "--pdna"), type = "character", default = NULL, help = "Path to the pDNA counts file"),
        make_option(c("-c", "--cdna"), type = "character", default = NULL, help = "Path to the cDNA counts file"),
        make_option(c("-d", "--design"), type = "character", default = NULL, help = "Path to the design file"),
        make_option(c("-o", "--output"), type = "character", default = NULL, help = "Path to the output directory"),
        make_option(c("-t", "--threads"), type = "integer", default = 1, help = "Number of threads to use"),
        make_option(c("-q", "--pval_threshold"), type = "numeric", default = 0.05, help = "P-value threshold for significance"),
        make_option(c("--contrast_condition"), type = "character", default = NULL, help = "Condition to contrast"),
        make_option(c("--reference_condition"), type = "character", default = NULL, help = "Condition to contrast"),
        make_option(c("--corrected_activity"), type = "character", default = NULL, help = "Path to the corrected activity file"),
        # flag to tell wheter to make a single model for all the TFs or not
        make_option(c("--single_model"), type = "logical", default = FALSE, help = "Whether to make a single model for all the TFs or not")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


##########################################################################################
# Preparing the data #####################################################################
##########################################################################################

pdna <- read.table(opt$pdna, header = TRUE, sep = "\t")
cdna <- read.table(opt$cdna, header = TRUE, sep = "\t")
name_of_the_pdna_replicate <- colnames(pdna) %>% setdiff(c("barcode", "negative_control"))
df <-
        merge(pdna, cdna %>% select(-negative_control), by = "barcode") %>%
        group_by(tf) %>%
        # reframe and keep all the other columns
        mutate(
                barcode_id = paste0("bc", row_number())
        ) %>%
        # melt the data to have the columns: tf, promomer, barcode_id, replicate, count
        melt(
                id.vars = c("tf", "barcode_id", "promoter", "barcode", "negative_control"),
                variable.name = "replicate",
                value.name = "count"
        )

# Parse the design DF
design_df <- read.table(opt$design, header = TRUE, sep = "\t", stringsAsFactors = T)
rownames(design_df) <- design_df$replicate
design_df %>%
        filter(!as.logical(pDNA)) %>%
        select(-pDNA) -> design_df

# make the matrices for MPRAanalyze
# =================================

# 1. Pares the cDNA df
parsed_cdna <-
        df %>%
        filter(replicate != name_of_the_pdna_replicate) %>%
        mutate(replicate_id = paste(replicate, barcode_id, sep = "_")) %>%
        select(-barcode_id, -barcode, -replicate) %>%
        # Transform back again to every column being a replicate, but now we have the Sample_RepID_bcID
        pivot_wider(names_from = replicate_id, values_from = count) %>%
        as.data.frame()

# Construct a annotation df from the columns of the parsed_cdna
annotation_df <-
        data.frame(
                obs = colnames(parsed_cdna)
        ) %>%
        as_tibble() %>%
        mutate(
                barcode = gsub(".*_", "", obs),
                replicate = gsub("_bc.*", "", obs),
                # condition for U2OS_Calcitriol_12W_1 is U2OS_Calcitriol_12W. So just remove the last number and underscore
                bio_replicate = gsub("_\\d$", "", gsub(".*_", "", gsub("_bc.*", "", replicate)))
        ) %>%
        merge(design_df, by = "replicate")

rownames(annotation_df) <- annotation_df$obs
annotation_df %>% select(-sample, -replicate, -obs) -> annotation_df

# 1. Make the pDNA df
df %>%
        filter(replicate == name_of_the_pdna_replicate) %>%
        mutate(replicate_id = paste(replicate, barcode_id, sep = "_")) %>%
        select(-barcode_id, -barcode, -replicate) %>%
        pivot_wider(names_from = replicate_id, values_from = count) -> tmp_parsed_pdna


# now get the values from tmp_parsed_pdna but with the names of the columns of
# parsed_cdna, matching the barcode id
# =================================
# the idea is that I will iterate over all the columns of the tmp_parsed_pdna (which are
# like pDNA_bc1, pDNA_bc2, etc) then: see which barcode it is (bc1, bc2, etc) and then get
# the column of parsed_cdna that has the same barcode. The output is a df with the same
# columns as parsed_cdna, but with the values of the pDNA
parsed_pdna <- list(negative_controls = tmp_parsed_pdna$negative_control, promoter = tmp_parsed_pdna$promoter)
columns_of_parsed_cdna <- rownames(annotation_df)
for (col in colnames(tmp_parsed_pdna) %>% setdiff(c("tf", "promoter", "negative_control"))) {
        barcode <- gsub(".*_", "", col)
        # print(barcode)
        # print('===========')
        # See which column of parsed_cdna has the same barcode and return them
        columns_w_same_bc <- columns_of_parsed_cdna[columns_of_parsed_cdna %>% grep(barcode, .)]
        # print(columns_w_same_bc)
        for (col_w_same_bc in columns_w_same_bc) {
                parsed_pdna[[col_w_same_bc]] <- as.numeric(tmp_parsed_pdna[[col]])
        }
}
parsed_pdna <- do.call(cbind, parsed_pdna) %>% as.data.frame()

# Convert treatment and bio_replicate to factors
annotation_df %>%
        # Add AAA_ in front of the row where treatment is reference condition
        mutate(
                treatment =
                        ifelse(treatment == opt$reference_condition,
                                paste0("000_", treatment),
                                ifelse(treatment == opt$contrast_condition,
                                        paste0("ZZZ_", treatment),
                                        treatment
                                )
                        )
        ) %>%
        # Add BBB_ in front of the row where treatment is contrast condition
        mutate(treatment = ifelse(treatment == opt$contrast_condition, paste0("ZZZ_", treatment), treatment)) %>%
        arrange(factor(treatment)) %>%
        mutate(
                treatment = as.factor(treatment),
                bio_replicate = as.factor(bio_replicate)
        ) -> annotation_df


# Checking if both treatments contain one replicate. If true, some changes might be needed
# downstream
one_rep_per_treatment <- (nrow(annotation_df %>%
        select(treatment, bio_replicate) %>%
        distinct()) == 2)
message("One replicate per treatment? ", one_rep_per_treatment)

if (one_rep_per_treatment) {
        annotation_df <- annotation_df %>% select(-bio_replicate)
}

# Remember, parsed_pdna has the promter column that should be removed
negative_controls <- data.frame(
        promoter = parsed_pdna$promoter,
        neg_ctrls = as.logical(parsed_pdna$negative_controls)
)



##########################################################################################
# Running MPRAanalyze ####################################################################
##########################################################################################

message("==== Performing comparative analysis\n")
run_differential_analysis <- function(pdna_matrix,
                                      cdna_matrix, annotation_df, negative_controls,
                                      contrast_condition = opt$contrast_condition,
                                      reference_condition = opt$reference_condition,
                                      output = opt$output, promoter_name,
                                      p_threshold = opt$pval_threshold) {
        message(paste0("==== Running MPRAnalyze for ", promoter_name))
        # Print how many TFs are being analyzed
        message(paste0("==== Number of TFs for ", promoter_name, ": ", nrow(cdna_matrix)))

        # message("pDNA matrix")
        # print(pdna_matrix %>% as_tibble() %>% head(10))

        # message("\n\ncDNA matrix")
        # print(cdna_matrix %>% as_tibble() %>% head(10))

        # message("Annotation df")
        # print(annotation_df)
        obj <- MpraObject(
                dnaCounts = pdna_matrix,
                rnaCounts = cdna_matrix,
                dnaAnnot = annotation_df,
                rnaAnnot = annotation_df,
                controls = negative_controls,
                BPPARAM = MulticoreParam(100)
        )

        # Normalize
        # -- Considers bio_replicate only if there is more than one replicate per treatment
        if (one_rep_per_treatment) {
                obj <- estimateDepthFactors(obj,
                        lib.factor = c("treatment"),
                        which.lib = "dna",
                        depth.estimator = "uq"
                )
        } else {
                obj <- estimateDepthFactors(obj,
                        lib.factor = c("bio_replicate", "treatment"),
                        which.lib = "dna",
                        depth.estimator = "uq"
                )
        }
        obj <- estimateDepthFactors(obj,
                lib.factor = c("treatment"),
                which.lib = "rna",
                depth.estimator = "uq"
        )

        # Now, perform the analysis
        # Again, check if there is only one replicate per treatment.
        # If true, do not consider bio_replicate
        if (one_rep_per_treatment) {
                message("One replicate per treatment: will not consider replicates for model")
                obj <- analyzeComparative(
                        obj = obj,
                        dnaDesign = ~ barcode + treatment,
                        rnaDesign = ~treatment,
                        correctControls = TRUE,
                        fit.se = T,
                        BPPARAM = MulticoreParam(100)
                )
        } else {
                obj <- analyzeComparative(
                        obj = obj,
                        dnaDesign = ~ barcode + bio_replicate + treatment,
                        rnaDesign = ~treatment,
                        correctControls = TRUE,
                        fit.se = T,
                        BPPARAM = MulticoreParam(100)
                )
        }

        rds_name <- paste0(promoter_name, "_MPRAnalyze_obj.rds")
        # Save obj as RDS
        saveRDS(obj, file = file.path(opt$output, rds_name))
        # mpranalyze_result <- testLrt(obj)
        mpranalyze_result <- testCoefficient(obj,
                factor = "treatment",
                contrast = paste0("ZZZ_", contrast_condition)
        ) %>%
                mutate(tf = rownames(.))

        return(mpranalyze_result)
}


# Check if the MPRAnalyse should be run for all promoters or not
if (opt$single_model == TRUE) {
        promoters <- TRUE
} else {
        # Run the analysis for the promoter
        promoters <- parsed_cdna %>%
                select(promoter) %>%
                distinct() %>%
                pull()
}

mpranalyze_result <- data.frame()
for (this_promoter in promoters) {
        if (this_promoter == TRUE) {
                # Use all data without filtering
                this_cdna <- parsed_cdna %>%
                        select(-promoter)
                this_pdna <- parsed_pdna %>%
                        select(-promoter)
                this_tmp_parsed_pdna <- tmp_parsed_pdna
                negative_controls_this_promoter <- negative_controls %>%
                        select(neg_ctrls) %>%
                        pull()
        } else {
                # Select data for this promoter
                this_cdna <- parsed_cdna %>%
                        filter(promoter == this_promoter) %>%
                        select(-promoter)
                this_pdna <- parsed_pdna %>%
                        filter(promoter == this_promoter) %>%
                        select(-promoter)
                this_tmp_parsed_pdna <- tmp_parsed_pdna %>%
                        filter(promoter == this_promoter)
                negative_controls_this_promoter <- negative_controls %>%
                        filter(promoter == this_promoter) %>%
                        select(neg_ctrls) %>%
                        pull()
        }

        # Updating rownames
        rownames(this_pdna) <- this_tmp_parsed_pdna$tf
        rownames(this_cdna) <- this_cdna$tf
        names(negative_controls_this_promoter) <- this_tmp_parsed_pdna$tf

        # Convert to matrix
        this_pdna %>%
                select(-negative_controls) %>%
                mutate_all(as.numeric) %>%
                as.matrix() -> this_pdna_matrix
        this_cdna %>%
                select(-tf, -negative_control) %>%
                mutate_all(as.numeric) %>%
                as.matrix() -> this_cdna_matrix

        # Run the analysis
        this_result <- run_differential_analysis(
                pdna_matrix = this_pdna_matrix,
                cdna_matrix = this_cdna_matrix,
                annotation_df = annotation_df,
                negative_controls = negative_controls_this_promoter,
                contrast_condition = opt$contrast_condition,
                reference_condition = opt$reference_condition,
                output = opt$output,
                promoter_name = ifelse(this_promoter == TRUE, "all_promoters", this_promoter),
                p_threshold = opt$pval_threshold
        )
        this_result$promoter <- ifelse(this_promoter == TRUE, "all_promoters", this_promoter)
        mpranalyze_result <- rbind(mpranalyze_result, this_result)
}

# Correcting the p-values after running the analysis for all the promoters
message("==== Correcting p-values")
mpranalyze_result <- mpranalyze_result %>%
        mutate(p_adjusted = p.adjust(pval, method = "BH"))

# ========================================================================================
# Getting the corrected_activity per treatment
# ========================================================================================
# This would be done by reading the corrected activity, getting the treatment info for them,
# then generating a activity value for every treatment. Then, the resulting DF will be merged
# with the MPRAnalyse results.

read.table(opt$corrected_activity, header = TRUE, sep = "\t") %>%
    as_tibble() %>%
    mutate(sample = replicate) %>%
    select(tf, corrected_activity_RPM, sample) %>%
    # Now, join with the design_df to get the treatment from the samples
    inner_join(design_df %>% select(sample, treatment) %>% distinct(), by = "sample") %>%
    # Columns: tf, corrected_activity_RPM, sample, condition
    group_by(tf, treatment) %>%
    summarize(
        corrected_activity_RPM = mean(corrected_activity_RPM)
    ) %>%
    # Now, pivot the data to have the condition as columns
    pivot_wider(names_from = treatment, values_from = corrected_activity_RPM) -> corrected_activity_per_treatment


print(corrected_activity_per_treatment %>% head())
# Now, merge the above with the results of MPRAnalyze
mpranalyze_result %>%
    select(tf, pval, p_adjusted) %>%
    inner_join(corrected_activity_per_treatment, by = "tf") %>%
    # Calculate the Log2FoldChange, by doing: log2(contrast_condition) - log2(reference_condition)
    mutate(logFC = log2(!!sym(opt$contrast_condition)) - log2(!!sym(opt$reference_condition))) -> mpranalyze_result

##########################################################################################
# Plot results ###########################################################################
##########################################################################################
p_threshold <- opt$pval_threshold
message(paste0("==== Applying signifficance thresholds: fdr <= ", p_threshold, "\n"))
pdf(file.path(opt$output, "MPRAnalyze_volcano.pdf"), width = 10, height = 10)

# mpranalyze_result %>%
#         head() %>%
#         print()

parsed_results <-
        mpranalyze_result %>%
        mutate(
                # Sig is true if logFC is greater than 1 and fdr is less than 0.05
                sig = ifelse(logFC > 0 & p_adjusted <= p_threshold, "Upregulated",
                        ifelse(logFC <= 0 & p_adjusted <= p_threshold, "Downregulated",
                                "NS"
                        )
                )
        )
message("==== Writing MPRAnalyze results")
message("Reference condition:", opt$reference_condition)
message("Contrast condition:", opt$contrast_condition)

plot_title <- paste(opt$contrast_condition, "vs.", opt$reference_condition, "(p.adjusted <=", p_threshold, ")")
message("Plot title:", plot_title)


print(parsed_results %>% head())
parsed_results %>%
        ggplot(aes(x = logFC, y = -log10(p_adjusted), color = sig, alpha = sig)) +
        geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        geom_point() +
        # sig will be tab-blue
        scale_color_manual(values = c("NS" = "grey", "Upregulated" = "#6495ed", "Downregulated" = "#f37f80")) +
        scale_alpha_manual(values = c("NS" = 0.4, "Upregulated" = 1, "Downregulated" = 1)) +
        ggpubr::theme_pubr(border = T) +
        # annotate the significant points
        geom_label_repel(data = . %>% filter(sig != "NS"), size=5, aes(label = tf), box.padding = 0.5) +
        # TeX on labs
        labs(
                x = expression(log[2]("Fold Change")),
                y = expression(-log[10]("Adjusted p-value")),
                title = plot_title,
        ) +
        guides(color = "none", alpha = "none") +
        theme(text = element_text(size = 18))
invisible(dev.off())

parsed_results %>%
        write.table(file.path(opt$output, "MPRAnalyze_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
