# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script overlaps the barcode counts with the barcode annotation.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================

suppressPackageStartupMessages({
    library("optparse")
    library(dplyr)
})
options(warn = -1)

option_list <- list(
    make_option(c("--bc_counts"),
        type = "character", default = NULL,
        help = "path to the barcode counts, after clustering",
        metavar = "path.cluster.txt"
    ),
    make_option(c("--bc_annotation"),
        type = "character", default = NULL,
        help = "Path to the barcode annotation file",
        metavar = "path.csv"
    ),
    make_option(c("--sample"),
        type = "character", default = NULL,
        help = "Sample name, for metadata",
        metavar = "sample_name"
    ),
    make_option(c("--replicate_number"),
        type = "numeric", default = NULL,
        help = "Replicate number, for metadata", metavar = "replicate_no"
    ),
    make_option(c("--treatment"),
        type = "numeric", default = NULL,
        help = "Treatment name, for metadata", metavar = "Treatment"
    ),
    make_option(c("--is_pdna"),
        type = "logical", default = NULL,
        help = "Is this sample pDNA?", metavar = "T/F"
    ),
    make_option(c("--output"),
        type = "character", default = NULL,
        help = "Path to the output file", metavar = "path.csv"
    )
)

opt_parser <- OptionParser(
    usage = "\tRscript %prog [options]",
    option_list = option_list,
    description = "Annotate the barcode counts with the barcode annotation"
)
opt <- parse_args(opt_parser)

message("==== Annotating barcodes")
read.table(opt$bc_counts, col.names = c("barcode", "count", "unclustered_bc")) %>%
    as_tibble() -> bc_counts
read.table(opt$bc_annotation, sep = ",", header = TRUE) %>%
    as_tibble() -> bc_annotation
inner_join(bc_counts, bc_annotation, by = join_by("barcode")) %>%
    select(-unclustered_bc) %>%
    mutate(
        sample = opt$sample,
        replicate = opt$replicate_number,
        pDNA = as.logical(opt$is_pdna),
        treatment = opt$treatment,
    ) %>%
    write.table(opt$output, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
