suppressPackageStartupMessages({
library("optparse")
library(dplyr)
})
options(warn=-1)

option_list <- list(
    make_option(c("--bc_counts"),
        type = "list", default = NULL,
        help = "path to the barcode counts, after clustering",
        metavar = "path.cluster.txt"
    ),
    make_option(c("--output"),
        type = "character", default = NULL,
        help = "Path to the output file", metavar = "path.csv"
    )
)

opt_parser = OptionParser(
    usage = "\tRscript %prog [options]",
    option_list = option_list,
    description = "Merge the barcode counts with the barcode annotation")

opt = parse_args(opt_parser)

list_of_bc_files <- strsplit(opt$bc_counts, " ")[[1]]
merged_df = data.frame()
for (bc_file in list_of_bc_files) {
    read.table(bc_file, sep = "\t", header = TRUE) %>%
        mutate(sample_ID = paste(sample, replicate, sep="_rep")) %>%
        select(-sample, -replicate, -pDNA, -treatment, -cell_line) -> data
    reshape2::dcast(data, barcode + tf + promoter + neg_ctrls ~ sample_ID,
        value.var = "count"
    ) %>%
        as_tibble() -> data
    if (nrow(merged_df) == 0) {
        merged_df = data
    } else {
        merged_df = merge(merged_df, data,
            by = c("barcode", "tf", "promoter", "neg_ctrls"),
            all = TRUE
        )
    }
}

write.table(merged_df, opt$output, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)