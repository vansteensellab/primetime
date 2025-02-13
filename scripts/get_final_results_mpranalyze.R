# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script overlays the result of MPRAnalyze with the activity of the TFs.
# The result will be stored as txt and pdf.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(ggrepel)
    library(ggbeeswarm)
    library(optparse)
    library(ggnewscale)
})

option_list <- list(
    make_option(c("-c", "--corrected_activity"),
        type = "character", default = NULL,
        help = "Path to the corrected activity file"
    ),
    make_option(c("-m", "--mpra_results"),
        type = "character", default = NULL,
        help = "Path to the MPRA results file"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Path to the output directory"
    ),
    make_option(c('-d', '--design'),
        type = 'character', default = NULL,
        help = 'Path to the design file'
    ),
    # contrast condition and reference conditions
    make_option(c('--reference_condition'),
        type = 'character', default = 'control',
        help = 'Reference condition for the MA plot'
    ),
    make_option(c('--contrast_condition'),
        type = 'character', default = 'treated',
        help = 'Contrast condition for the MA plot'
    ),
    make_option(c('--plot_output'),
        type = 'character', default = NULL,
        help = 'Path to the output directory'
    )
)
opt <- parse_args(OptionParser(option_list = option_list))
df = read.table(opt$mpra_results, header = T)

reference_condition <- opt$reference_condition
contrast_condition <- opt$contrast_condition

# The corrected activity will serve just to filter out the negative controls
aux_df = read.table(opt$corrected_activity, header = T) %>%
    select(tf, negative_control) %>%
    as_tibble()
pdf(
    file = file.path(opt$plot_output, "primetime_lollipop.pdf"),
    width = 21, height = 7
)
message("==== Creating lollipop plot")
       df %>%
           inner_join(aux_df, by = "tf") %>%
           filter(!negative_control) %>%
           mutate(color_axis = ifelse(sig == "Upregulated", "#f37f80",
               ifelse(sig == "NS", "gray30",
                   "#6495ed"
               )
           )) %>%
           # Remove duplicates
           distinct(tf, .keep_all = T) %>%
           arrange(desc(!!sym(reference_condition))) %>%
           mutate(tf = factor(tf, levels = unique(tf))) -> plot_df

print(plot_df)
        
        plot_df %>%
            ggplot() +
            geom_segment(aes(
                y = !!sym(reference_condition), 
                yend = !!sym(contrast_condition),
                x = tf, xend = tf,
                color = sig
            ), size = 1) +
            scale_color_manual(
                values = c("NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80"),
                name = "Result of\nComparative\nAnalysis"
            ) +
            new_scale_color() +
            geom_point(aes(x = tf, y = !!sym(reference_condition), color = "C"), size = 3) +
            geom_point(aes(x = tf, y = !!sym(contrast_condition), color = sig), size = 3) +
            scale_color_manual(values = c('C'= "black", "NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
            guides(color = "none") +
            theme_bw() +
            theme(
                axis.text.x = element_text(
                    angle = 90, hjust = 1, vjust = 0.5,
                    color = plot_df$color_axis
                ),
                # Remove inner lines of y-axis
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                text = element_text(size = 14)
            ) +
            labs(
                y = "Activity (log2(RPM+1))",
                x = "",
                title= paste(contrast_condition, "vs.", reference_condition)
            )
invisible(dev.off())
