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
    )
)
opt <- parse_args(OptionParser(option_list = option_list))

# ==============================================================================
# Read the activity and corrected activity files
# The background activity is calculated as the mean of the negative controls
message("==== Merging corrected activity with MPRAnalyze results")
read.table(opt$corrected_activity, header = T) %>% 
        as_tibble() %>%
        filter(negative_control) %>%
        group_by(replicate, promoter) %>%
        summarize(
                background_activity = mean(activity_RPM),
                .groups = 'drop') -> controls
# Read the TF activity
read.table(opt$corrected_activity, header = T) %>%
    as_tibble() %>%
    filter(!negative_control) -> tfs

# Normalize the activity: corrected_activity_RPM / background_activity
# VF241125: Trying without the correction
tfs %>%
        left_join(controls, by = c('replicate', 'promoter')) %>%
        mutate(
                log2corrected_activity_RPM = log2(activity_RPM + 1) # / background_activity
        ) %>%
        select(replicate, tf, promoter, log2corrected_activity_RPM) -> log2corrected_activity_RPM

df = read.table(opt$mpra_results, header = T) %>% as_tibble() %>%
        left_join(log2corrected_activity_RPM, by = c('tf')) %>% 
        na.omit()
# ==============================================================================

# Save the results
write.table(df,
    file = file.path(opt$output, "primetime_results.txt"),
    sep = "\t", quote = F,
    row.names = F
)
# ==============================================================================
# Merge with the design DF for the MA plot
design = read.table(opt$design, header = T) %>% as_tibble()
# replicate in df, sample in design
df_MA = df %>%
    left_join(design, by = c("replicate" = "sample")) %>%
    filter(treatment == opt$reference_condition) %>%
    mutate(log_activity = log2corrected_activity_RPM)  %>%
    # Remove duplicated rows
    distinct(tf, .keep_all = T) 

# print on screen
print(df_MA, nrows = 50)
# ==============================================================================
df_ns = df %>% filter(sig == "NS")
df_up = df %>% filter(sig == "Upregulated")
df_down = df %>% filter(sig == "Downregulated")

df_MA_ns = df_MA %>% filter(sig == "NS")
df_MA_up = df_MA %>% filter(sig == "Upregulated")
df_MA_down = df_MA %>% filter(sig == "Downregulated")

pdf(
    file = file.path(opt$output, "primetime_results.pdf"),
    width = 8, height = 10
)
ggplot() +
    geom_quasirandom(
        data = df_ns,
        mapping = aes(x = replicate, y = log2corrected_activity_RPM),
        color = "grey", alpha = 0.5
    ) +
    geom_label_repel(
        data = df_ns,
        mapping = aes(x = replicate, y = log2corrected_activity_RPM, label = tf),
        box.padding = 1,
        max.overlaps = 3,
        color = "grey"
    ) +
    geom_point(
        data = df_up,
        mapping = aes(x = replicate, y = log2corrected_activity_RPM),
        colour = "#6495ed", alpha = 1
    ) +
    geom_label_repel(
        data = df_up,
        mapping = aes(x = replicate, y = log2corrected_activity_RPM, label = tf),
        box.padding = 1,
        color = "#6495ed"
    ) +
    geom_point(
        data = df_down,
        mapping = aes(x = replicate, y = log2corrected_activity_RPM),
        colour = "#f37f80", alpha = 1
    ) +
    geom_label_repel(
        data = df_down,
        mapping = aes(x = replicate, y = log2corrected_activity_RPM, label = tf),
        box.padding = 1,
        color = "#f37f80"
    ) +
    ggpubr::theme_pubr(border = T) +
    labs(
        x = "",
        y = "Activity (log2(RPM+1))",
        title = "Activity of MPRAnalyze hits"
    ) +
    theme(text = element_text(size = 14))

#  Now, the MA plot
ggplot(df_MA, aes(y = logFC, x = log_activity, color = sig, alpha = sig)) +
    geom_point() +
    scale_color_manual(values = c("NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80"), name = "") +
    scale_alpha_manual(values = c("NS" = 0.5, "Upregulated" = 1, "Downregulated" = 1), name = "") +
    geom_label_repel(
        data = df_MA_up,
        mapping = aes(label = tf, y = logFC, x = log_activity),
        box.padding = 1,
        color = "#f37f80"
    ) +
    geom_label_repel(
        data = df_MA_down,
        mapping = aes(label = tf, y = logFC, x = log_activity),
        box.padding =0.5,
        max.overlaps = Inf,
        color = "#6495ed"
    ) +
    geom_label_repel(
        data = df_MA_ns,
        mapping = aes(label = tf, y = logFC, x = log_activity),
        box.padding =0.5,
        max.overlaps = 3,
        color = "grey"
    ) +
    # facet_wrap(~treatment, ncol=1) +
    ggpubr::theme_pubr(border = T) +
    labs(
        y = "Log2FC",
        x = paste("Activity in", opt$reference_condition, "condition (log2(RPM+1))"),
        title="MA plot of MPRAnalyze hits"
    ) +
    theme(text = element_text(size = 14))
invisible(dev.off())
# ==============================================================================
# Lollipop plot: diff in activity
library(ggnewscale)
pdf(
    file = file.path(opt$output, "primetime_results_lollipop.pdf"),
    width = 21, height = 7
)
message("==== Creating lollipop plot")
       df %>%
           left_join(design, by = c("replicate" = "sample")) -> tmp_df
        tmp_df %>%
           mutate(log_activity = log2corrected_activity_RPM) %>%
           group_by(tf, treatment) %>%
           summarise(
               mean_act = mean(log2corrected_activity_RPM),
               sig = first(sig),
               .groups = "drop"
           ) %>%
                # Order by mean
                pivot_wider(names_from = treatment, values_from = mean_act) %>%
                mutate(color_axis = "orange") %>%
                arrange(desc(!!sym(opt$reference_condition))) %>%
                mutate(tf = factor(tf, levels = unique(tf))) -> plot_df
                
        
        plot_df %>%
            ggplot() +
            geom_segment(aes(
                y = !!sym(opt$reference_condition), 
                yend = !!sym(opt$contrast_condition),
                x = tf, xend = tf,
                color = sig
            ), size = 1) +
            scale_color_manual(
                values = c("NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80"),
                name = "Result of\nComparative\nAnalysis"
            ) +
            new_scale_color() +
            geom_point(aes(x = tf, y = !!sym(opt$reference_condition), color = "C"), size = 3) +
            geom_point(aes(x = tf, y = !!sym(opt$contrast_condition), color = sig), size = 3) +
            scale_color_manual(values = c('C'= "black", "NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
            guides(color = F) +
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
                title= paste(tmp_df$replicate %>% as.factor() %>% levels(), collapse=' vs. ')
            )
invisible(dev.off())
