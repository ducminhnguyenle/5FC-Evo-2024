library(tidyverse)
library(ggprism)
library(ggtext)
library(patchwork)

# Define function to read in copy number ratio (cnr) data and plot cnr scatter plot
plot_cnr = function(input_cnr, output_cnr_plot, title) {
    # Read in cnr data
    cnr_data = read_tsv(input_cnr) %>%
        mutate(chromosome = str_replace(chromosome, "\\..*", ""))
    
    # Identify unique scaffolds/contigs for each sample for coloring
    scaffolds = cnr_data %>%
        arrange(chromosome) %>% distinct(chromosome) %>% pull(chromosome)
    scaffold_cols = rep(c("gray", "orange"), length.out = length(scaffolds))
    
    # cnr jitter plot for each scaffold/contig
    cnr_plot = cnr_data %>%
        ggplot(
            aes(x = chromosome, y = log2, color = chromosome)
        ) +
            geom_jitter(
                show.legend = FALSE,
                stroke = 0
            ) +
            geom_hline(
                yintercept = 0,
                color = "red",
                linewidth = 1
            ) +
            scale_y_continuous(
                expand = c(0, 0)
            ) +
            scale_color_manual(
                breaks = scaffolds,
                values = scaffold_cols
            ) +
            ylim(c(-1, 1)) +
            labs(
                x = NULL, y = "log<sub>2</sub>(mean coverage)",
                subtitle = title
            ) +
            theme_prism(
                base_size = 16,
                base_fontface = "plain"
            ) +
            theme(
                axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
                axis.title.y = element_markdown(),
                plot.subtitle = element_text(size = 18, vjust = -1),
                plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm")
            )
    
    # Save the plot
    ggsave(filename = output_cnr_plot, cnr_plot, unit = "in", width = 10, height = 6, dpi = 600)
    ggsave(filename = gsub("pdf$", "png", output_cnr_plot), cnr_plot, unit = "in", width = 10, height = 6, dpi = 300)

    return(cnr_plot)
}
# S48 (WT)
s48_scatter_p = plot_cnr(
    input_cnr = "S48/S48.cnr",
    output_cnr_plot = "S48/S48_cnr_scatter.pdf",
    title = "WT"
)
# S187 (R6)
s187_scatter_p = plot_cnr(
    input_cnr = "S187/S187.cnr",
    output_cnr_plot = "S187/S187_cnr_scatter.pdf",
    title = "R6"
)
# S191 (T1)
s191_scatter_p = plot_cnr(
    input_cnr = "S191/S191.cnr",
    output_cnr_plot = "S191/S191_cnr_scatter.pdf",
    title = "T1"
)

# Combine plots
p = s48_scatter_p / s187_scatter_p / s191_scatter_p
ggsave(filename = "Figure_4B.pdf", width = 10, height = 12, dpi = 600)