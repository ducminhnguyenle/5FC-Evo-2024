library(ggplot2)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(gghalves)
library(ggpol)
library(ggtext)
setwd("/mnt/hdd/dminh/Cauris/AlphaFold/results/")

# Boxplot for 9 poses binding affinity comparison between WT and R1 (R214T)
affi = read.table("PRPP_5FU_simultaneous_docked_result.csv", sep = ",", header = T, stringsAsFactors = F, check.names = F)

l_affi = affi |> pivot_longer(
    cols = c("WT", "R1"),
    names_to = "strain",
    values_to = "affinity"
) |>
    mutate(
        strain = fct_relevel(
            strain,
            c("WT", "R1")
        )
    )
ggplot(
    data = l_affi,
    aes(
        y = abs(affinity),
        x = strain,
        fill = strain
    )
) +
    geom_boxjitter(
        color = "black",
        jitter.shape = 21, jitter.alpha = 0.8,
        jitter.size = 4, jitter.params = list(width = 0.04, height = 0),
        errorbar.draw = T,
        outlier.intersect = T, outlier.shape = 24, outlier.size = 4
    ) +
    scale_y_continuous(
        labels = function(x){paste0("-", x)}
    ) +
    labs(
        y = "binding affinity (kcal/mol)\n5FU + PRPP"
    ) +
    scale_fill_npg() +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
    ) +
    stat_compare_means(
        method = "wilcox.test",
        label.y = 9.8,
        label.x = 1.3
    )
ggsave("binding_affinity.pdf", unit = "in", width = 4, height = 4, device = "pdf", dpi = 300)

# Barplot for best binding affinity comparison between WT and R1 (R214T)
top_affi = l_affi |>
    filter(mode == 1) |>
    mutate(
        strain = fct_relevel(
            strain,
            c("WT", "R1")
        )
    )
ggplot(
    data = top_affi,
    aes(
        x = strain,
        y = abs(affinity),
        fill = strain
    )
) +
    geom_col(width = 0.5) +
    geom_text(
        aes(label = affinity),
        size = 4.5,
        nudge_y = 0.4
    ) +
    scale_fill_npg() +
    scale_y_continuous(
        labels = function(x){paste0("-", x)}
    ) +
    labs(
        y = "Best binding affinity (kcal/mol)\n5FU + PRPP"
    ) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
    )
ggsave("best_binding_affinity.pdf", unit = "in", width = 4, height = 4, device = "pdf", dpi = 300)

# Differences between single and multiple ligands docking
comp = read.table("comparisons.csv", sep = ",", header = T, stringsAsFactors = F, check.names = F)
comp$affinity = abs(comp$affinity)
comp$receptor = factor(comp$receptor, levels = c("WT", "R1"))
p = ggboxplot(
    data = comp, x = "receptor", y = "affinity",
    color = "docking", palette = "jco",
    add = "dotplot",
    add.params = list(dotsize = 0.5), width = 1,
    ylab = "Binding affinity (kcal/mol)",
    bxp.errorbar = T,
    ggtheme = theme_classic()
) +
    scale_fill_npg() +
    scale_color_npg() +
    scale_y_continuous(
        labels = function(x){paste0("-", x)}
    ) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
    )

stat_test = comp |>
    group_by(receptor) |>
    t_test(affinity ~ docking) |>
    add_xy_position(x = "receptor", dodge = 0.8)

p1 = p + stat_pvalue_manual(
    stat_test,
    label = "p.adj",
    tip.length = 0.01
)

ggsave("single_vs_multiple_ligands_docking_comparison.pdf", unit = "in", width = 6, height = 5, device = "pdf", dpi = 300)

# WT vs R1 (R214T) receptor docking
ggplot(
    data = comp,
    aes(
        y = affinity,
        x = docking,
        fill = receptor
    )
) +
    geom_boxjitter(
        color = "black",
        jitter.shape = 21, jitter.alpha = 0.8,
        jitter.size = 2, jitter.params = list(width = 0.04, height = 0),
        errorbar.draw = T,
        outlier.intersect = T, outlier.shape = 24, outlier.size = 2
    ) +
    scale_y_continuous(
        labels = function(x){paste0("-", x)}
    ) +
    labs(
        y = "Binding affinity (kcal/mol)"
    ) +
    scale_fill_npg() +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
    ) +
    stat_compare_means(
        method = "kruskal.test",
        aes(group = docking),
        label.y = 10.5,
        # label.x = 1.3
    ) +
    stat_compare_means(
        method = "wilcox.test",
        label = "p",
        label.y = 9.8,
        aes(group = receptor)
    )
ggsave("single_vs_multiple_ligands_docking_comparison2.pdf", unit = "in", width = 6, height = 5, device = "pdf", dpi = 300)
