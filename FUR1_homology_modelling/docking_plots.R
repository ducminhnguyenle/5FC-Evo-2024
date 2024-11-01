library(ggplot2)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(gghalves)
library(ggpol)
library(ggtext)
library(rstatix)

setwd("autodock_vina/results/")

affi = read.table("binding_affinity.csv", sep = ",", header = T, stringsAsFactors = F, check.names = F)
affi$affinity = abs(affi$affinity)
affi$receptor = factor(affi$receptor, levels = c("WT", "R1", "R4"))

# Barplot for best binding affinity comparison between WT, R1 and R4
affi |>
    group_by(receptor, ligand) |>
    filter(affinity == max(affinity)) |>
    ggplot(
        aes(x = ligand, y = affinity, fill = receptor)
    ) +
        geom_bar(
            position = "dodge",
            stat = "identity"
        ) +
        geom_text(
            aes(label = round(affinity,2)),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 4,
        ) +
        scale_fill_jama() +
        scale_y_continuous(
            labels = function(x){paste0("-", x)}
        ) +
        labs(
            y = "Best binding affinity (kcal/mol)"
        ) +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12, color = "black"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 12)
        )
ggsave("best_binding_affinity.pdf", unit = "in", width = 6, height = 5, device = "pdf", dpi = 300)

# Differences between single and multiple ligands docking
## Calculate pairwise pvalue between ligands
stat_test = affi |>
    group_by(receptor) |>
    wilcox_test(affinity ~ ligand) |>
    add_xy_position(x = "receptor", dodge = 0.8)
## Grouped Boxplot
p = ggboxplot(
    data = affi, x = "receptor", y = "affinity",
    color = "ligand", palette = "jco",
    add = "jitter",
    add.params = list(dotsize = 0.5, alpha = 0.8), width = 1,
    ylab = "Binding affinity (kcal/mol)",
    bxp.errorbar = T,
    ggtheme = theme_classic()
) +
    scale_fill_jama() +
    scale_color_jama() +
    scale_y_continuous(
        labels = function(x){paste0("-", x)}
    ) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12)
    ) +
    stat_pvalue_manual(
        stat_test,
        label = "p.adj.signif",
        tip.length = 0.01
    )
ggsave("single_vs_multiple_ligands_docking_comparison.pdf", unit = "in", width = 6, height = 5, device = "pdf", dpi = 300)

# WT vs R1 (R214T) vs R4 receptor between different docking molecules
## Calculate pairwise pvalue between receptors
stat_test2 = affi |>
    group_by(ligand) |>
    wilcox_test(affinity ~ receptor) |>
    add_xy_position(x = "ligand", dodge = 0.8)
## Grouped Boxplot with pvalue
fig = ggboxplot(
    data = affi, x = "ligand", y = "affinity",
    color = "receptor", add = "jitter",
    add.params = list(dotsize = 0.5, alpha = 0.8),
    ylab = "Binding affinity (kcal/mol)",
    bxp.errorbar = T,
    ggtheme = theme_classic()
) +
    scale_fill_jama() +
    scale_color_jama() +
    scale_y_continuous(
        labels = function(x){paste0("-", x)}
    ) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12)
    ) +
    stat_compare_means(
        method = "kruskal.test",
        aes(group = ligand),
        label.y = 12,
        # label.x = 1.3
    ) +
    stat_pvalue_manual(
        stat_test2, label = "p.adj.signif",
        tip.length = 0.01
    )
ggsave("WT_vs_R1_vs_R4_docking_comparisons.pdf", unit = "in", width = 6, height = 5, device = "pdf", dpi = 300)

