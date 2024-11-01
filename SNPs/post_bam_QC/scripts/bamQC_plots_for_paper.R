library(ggplot2)
library(tidyverse)
library(ggtext)
library(ggsci)
library(ggpubr)
library(patchwork)

setwd("/mnt/hdd/dminh/Cauris/SNPs/post_bam_QC")

#  STAR alignment plot + Dedup stats
star = read.csv("plots_data/star_alignment_plot.csv", header = T, stringsAsFactors = F)
colnames(star) = c("strain", "Uniquely mapped", "Mapped to multiple loci", "Mapped to too many loci", "Unmapped: too many mismatches", "Unmapped: too short", "Unmapped: other")
colors = c("#A91D3A", "#D04848", "#FF407D", "#FBA834", "#378CE7", "#10439F")
l_star = star |> pivot_longer(
    cols = c(2:7),
    names_to = "type",
    values_to = "percent"
) |>
  mutate(
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
    )
  )
l_star$type = factor(
    l_star$type,
    levels = c("Unmapped: other", "Unmapped: too short", "Unmapped: too many mismatches", "Mapped to too many loci", "Mapped to multiple loci", "Uniquely mapped")
)
p_align_stat = ggplot(
    data = l_star,
    aes(
        x = strain,
        y = percent,
        fill = type
    )
) +
  geom_bar(
    position = "fill",
    stat = "identity"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x){paste0(x * 100, "%")}
  ) +
  coord_flip() +
  scale_fill_manual(values = colors) +
  labs(
    y = "Percentages", x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank()
  )
ggsave("img/star_align_stats.pdf", units = "in", width = 10, height = 5, dpi = 300, device = "pdf")

## Dedup plot
dedup_stat = read.csv("plots_data/picard_deduplication.csv", header = T, stringsAsFactors = F, check.names = F)
l_dedup_stat = dedup_stat |>
  rename(strain = colnames(dedup_stat)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "type"
  ) |>
  mutate(
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
    )
  )
p_dedup_stat = ggplot(
    data = l_dedup_stat,
    aes(
        x = strain,
        y = counts,
        fill = type
    )
) +
  geom_bar(
    position = "fill",
    stat = "identity"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x){paste0(x * 100, "%")}
  ) +
  coord_flip() +
  scale_fill_npg() +
  labs(
    y = "Percentages"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

align_fig = p_align_stat / p_dedup_stat
align_fig = align_fig +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 18))
ggsave("img/fig/align_figs.png", plot = align_fig, units = "in", width = 10, height = 6, dpi = 300, device = "png")

# Feature + Biotype count
## Biotype counts
bio_count = read.csv("plots_data/featureCounts_biotype_plot.csv", header = T, stringsAsFactors = F, check.names = F)
l_bio_count = bio_count |>
  rename(strain = colnames(bio_count)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "type"
  ) |>
  mutate(
    type = factor(str_replace(type, "_", " ")),
    type = fct_relevel(
      type,
      c("Unassigned Ambiguity", "Unassigned NoFeatures", "tRNA pseudogene", "rRNA", "tRNA", "protein coding")
    ),
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
    )
  )
p_biotype = ggplot(
    data = l_bio_count,
    aes(
        x = strain,
        y = counts,
        fill = type
    )
) +
  geom_bar(
    position = "fill",
    stat = "identity"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("#A91D3A", "#7469B6", "orange", "green", "black", "#7BC9FF")
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = function(x){paste0(x * 100, "%")}) +
  labs(
    y = "Percentages"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
## Feature assignment
feature_count = read.csv("plots_data/featureCounts_assignment_plot.csv", header = T, stringsAsFactors = F, check.names = F)
l_feature_count = feature_count |>
  rename(strain = colnames(feature_count)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "type"
  ) |>
  mutate(
    type = factor(str_replace(type, "_", " ")),
    type = fct_relevel(
      type,
      c("Unassigned Ambiguity", "Unassigned NoFeatures", "Unassigned MultiMapping", "Assigned")),
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
      )
  )
p_feature = ggplot(
    data = l_feature_count,
    aes(
        x = strain,
        y = counts,
        fill = type
    )
) +
  geom_bar(
    position = "fill",
    stat = "identity"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = function(x){paste0(x*100, "%")}) +
  coord_flip() +
  scale_fill_manual(
    values = c("orange", "#90ed7d", "#322C2B", "#7BC9FF")
  ) +
  labs(
    y = "Percentages"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
feature_fig = p_biotype / p_feature
feature_fig = feature_fig +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 18))
ggsave("img/fig/feature_figs.png", plot = feature_fig, units = "in", width = 10, height = 6, dpi = 300, device = "png")

# Complexity curve and DupRadar
## Complexity curve
preseq_stat = read.csv("plots_data/preseq_plot.csv", header = T, stringsAsFactors = F, check.names = F)
l_preseq_stat = preseq_stat |>
  rename(total_mol = colnames(preseq_stat)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "strain"
  ) |>
  mutate(
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
      )
  )
complex_curve = ggplot(
    data = l_preseq_stat,
    aes(
        x = total_mol,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain),
    linewidth = 1
  ) +
  geom_segment(
    aes(
      x = 0, xend = 6.88,
      y = 0, yend = 6.88
    ),
    linetype = 2
  ) +
  annotate(
    geom = "richtext",
    size = 6,
    x = 40, y = 6.80,
    label = "_a perfect library where each read is unique_<br>**6.88 M unique molecules**"
  ) +
  scale_x_continuous(breaks = seq(0, 140, 20),
                     labels = function(x){paste0(x, "M")}) +
  scale_y_continuous(breaks = seq(0, 8, 1),
                     labels = function(x){paste0(x, "M")}) +
  scale_color_npg() +
  labs(
    x = "Total molecules (including duplicates)",
    y = "Unique molecules"
  ) +
  theme_pubclean() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
## DupRadar
dup_stat = read.csv("plots_data/mqc_hcplot_bsqlomaevp.csv", header = T, stringsAsFactors = F)
l_dup_stat = dup_stat |>
  rename(expression = colnames(dup_stat)[1]) |>
  pivot_longer(
    cols = starts_with("strain"),
    values_to = "counts",
    names_to = "strain"
  ) |>
  mutate(
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
      )
  ) |>
  filter(!is.na(counts))
dupradar_p = ggplot(
    data = l_dup_stat,
    aes(
        x = expression,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain),
    linewidth = 1
  ) +
  geom_vline(
    xintercept = c(0.5, 1000),
    linetype = "dashed",
    color = c("#367E18", "#DF2E38")
    ) +
  # log10 scaling and label with exponential format
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_color_npg() +
  labs(
    x = "expression (reads/kbp)",
    y = "% duplicate reads"
  ) +
  theme_pubclean() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  )
complex_dupradar = complex_curve / dupradar_p
complex_dupradar = complex_dupradar +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 18))
ggsave("img/fig/complex_dupradar_figs.png", plot = complex_dupradar, units = "in", width = 14, height = 14, dpi = 300, device = "png")

# BAM coverage
## Cumulative coverage plot
cum_cov = read.delim("/mnt/hdd/dminh/Cauris/SNPs/post_bam_QC/plots_data/mosdepth-cumcoverage-dist-id.txt", sep = "\t", stringsAsFactors = F, check.names = F)
l_cum_cov = cum_cov |>
  pivot_longer(
    cols = !Sample,
    names_to = "coverage",
    values_to = "cov_pct"
  ) |>
  as.data.frame() |>
  mutate(
    pct = as.numeric(str_extract(cov_pct, "\\d+\\.\\d+")),
    coverage = as.numeric(coverage),
    Sample = case_when(
      Sample == "strain_185" ~ "185-R",
      Sample == "strain_188" ~ "188-S",
      Sample == "strain_189" ~ "189-R",
      Sample == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
    )
  )
cum_cov_p = ggplot(
  data = l_cum_cov,
  aes(
    x = coverage,
    y = pct,
    color = Sample
  )
) +
  geom_line(size = 1) +
  geom_vline(
    xintercept = c(100, 400),
    linetype = "dashed",
    alpha = 0.4
    ) + 
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = function(x){paste0(x, "%")}
  ) +
  scale_x_continuous(
    breaks = seq(0, 600, 100),
    labels = function(x){paste0(x, "X")}
  ) +
  theme_minimal() +
  scale_color_npg() +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 24, vjust = 10),
    axis.title.y = element_text(size = 20),
    legend.position = "top",
    legend.justification = c(4,4),
    legend.text = element_text(size = 22)
  ) +
  labs(
    x = "Cumulative Coverage (X)",
    y = "% bases in genome\ncovered by at least X reads"
  )

## Average coverage per contig
avg_cov_contig = read.delim("/mnt/hdd/dminh/Cauris/SNPs/post_bam_QC/plots_data/mosdepth-coverage-per-contig-multi.txt", sep = "\t", header = TRUE, stringsAsFactors = F, check.names = F)
l_avg_cov_contig = avg_cov_contig |>
  pivot_longer(
    cols = !Sample,
    names_to = "no_scaffold",
    values_to = "scaff_cov"
  ) |>
  mutate(
    scaff_cov = str_replace_all(scaff_cov, "[()']", "")
  ) |>
  separate(
    scaff_cov,
    into = c("scaffold", "coverage"),
    sep = ", ",
    convert = TRUE
  ) |>
  mutate(
    scaffold = str_replace(scaffold, "\\.1", ""),
    # Refactor level using natural sorting
    scaffold = fct_relevel(scaffold, unique(str_sort(scaffold, numeric = T))),
    Sample = case_when(
      Sample == "strain_185" ~ "185-R",
      Sample == "strain_188" ~ "188-S",
      Sample == "strain_189" ~ "189-R",
      Sample == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
    )
  ) |>
  as.data.frame()
contig_cov_p = ggplot(
  data = l_avg_cov_contig,
  aes(
    x = scaffold,
    y = coverage,
    group = Sample,
    color = Sample,
  )
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(
    breaks = seq(50, 350, 50),
    labels = function(x){paste0(x, "X")}
  ) +
  theme_minimal() +
  scale_color_npg() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
    legend.position = "none",
    axis.ticks.y = element_blank()
  ) +
  labs(
    y = "Average Coverage"
  )
cov_fig = (cum_cov_p + contig_cov_p) / variant_DP_p
cov_fig = cov_fig +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 24))
cov_fig
ggsave("img/fig/coverage_summary_figs.pdf", plot = cov_fig, units = "in", width = 18, height = 10, dpi = 300, device = "pdf")
