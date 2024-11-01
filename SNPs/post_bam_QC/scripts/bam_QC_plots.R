library(ggplot2)
library(tidyverse)
library(ggtext)
library(ggsci)
library(ggpubr)
library(patchwork)

setwd("/mnt/hdd/dminh/Cauris/SNPs/post_bam_QC")
#  STAR alignment plot
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
ggplot(
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

# RSeQC read distribution plot
read_dist = read.csv("./rseqc_read_distribution_plot.csv", header = T, stringsAsFactors = F)
colnames(read_dist)[1] = "strain"
l_read_dist = read_dist |> pivot_longer(
    cols = c(2:10),
    names_to = "type",
    values_to = "reads"
)
l_read_dist$type = factor(
    l_read_dist$type,
    levels = c("Other_intergenic", "TES_down_10kb", "TES_down_5kb", "TES_down_1kb", "TSS_up_10kb", "TSS_up_5kb", "TSS_up_1kb", "Introns", "CDS_Exons")
)
ggplot(
    data = l_read_dist,
    aes(
        x = strain,
        y = reads,
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
  #scale_fill_manual(values = colors) +
  labs(
    title = "RSeQC: Read Distribution", y = "Percentages", x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank()
  )
ggsave("img/rseqc_read_distribution_stats.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# RSeQC inner distance plot
inner_dist = read.csv("./rseqc_inner_distance_plot.csv", header = T, stringsAsFactors = F)
colnames(inner_dist)[1] = "inner_distance"
l_inner_dist = inner_dist |> pivot_longer(
    cols = starts_with("strain"),
    values_to = "Counts",
    names_to = "strain"
)
ggplot(
    data = l_inner_dist,
    aes(
        x = inner_distance,
        y = Counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain)
  ) +
  scale_x_continuous(breaks = seq(-250, 250, 50)) +
  scale_y_continuous(breaks = seq(0, 60000, 10000)) +
  scale_color_npg() +
  labs(
    x = "Inner Distance (bp)",
    y = "Counts",
    title = "RSeQC: Inner distance"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
    )
  )
ggsave("img/rseqc_inner_distance_stats.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# RSeQC Junction annotation - Splicing events
spli_event = read.csv("./rseqc_junction_annotation_events_plot.csv", header = T, stringsAsFactors = F)
colnames(spli_event) = c("strain", "Known Splicing Events", "Partial Novel Splicing Events", "Novel Splicing Events")
l_spli_event = spli_event |> pivot_longer(
    cols = contains("Splicing"),
    names_to = "events",
    values_to = "counts"
)
l_spli_event$events = factor(
  l_spli_event$events,
  levels = c("Novel Splicing Events", "Partial Novel Splicing Events", "Known Splicing Events")
)
ggplot(
    data = l_spli_event,
    aes(
        x = strain,
        y = counts,
        fill = events
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
  scale_fill_manual(values = c("#9BCF53", "black", "#3AA6B9")) +
  labs(
    title = "Splicing events", y = "Percentages", x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank()
  )
ggsave("img/rseqc_splicing_events_stats.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# RSeQC Junction annotation - Splicing junctions
spli_junc = read.csv("./rseqc_junction_annotation_junctions_plot.csv", header = T, stringsAsFactors = F)
## Convert to long format
l_spli_junc = spli_junc |>
  rename(strain = colnames(spli_junc)[1]) |>
  pivot_longer(
    cols = contains("splicing"),
    names_to = "events",
    values_to = "counts"
  ) |>
  mutate(
    events = str_replace_all(events,
                             "_", " "),
    events = factor(str_to_title(events)),
    events = fct_relevel(events, c("Novel Splicing Junctions", "Partial Novel Splicing Junctions", "Known Splicing Junctions"))
  )
ggplot(
    data = l_spli_junc,
    aes(
        x = strain,
        y = counts,
        fill = events
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
  scale_fill_manual(
    values = c("#9BCF53", "black", "#3AA6B9")
  ) +
  labs(
    title = "Splicing junctions", y = "Percentages", x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank()
  )
ggsave("img/rseqc_splicing_junctions_stats.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# RSeQC Junction saturation
junc_sat = read.csv("./rseqc_junction_saturation_plot.csv", header = T, stringsAsFactors = F)
colnames(junc_sat)[1] = "percent_of_reads"
l_junc_sat = junc_sat |> pivot_longer(
  cols = starts_with("strain"),
  values_to = "counts",
  names_to = "strain"
)
ggplot(
    data = l_junc_sat,
    aes(
        x = percent_of_reads,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain)
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 45000, 5000)) +
  scale_color_nejm() +
  labs(
    x = "Percent of reads",
    y = "All junctions",
    title = "Junction Saturation"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
    )
  )
ggsave("img/rseqc_junction_saturation_stats.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# RSeQC: Infer experiments
infer_stat = read.csv("./rseqc_infer_experiment_plot.csv", header = T, stringsAsFactors = F)
l_infer_stat = infer_stat |>
  rename(strain = colnames(infer_stat)[1]) |>
  pivot_longer(
    cols = c(2:4),
    values_to = "percentage",
    names_to = "type"
  ) |>
  mutate(type = fct_relevel(type, c("Sense", "Antisense", "Undetermined")))
ggplot(
    data = l_infer_stat,
    aes(
        x = strain,
        y = percentage,
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
  scale_fill_manual(
    values = c("#3AA6B9", "black", "#9BCF53")
  ) +
  labs(
    title = "Infer experiment", y = "Percentages", x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank()
  )
ggsave("img/rseqc_infer_exp_stats.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# Qualimap: Gene coverage profile
qual_stat = read.csv("./qualimap_gene_coverage_profile.csv", header = T, stringsAsFactors = F)
l_qual_stat = qual_stat |>
  rename(trans_pos = colnames(qual_stat)[1]) |>
  pivot_longer(
    cols = starts_with("strain"),
    values_to = "counts",
    names_to = "strain"
  )
ggplot(
    data = l_qual_stat,
    aes(
        x = trans_pos,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain)
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 22500, 2500)) +
  scale_color_npg() +
  labs(
    x = "Transcript position (bp)",
    y = "Counts",
    title = "Coverage Profile Along Genes (total)"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
    )
  )
ggsave("img/qualimap_gene_coverage_profile.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# Deduplication stats
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
ggplot(
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("img/dedup_stats.pdf", units = "in", width = 10, height = 5, dpi = 300, device = "pdf")

# Complexity curve: complexity of a library
preseq_stat = read.csv("./preseq_plot.csv", header = T, stringsAsFactors = F, check.names = F)
l_preseq_stat = preseq_stat |>
  rename(total_mol = colnames(preseq_stat)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "strain"
  )
ggplot(
    data = l_preseq_stat,
    aes(
        x = total_mol,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain)
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
    x = 36, y = 6.88,
    label = "_a perfect library where each read is unique_<br>**6.88 M unique molecules**"
  ) +
  scale_x_continuous(breaks = seq(0, 140, 20),
                     labels = function(x){paste0(x, "M")}) +
  scale_y_continuous(breaks = seq(0, 8, 1),
                     labels = function(x){paste0(x, "M")}) +
  # ylim(c(0, 8)) +
  scale_color_npg() +
  labs(
    x = "Total molecules (including duplicates)",
    y = "Unique molecules",
    title = "Complexity curve"
  ) +
  theme_pubclean() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
    )
  )
ggsave("img/complexity_curve.tiff", units = "in", width = 10, height = 6.5, dpi = 300, device = "tiff", compression = "lzw")

# DupRadar
dup_stat = read.csv("./mqc_hcplot_bsqlomaevp.csv", header = T, stringsAsFactors = F)
l_dup_stat = dup_stat |>
  rename(expression = colnames(dup_stat)[1]) |>
  pivot_longer(
    cols = starts_with("strain"),
    values_to = "counts",
    names_to = "strain"
  ) |>
  filter(!is.na(counts))
  # filter_at(vars(counts), all_vars(!is.na(.)))
  # drop_na(counts)
ggplot(
    data = l_dup_stat,
    aes(
        x = expression,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain)
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
    y = "% duplicate reads",
    title = "A summary of the gene duplication distributions"
  ) +
  theme_pubclean() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
    )
  )
ggsave("img/dupradar_dup_dist.tiff", units = "in", width = 12, height = 6.5, dpi = 300, device = "tiff", compression = "lzw")

# Biotype counts
bio_count = read.csv("./featureCounts_biotype_plot.csv", header = T, stringsAsFactors = F, check.names = F)
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
    )
  )
ggplot(
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
    title = "Genomic features of different biotypes", y = "Percentages", x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("img/biotype_count.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# Feature counts assignment
feature_count = read.csv("./featureCounts_assignment_plot.csv", header = T, stringsAsFactors = F, check.names = F)
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
      c("Unassigned Ambiguity", "Unassigned NoFeatures", "Unassigned MultiMapping", "Assigned"))
  )
ggplot(
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
    title = "Mapped reads summarization for genomic features", y = "Percentages", x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("img/feature_count.tiff", units = "in", width = 10, height = 5, dpi = 300, device = "tiff", compression = "lzw")

# bamQC: Cumulative Coverage
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
head(l_cum_cov)
ggplot(
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
  # annotate(
  #   geom = "richtext",
  #   x = 36, y = 10,
  #   label = "Median coverage: **100X**<br>_Fraction of genome with at least 30X coverage_: **95.1%**"
  # ) +  
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.title = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = "Cumulative Coverage (X)",
    y = "% bases in genome covered by at least X reads",
    title = "Cumulative coverage distribution"
  )
ggsave("img/cum_coverage_dist.tiff", units = "in", width = 14, height = 8, dpi = 300, device = "tiff", compression = "lzw")

# Coverage distribution
cov_dist = read.delim("/mnt/hdd/dminh/Cauris/SNPs/post_bam_QC/plots_data/mosdepth-coverage-dist-id.txt", sep = "\t", header = TRUE, stringsAsFactors = F, check.names = F)
l_cov_dist = cov_dist |>
  pivot_longer(
    cols = !Sample,
    names_to = "coverage",
    values_to = "cov_pct"
  ) |>
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
head(l_cov_dist)
ggplot(
  data = l_cov_dist,
  aes(
    x = coverage,
    y = pct,
    color = Sample
  )
) +
  geom_line(size = 1) +
  # annotate(
  #   geom = "richtext",
  #   x = 36, y = 10,
  #   label = "Median coverage: **100X**<br>_Fraction of genome with at least 30X coverage_: **95.1%**"
  # ) +  
  scale_y_continuous(
    breaks = seq(0, 45, 5),
    labels = function(x){paste0(x, "%")}
  ) +
  scale_x_continuous(
    breaks = seq(0, 500, 100),
    labels = function(x){paste0(x, "X")}
  ) +
  theme_minimal() +
  scale_color_npg() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.title = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = "Coverage (X)",
    y = "% bases in genome covered by X reads",
    title = "Coverage distribution"
  )
ggsave("img/coverage_dist.tiff", units = "in", width = 14, height = 8, dpi = 300, device = "tiff", compression = "lzw")

# Average coverage per contig
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
head(l_avg_cov_contig)
ggplot(
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.title = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    y = "Average Coverage",
    title = "Coverage per contig"
  )
ggsave("img/avg_cov_per_contig.tiff", units = "in", width = 16, height = 8, dpi = 300, device = "tiff", compression = "lzw")
