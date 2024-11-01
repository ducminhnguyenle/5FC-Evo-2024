library(vcfR)
library(tidyverse)
library(ggsci)
library(ggtext)
library(gghalves)
library(ggpubr)
library(RColorBrewer)

setwd("/mnt/hdd/dminh/Cauris/SNPs/annotation/marked_filtered/")
filtered_vcf_path = "filtered/cauris.snps.PASS.DP30.PEKT.ann.vcf"
ref = ape::read.dna("/mnt/hdd/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta/genome.fa", format = "fasta")
gff = read.table("/mnt/hdd/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Annotation/Genes/genes.gff", sep = "\t", quote = "")

raw_vcf = read.vcfR("/mnt/hdd/dminh/Cauris/SNPs/annotation/marked_filtered/cauris.snp.filtered.PEKT.ann.vcf", verbose = FALSE)

# Create a full chromR object
full_raw_chrom = create.chromR(name = "", vcf = raw_vcf, seq = ref, ann = gff, verbose = TRUE)
full_raw_chrom = masker(full_raw_chrom, min_QUAL = 1, min_DP = 30)
full_raw_chrom = proc.chromR(full_raw_chrom, verbose = TRUE)

# Overall vcf QC
tiff("SNP_overallqc_vcf_minDP30.tiff", units="in", width=8, height=6, res=300, compression = "lzw")
plot(full_raw_chrom)
dev.off()

# Composite vcf plot
tiff("SNP_vcf_composite_plot.tiff", units="in", width=14, height=10, res=300, compression = "lzw")
chromoqc(full_raw_chrom)
dev.off()

# Filtering
chrom2 = create.chromR(name = "", vcf = raw_vcf, seq = ref, ann = gff, verbose = TRUE)
chrom2 = masker(chrom2, min_QUAL = 1, min_DP = 400, max_DP = 4000, min_MQ = 59.5, max_MQ = 60.5)
chrom2 = proc.chromR(chrom2, verbose = TRUE)
tiff("SNP_overallqc_vcf_DP40-4000.tiff", units="in", width=8, height=6, res=300, compression = "lzw")
plot(chrom2)
dev.off()
tiff("SNP_vcf_composite_plot_DP40-4000.tiff", units="in", width=14, height=10, res=300, compression = "lzw")
chromoqc(chrom2)
dev.off()

# Summarize chromR object
head(full_raw_chrom)

# Extract DP stats in FORMAT column (extract.gt)
dp <- extract.gt(full_raw_chrom, element = "DP", as.numeric = TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)
is.na(dp[na.omit(dp == 0)]) <- TRUE

# Heatmap for DP
tiff("DP_heatmap.tiff", units="in", width=10, height=12, res=300, compression = "lzw")
heatmap.bp(dp, rlabels = F)
dev.off()

# The sequence variant depth for each variant per sample
variant_dp_per_sample_df = dp |>
    as.data.frame() |>
    pivot_longer(
        cols = starts_with("strain"),
        values_to = "DP",
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
# Boxplot
ggplot(
    data = variant_dp_per_sample_df,
    aes(
        x = strain,
        y = DP,
        fill = strain,
    )
) +
    geom_boxplot() +
    theme_minimal() +
    labs(
        y = "Variant depth", title = "Low variant depth between wildtype strain and others"
    ) +
    scale_fill_npg() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5, face ="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank()
    )
ggsave("variant_depth_per_sample_box.png", unit = "in", width = 12, height = 8, device = "png", dpi = 300)
# Violin plot
ggplot(
    data = variant_dp_per_sample_df,
    aes(
        x = strain,
        y = DP,
        fill = strain
    )
) +
    geom_violin(
        adjust = 1.0, scale = "count", trim = F,
        width = 0.8
    ) +
    scale_y_continuous(
        transform = scales::log2_trans(),
        breaks = c(1, 10, 100, 1000, 8000),
        minor_breaks = c(1:10, 2:10*10, 2:10*100, 2:8*1000)
    ) +
    stat_summary(
        fun.y = median,
        geom = "point",
        shape = 23,
        size = 8
    ) +
    theme_minimal() +
    labs(
        y = "Variant depth", title = "Low variant depth in wildtype strain"
    ) +
    scale_fill_npg() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5, face ="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank()
    )
ggsave("variant_depth_per_sample.png", unit = "in", width = 12, height = 8, device = "png", dpi = 300)

# Hybrid plot style
variant_DP_p = ggplot(
    data = variant_dp_per_sample_df,
    aes(
        x = strain,
        y = DP,
        fill = strain
    )
) +
    geom_half_violin(side = "l") +
    geom_half_boxplot(
        side = "r", width = 0.2,
        position = position_nudge(x = 0.05)
    ) +
    geom_half_point_panel(
        aes(color = strain),
        position = position_nudge(x = 0.05)
    ) +
    scale_y_continuous(
        transform = scales::log2_trans(),
        breaks = c(1, 10, 100, 1000, 8000),
        minor_breaks = c(1:10, 2:10*10, 2:10*100, 2:8*1000)
    ) +
    stat_compare_means(
        aes(group = strain),
        label = "p.signif",
        comparisons = list(
            c("185-R", "WT-S"),
            c("188-S", "WT-S"),
            c("189-R", "WT-S"),
            c("191-I", "WT-S")
        )
    ) +
    labs(
        x = "strain", y = "Variant depth"
    ) +
    theme_minimal() +
    scale_fill_npg() +
    scale_color_npg() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5, face ="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
ggsave("img/variant_depth_across_samples.png", unit = "in", width = 20, height = 8, device = "png", dpi = 300)

# Missingness across all samples
# missingness_per_sample = apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
# missingness_per_sample = missingness_per_sample/nrow(raw_vcf)
missingness_per_sample = dp |>
    as_tibble() |>
    summarise(
        across(everything(), ~ sum(is.na(.x)))
    ) |>
    pivot_longer(
        cols = starts_with("strain"),
        values_to = "missingness",
        names_to = "strain"
    ) |>
    mutate(
        miss_pct = missingness/nrow(raw_vcf) *100,
        strain = case_when(
            strain == "strain_185" ~ "185-R",
            strain == "strain_188" ~ "188-S",
            strain == "strain_189" ~ "189-R",
            strain == "strain_191" ~ "191-I",
            TRUE ~ "WT-S"
        )
    )
missingness_per_sample_p = ggplot(
    data = missingness_per_sample,
    aes(
        x = strain, y = miss_pct,
        fill = strain,
        label = miss_pct
    )
) +
    geom_col() +
    geom_text(
        aes(label = round(miss_pct, 2)),
        nudge_y = 0.05,
        size = 8
    ) +
    labs(
        y = "Missingness percent (%)"
    ) +
    theme_classic() +
    scale_fill_npg() +
    theme(
        plot.title = element_text(size = 20, hjust = 0.5, face ="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
# ggsave("img/missingness_per_sample.png", unit = "in", width = 10, height = 6, device = "png", dpi = 300)

# Missingness across variants
# missingness_per_variant <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) })
# missingness_per_variant <- myMiss/ncol(raw_vcf@gt[,-1])
missingness_per_variant = dp |>
    as_tibble() |>
    rowwise() |>
    mutate(
        missingness = sum(is.na(c_across(everything())))
    ) |>
    ungroup() |>
    mutate(
        miss_pct = missingness/ncol(raw_vcf@gt[, -1]) * 100
    ) |>
    select(!contains("strain"))
missingness_per_variant_dp = ggplot(
    data = missingness_per_variant,
    aes(miss_pct)
) +
    geom_histogram(
        aes(y = ..density..),
        fill = "steelblue",
        color = "black", alpha = 0.6
    ) +
    geom_density(
        alpha = 0.6,
        color = "#a063b7",
        fill = "#a063b7"
    ) +
    labs(
        x = "Missingness (%)"
    ) +
    theme_classic() +
    scale_fill_npg() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(hjust = 1, size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, vjust = 10),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
# ggsave("img/missingness_across_variants.png", unit = "in", width = 10, height = 6, device = "png", dpi = 300)

# Combined plots
varDP_missingness_fig = variant_DP_p / (missingness_per_sample_p | missingness_per_variant_dp)
varDP_missingness_fig = varDP_missingness_fig +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 24))
varDP_missingness_fig
ggsave("img/varDP_missingness_figs.png", plot = varDP_missingness_fig, units = "in", width = 16, height = 10, dpi = 300, device = "png")
