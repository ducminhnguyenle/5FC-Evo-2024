library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(ggforce)
library(ggrepel)
library(ggplotify)
library(pheatmap)
library(patchwork)

# Input multiple samples vcf file
## Using "PASS" hard filtering and overall DP >= 30
snp_filtered_vcf = "cauris.snps.PASS.DP30.PEKT.ann.vcf"

## Reformat from vcf to gds object
snpgdsVCF2GDS(snp_filtered_vcf, "cauris.snps.PASS.DP30.PEKT.ann.gds", method = "biallelic.only")

## Summary
snpgdsSummary("cauris.snps.PASS.DP30.PEKT.ann.gds")

# Open a GDS file
snp_filtered_gds = snpgdsOpen("cauris.snps.PASS.DP30.PEKT.ann.gds")

# ---LD-based SNP pruning---

## Different LD thresholds for sensitivity analysis
set.seed(134)   # For reproducibility
filtered_snpset = snpgdsLDpruning(snp_filtered_gds, ld.threshold = 0.2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)
str(filtered_snpset)
names(filtered_snpset)

## Get all selected snp id
filtered_snpset_id = unlist(unname(filtered_snpset))
filtered_snpset_id

# ---Principal Component Analysis (PCA)---

## Pruning SNPs with LD threshold = 0.2 to exclude highly correlated SNPs
filtered_pca <- snpgdsPCA(snp_filtered_gds, snp.id = filtered_snpset_id, num.thread = 2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)

filtered_pca_df <- data.frame(sample.id = filtered_pca$sample.id,
    EV1 = filtered_pca$eigenvect[,1],
    EV2 = filtered_pca$eigenvect[,2],
    stringsAsFactors = FALSE)
filtered_pca_df = filtered_pca_df |> mutate(
    phenotype = case_when(
        sample.id == "strain_188" ~ "188-S",
        sample.id == "strain_185" ~ "185-R",
        sample.id == "strain_189" ~ "189-R",
        sample.id == "strain_191" ~ "191-I",
        TRUE ~ "WT-S"
    )
) |> mutate(
    sample.id = str_replace(sample.id, "_", " ")
)
head(filtered_pca_df)

## Variance proportion (%)
filtered_pc_percent <- filtered_pca$varprop*100
head(round(filtered_pc_percent, 2))

## PCA plot using pruned SNPs
set.seed(123)
PCA_vcf_p = ggplot(
    data = filtered_pca_df,
    aes(
        x = EV1, y = EV2,
    )
) +
    geom_point(
        aes(color = factor(phenotype)),
        position = position_jitter(width = 0.15, height = 0.15),
        size = 6, alpha = 1
    ) +
    # xlim(c(-0.8, 0.9)) +
    scale_x_continuous(limits = c(-0.75, 0.8),
                       breaks = seq(-0.75, 0.8, by = 0.5)) +
    labs(
        x = "PC1 (62.5%)", y = "PC2 (37.5%)",
    ) +
    theme_classic() +
    scale_color_npg() +
    scale_fill_npg() +
    theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "top",
        # legend.justification = c(1,1)
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)  # Adjust the size of legend point label
        )
    )
# ggsave("PCA_prunedSNPs.pdf", plot = PCA_vcf_p, units = "in", width = 6, height = 5, dpi = 300, device = "pdf")

#---Correlation matrix heatmap---#

## Parsing data
custom_color = list(
    strain = c(S191 = "#3C5488FF", S185 = "#DC0000FF", S188 = "#4DBBD5FF", S189 = "#00A087FF", S48 = "#F39B7FFF"),
    phenotype = c('5FC-R' = "#95D2B3", '5FC-S' = "#FBA834", 'WT-S' = "#FF76CE", '5FC-I' = "#45ba36")
)
pairwise_snp_dist = read.table("/mnt/hdd/dminh/Cauris/SNPs/snpdists/pairwise_snp_dist.tsv", sep = "\t", header = T, stringsAsFactors = F, check.names = F, row.names = 1)
cor_dist_matrix = pairwise_snp_dist |>
    filter(!row_number() %in% c(1)) |>
    select(!reference) |>
    cor(method = "pearson") |>
    data.frame() |>
    rename_with(~ str_replace(., "strain_", "S"))
rownames(cor_dist_matrix) = c("S185", "S188", "S189", "S191", "S48")

## Annotation label
group = data.frame(
    strain = c("S185", "S188", "S189", "S191", "S48"),
    phenotype = c("5FC-R", "5FC-S", "5FC-R", "5FC-I", "WT-S")
)
rownames(group) = colnames(cor_dist_matrix)

## callback function
callback = function(hc, mat){
  sv = svd(t(mat))$v[,c(1)]
  dend = reorder(as.dendrogram(hc), wts = sv^2)
  as.hclust(dend)
}

## Plot pairwise distance matrix heatmap
# jpeg("SNP_profile_distance_matrix_heatmap.jpg", width=10, height=8, unit="in", res=300)
heatmap_vcf = pheatmap(
    cor_dist_matrix,
    clustering_callback = callback,
    annotation_col = group,
    annotation_colors = custom_color,
    fontsize = 14
)
# dev.off()
