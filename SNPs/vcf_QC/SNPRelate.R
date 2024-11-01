# https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(ggforce)
library(ggrepel)
library(ggplotify)

setwd("/mnt/hdd/dminh/Cauris/SNPs/annotation/marked_filtered/filtered/")
# Input multiple samples vcf file
## Using "PASS" hard filtering and overall DP >= 30
snp_filtered_vcf = "cauris.snps.PASS.DP30.PEKT.ann.vcf"

# Reformat from vcf to gds object
snpgdsVCF2GDS(snp_filtered_vcf, "cauris.snps.PASS.DP30.PEKT.ann.gds", method = "biallelic.only")

# Summary
snpgdsSummary("cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.gds")     # 5 samples with 926 SNPs with masked genotype

# Open a GDS file
snp_filtered_gds = snpgdsOpen("cauris.snps.PASS.DP30.PEKT.ann.gds")
# snp_geno_filtered_gds = snpgdsOpen("cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.gds")
# snp_subset_filtered_gds = snpgdsOpen("test/cauris.subset.snps.PASS.DP30.PEKT.ann.gds")

## 1. Inspect all samples in vcf file
head(read.gdsn(index.gdsn(snp_filtered_gds, "sample.id")))
## 2. Inspect pos of all variants
head(read.gdsn(index.gdsn(snp_filtered_gds, "snp.position")))
## 3. Inspect the contigs/scaffolds/chr of each variant
head(read.gdsn(index.gdsn(snp_filtered_gds, "snp.chromosome")))
## 4. Inspect genotype data
# Haploid: rows -> samples, columns -> variants
# vcf: 0 (REF) -> gds: 1, vcf: 1 (ALT) -> gds: 0, vcf: missing -> gds: 3
read.gdsn(index.gdsn(snp_filtered_gds, "genotype"), start = c(1,1), count = c(5, 100))

# ---LD-based SNP pruning---

## Different LD thresholds for sensitivity analysis
set.seed(134)
filtered_snpset = snpgdsLDpruning(snp_filtered_gds, ld.threshold = 0.2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)
str(filtered_snpset)
names(filtered_snpset)

## Get all selected snp id
filtered_snpset_id = unlist(unname(filtered_snpset))
filtered_snpset_id

# ---Principal Component Analysis (PCA)---

## Run PCA
filtered_pca <- snpgdsPCA(snp_filtered_gds, snp.id = filtered_snpset_id, num.thread = 2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)
# filtered_pca <- snpgdsPCA(snp_geno_filtered_gds, num.thread = 2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)
# filtered_pca_all <- snpgdsPCA(snp_filtered_gds, num.thread = 2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)       # all 926 SNPs that pass hard filtering and overall DP >= 30
# filtered_subset_pca <- snpgdsPCA(snp_subset_filtered_gds, num.thread = 2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)

## Variance proportion (%)
filtered_pc_percent <- filtered_pca$varprop*100
head(round(filtered_pc_percent, 2))

## Make a data.frame
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
filtered_pca_df2 = filtered_pca_df |> filter(
    phenotype != "5FC-R"
)
## Draw simple plot
# plot(filtered_pca_df$EV2, filtered_pca_df$EV1, xlab = "eigenvector 2", ylab = "eigenvector 1")

## ggplot2
set.seed(123)
ggplot(
    data = filtered_pca_df,
    aes(
        x = EV1, y = EV2,
        # color = sample.id
    )
) +
    geom_point(
        aes(color = factor(phenotype)),
        position = position_jitter(width = 0.1, height = 0.1),
        size = 12, alpha = 0.8
    ) +
    # geom_mark_ellipse(
    #     aes(
    #         fill = phenotype,
    #         label = phenotype,
    #         filter = phenotype == "5FC-R"
    #     ),
    #     expand = unit(14, "mm"),
    #     show.legend = FALSE,
    #     label.fontsize = 14
    # ) +
    # xlim(c(-0.5, 1)) +
    # ylim(c(-0.75, 0.5)) +
    # geom_label_repel(
    #     data = filtered_pca_df2,
    #     aes(
    #         label = phenotype,
    #         fontface = "bold",
    #         size = 14
    #     ),
    #     min.segment.length = Inf,
    #     box.padding = 0.4,
    #     show.legend = FALSE,
    #     nudge_x = 0.1,
    #     direction = "y"
    # ) +
    labs(
        title = "Candida auris 5FC-resistance strains\nare separated from others", x = "PC1 (62.5%)", y = "PC2 (37.5%)",
    ) +
    theme_classic() +
    scale_color_manual(values = c("#cc4ee0", "#9e82ed", "#e89829", "#82ed82", "#5977ff")) +
    # scale_color_npg() +
    # scale_fill_npg() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "American Typewriter", lineheight = 1.2 ),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
    )
# ggsave("PCA_SNP_profiles.tiff", unit = "in", width = 12, height = 8, device = "tiff", dpi = 300, compression = "lzw")

# Plot the principal component pairs for the first four PCs:

## Get sample id
sample.id <- read.gdsn(index.gdsn(filtered_gds, "sample.id"))

## Get group information
group_code <- scan("/mnt/hdd/dminh/Cauris/SNPs/annotation/marked_filtered/filtered/group.txt", what=character())
group = as.data.frame(cbind(sample.id, group_code))

## Merge group with filtered_pca_df
merged_df = merge(filtered_pca_df, group, by = "sample.id")

## Plot the principal component pairs for the first four PCs:
lbls = paste("PC", 1:4, "\n", format(filtered_pc_percent[1:4], digits=2), "%", sep="")
pairs(filtered_pca$eigenvect[,1:4], col=merged_df$group_code, labels=lbls)

# # To calculate the SNP correlations between eigenvectors and SNP genotypes:

# ## Get scaffold index
# scaf <- read.gdsn(index.gdsn(filtered_gds, "snp.chromosome"))
# scaf2 = gsub("PEKT020000|\\.1", "", scaf)
# CORR <- snpgdsPCACorr(filtered_pca, filtered_gds, eig.which=1:4)

# ## Plot
# savepar <- par(mfrow=c(2,1), mai=c(0.45, 0.55, 0.1, 0.25))
# for (i in 1:2) {
#     plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i), col = scaf2,
#     pch="+")
# }
# par(savepar)

# test = CORR$snpcorr[1, ] |>
#     as.data.frame() |>
#     rownames_to_column() |>
#     cbind(scaf)

# ggplot(data = test, aes(x = as.numeric(rowname), y = CORR$snpcorr[1, ])) +
#     geom_point(aes(color = factor(scaf)), size = 4, alpha = 1, shape = "+") +
#     ylim(c(0, 1))
