library(VennDiagram)
library(tidyverse)
library(patchwork)

rna_calling = read_tsv("Supplementary_Figure_4/S48_B8441_WGS_RNA/S48.RNA.snps.ann.tsv")
wgs_calling = read_tsv("Supplementary_Figure_4/S48_B8441_WGS_RNA/S48.WGS.snps.ann.tsv")

# Overlapped variants: 660
overlapped_vars_pos = inner_join(rna_calling, wgs_calling, by = "CHROM:POS") %>% pull("CHROM:POS")
wgs_calling %>%
    filter(`CHROM:POS` %in% overlapped_vars_pos) %>%
    rename(
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        Impact = `ANN[0].IMPACT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        GT = `GEN[0].GT`,
        AD = `GEN[0].AD`
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_B8441_WGS_RNA/overlapped-variants_S48-B8441_WGS-RNA.tsv")

# Non-overlapped variants in RNA calling without a match in WGS calling: 7
non_overlapped_vars_pos_RNA = anti_join(rna_calling, wgs_calling, by = "CHROM:POS") %>% pull("CHROM:POS")
rna_calling %>%
    filter(`CHROM:POS` %in% non_overlapped_vars_pos_RNA) %>%
    rename(
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        Impact = `ANN[0].IMPACT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        GT = `GEN[0].GT`,
        AD = `GEN[0].AD`
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_B8441_WGS_RNA/non-overlapped-variants_S48-B8441_RNA.tsv")

# Non-overlapped variants in WGS calling without a match in RNA calling: 309
non_overlapped_vars_pos_wgs = anti_join(wgs_calling, rna_calling, by = "CHROM:POS") %>% pull("CHROM:POS")
wgs_calling %>%
    filter(`CHROM:POS` %in% non_overlapped_vars_pos_wgs) %>%
    rename(
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        Impact = `ANN[0].IMPACT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        GT = `GEN[0].GT`,
        AD = `GEN[0].AD`
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_B8441_WGS_RNA/non-overlapped-variants_S48-B8441_WGS.tsv")

# Create venn object
var_lst = list(
    RNA = rna_calling %>% pull("CHROM:POS"),
    WGS = wgs_calling %>% pull("CHROM:POS")
)

# Proportional venn plot with hypergeometric test
total_pop_size = full_join(rna_calling, wgs_calling, by = "CHROM:POS") %>% nrow()
venn.diagram(
    var_lst,
    fill = c("#9467bd", "#2ca02c"),
    alpha = c(0.6, 0.6),
    lwd = 0,
    sub.cex = 1.5,
    cat.fontface = "bold", cat.cex = c(1.5, 1.5),
    cex = 2,
    hyper.test = TRUE, total.population = total_pop_size,
    sub.pos = c(0.5, 1.0),
    # filename = "Supplementary_Figure_4/S48_B8441_WGS_RNA/S48-B8441_WGS-RNA_alt80_venn.tiff", compression = "lzw"
    filename = "Supplementary_Figure_4/S48_B8441_WGS_RNA/S48-B8441_WGS-RNA_alt80_venn.png", imagetype = "png"     # PNG format
)

#--------------------------------------------------------#
# SNPs in S48 vs S48-YPD (WGS)
s48_YPD_control_df = read_tsv("Supplementary_Figure_4/S48_S48-YPD_WGS/S48-YPD.WGS.snps.ann.tsv")

## Overlapped SNPs between S48 and S48-YPD
overlapped_control_snps = inner_join(s48_YPD_control_df, wgs_calling, by = "CHROM:POS") %>% pull("CHROM:POS")
s48_YPD_control_df %>%
    filter(`CHROM:POS` %in% overlapped_control_snps) %>%
    rename(
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        Impact = `ANN[0].IMPACT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        GT = `GEN[0].GT`,
        AD = `GEN[0].AD`
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_S48-YPD_WGS/overlapped-snps_S48-S48YPD_WGS.tsv")

## Non-overlapped SNPs in S48-YPD without a match in S48 (WGS)
non_overlapped_S48_YPD_snps = anti_join(s48_YPD_control_df, wgs_calling, by = "CHROM:POS") %>% pull("CHROM:POS")
s48_YPD_control_df %>%
    filter(`CHROM:POS` %in% non_overlapped_S48_YPD_snps) %>%
    rename(
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        Impact = `ANN[0].IMPACT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        GT = `GEN[0].GT`,
        AD = `GEN[0].AD`
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_S48-YPD_WGS/non-overlapped-snps_S48-YPD_WGS.tsv")

## Non-overlapped SNPs in S48 without a match in S48-YPD (WGS)
non_overlapped_S48_snps = anti_join(wgs_calling, s48_YPD_control_df, by = "CHROM:POS") %>% pull("CHROM:POS")
wgs_calling %>%
    filter(`CHROM:POS` %in% non_overlapped_S48_snps) %>%
    rename(
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        Impact = `ANN[0].IMPACT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        GT = `GEN[0].GT`,
        AD = `GEN[0].AD`
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_S48-YPD_WGS/non-overlapped-snps_S48_WGS.tsv")

# Create venn object
snps_control_lst = list(
    S48 = wgs_calling %>% pull("CHROM:POS"),
    "S48-YPD" = s48_YPD_control_df %>% pull("CHROM:POS")
)

# Proportional venn plot with hypergeometric test
total_pop_size_control_snp = full_join(wgs_calling, s48_YPD_control_df, by = "CHROM:POS") %>% nrow()
venn.diagram(
    snps_control_lst,
    fill = c("#2ca02c", "#9467bd"),
    alpha = c(0.6, 0.6),
    lwd = 0,
    sub.cex = 1.5,
    cat.fontface = "bold", cat.cex = c(1.5, 1.5),
    cex = 2,
    hyper.test = TRUE, total.population = total_pop_size_control_snp,
    sub.pos = c(0.5, 1.02),
    cat.pos = c(226, 136),
    filename = "Supplementary_Figure_4/S48_S48-YPD_WGS/S48-YPD_WGS_alt80_venn.tiff", compression = "lzw"
    # filename = "Supplementary_Figure_4/S48_S48-YPD_WGS/S48-YPD_WGS_alt80_venn.png", imagetype = "png"     # PNG format
)

