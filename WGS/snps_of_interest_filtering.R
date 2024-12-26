library(tidyverse)

strains = c("185", "186", "187", "188", "189", "190", "191", "192", "48", "48-YPD")
GT_transform = function(GT_string, strains) {
    GTs = strsplit(GT_string, ",")[[1]]
    specific_strains = strains[GTs == 1]
    paste(specific_strains, collapse = ",")
}

# ALT allele cut-off 40% in AD
snp_alt40_df = read_tsv("WGS/alt_AD_40/results/cauris.wgs.snps.of.interest.single-ann.tsv")
snp_alt40_df %>%
    select(`GEN[*].GT`, `CHROM:POS`, `ANN[0].GENE`, `ANN[0].EFFECT`, `ANN[0].HGVS_C`, `ANN[0].HGVS_P`) %>%
    rename(
        Strain = `GEN[*].GT`,
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`
    ) %>%
    mutate(
        Strain = sapply(Strain, GT_transform, strains = strains),
        Effect = recode(
            Effect,
            upstream_gene_variant = "Upstream",
            missense_variant = "Missense",
            downstream_gene_variant = "Downstream",
            stop_gained = "Stop gained",
            synonymous_variant = "Synonymous"
        )
    ) %>%
    write_tsv("WGS/alt_AD_40/results/cauris.wgs.snps.of.interest.alt40.tsv")

## CDS SNPs in S48 vs B8441
s48_snps_alt40_df = read_tsv("WGS/alt_AD_40/results/S48.snps.ann.tsv")
s48_snps_alt40_df %>%
    select(`GEN[0].GT`, `CHROM:POS`, `ANN[0].GENE`, FILTER, `ANN[0].EFFECT`, `ANN[0].HGVS_C`, `ANN[0].HGVS_P`, `GEN[0].AD`) %>%
    rename(
        Strain = `GEN[0].GT`,
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        Allelic_depths = `GEN[0].AD`
    ) %>%
    filter(
        !Effect %in% c("downstream_gene_variant", "upstream_gene_variant", "synonymous_variant")
    ) %>%
    mutate(
        Ref = "B8441",
        Strain = as.character(if_else(Strain == 1, 48, Strain))
    ) %>%
    select(
        Strain, Ref, `Scaffold:Pos`, GeneID, FILTER, Effect, cDNA, aa, Allelic_depths
    ) %>%
    write_tsv("WGS/alt_AD_40/results/S48_vs_B8441.snps.CDS.PASS.alt40.ann.tsv")

# ALT allele cut-off 80% in AD
snp_alt80_df = read_tsv("WGS/alt_AD_80/results/cauris.wgs.snps.of.interest.single-ann.tsv")
snp_alt80_df %>%
    select(`GEN[*].GT`, `CHROM:POS`, `ANN[0].GENE`, `ANN[0].EFFECT`, `ANN[0].HGVS_C`, `ANN[0].HGVS_P`) %>%
    rename(
        Strain = `GEN[*].GT`,
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`
    ) %>%
    mutate(
        Strain = sapply(Strain, GT_transform, strains = strains),
        Effect = recode(
            Effect,
            upstream_gene_variant = "Upstream",
            missense_variant = "Missense",
            downstream_gene_variant = "Downstream",
            stop_gained = "Stop gained",
            synonymous_variant = "Synonymous"
        )
    ) %>%
    write_tsv("WGS/alt_AD_80/results/cauris.wgs.snps.of.interest.alt80.tsv")

## CDS SNPs in S48 vs B8441
s48_snps_alt80_df = read_tsv("WGS/alt_AD_80/results/S48.snps.ann.tsv")
s48_snps_alt80_df %>%
    select(`GEN[0].GT`, `CHROM:POS`, `ANN[0].GENE`, FILTER, `ANN[0].EFFECT`, `ANN[0].HGVS_C`, `ANN[0].HGVS_P`, `GEN[0].AD`) %>%
    rename(
        Strain = `GEN[0].GT`,
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        Allelic_depths = `GEN[0].AD`
    ) %>%
    filter(
        !Effect %in% c("downstream_gene_variant", "upstream_gene_variant", "synonymous_variant")
    ) %>%
    mutate(
        Ref = "B8441",
        Strain = as.character(if_else(Strain == 1, 48, Strain))
    ) %>%
    select(
        Strain, Ref, `Scaffold:Pos`, GeneID, FILTER, Effect, cDNA, aa, Allelic_depths
    ) %>%
    write_tsv("WGS/alt_AD_80/results/S48_vs_B8441.snps.CDS.PASS.alt80.ann.tsv")

# S48 vs B8441 filtering on CDS SNPs
read_tsv("Supplementary_Figure_4/S48_B8441_WGS_RNA/S48.RNA.snps.ann.tsv") %>%
    select(`GEN[0].GT`, `CHROM:POS`, `ANN[0].GENE`, FILTER, `ANN[0].EFFECT`, `ANN[0].HGVS_C`, `ANN[0].HGVS_P`, `GEN[0].AD`) %>%
    rename(
        Strain = `GEN[0].GT`,
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        Allelic_depths = `GEN[0].AD`
    ) %>%
    filter(
        !Effect %in% c("downstream_gene_variant", "upstream_gene_variant", "synonymous_variant")
    ) %>%
    mutate(
        Ref = "B8441",
        Datasets = "RNA-seq",
        Strain = as.character(if_else(Strain == 1, 48, Strain))
    ) %>%
    select(
        Strain, Ref, `Scaffold:Pos`, GeneID, FILTER, Effect, cDNA, aa, Allelic_depths, Datasets
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_B8441_WGS_RNA/S48-B8441_CDS_snps_RNA.tsv")

# S48-YPD vs B8441 filtering on CDS SNPs
read_tsv("Supplementary_Figure_4/S48_S48-YPD_WGS/S48-YPD.WGS.snps.ann.tsv") %>%
    select(`GEN[0].GT`, `CHROM:POS`, `ANN[0].GENE`, FILTER, `ANN[0].EFFECT`, `ANN[0].HGVS_C`, `ANN[0].HGVS_P`, `GEN[0].AD`) %>%
    rename(
        Strain = `GEN[0].GT`,
        "Scaffold:Pos" = `CHROM:POS`,
        GeneID = `ANN[0].GENE`,
        Effect = `ANN[0].EFFECT`,
        cDNA = `ANN[0].HGVS_C`,
        aa = `ANN[0].HGVS_P`,
        Allelic_depths = `GEN[0].AD`
    ) %>%
    filter(
        !Effect %in% c("downstream_gene_variant", "upstream_gene_variant", "synonymous_variant")
    ) %>%
    mutate(
        Ref = "B8441",
        Datasets = "WGS",
        Strain = as.character(if_else(Strain == 1, "48YPD-178", "MismatchGT"))
    ) %>%
    select(
        Strain, Ref, `Scaffold:Pos`, GeneID, FILTER, Effect, cDNA, aa, Allelic_depths, Datasets
    ) %>%
    write_tsv("Supplementary_Figure_4/S48_S48-YPD_WGS/S48YPD-B8441_CDS_snps_WGS.tsv")
