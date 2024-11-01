library(pheatmap)
library(tidyverse)
library(patchwork)
library(ggplotify)
# https://stackoverflow.com/questions/73872992/pheatmap-manually-re-order-leaves-in-dendogram

setwd("/mnt/hdd/dminh/Cauris/SNPs/annotation/new/filtered/")
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
# Annotation label
group = data.frame(
    strain = c("S185", "S188", "S189", "S191", "S48"),
    phenotype = c("5FC-R", "5FC-S", "5FC-R", "5FC-I", "WT-S")
)
rownames(group) = colnames(cor_dist_matrix)

# callback function
callback = function(hc, mat){
  sv = svd(t(mat))$v[,c(1)]
  dend = reorder(as.dendrogram(hc), wts = sv^2)
  as.hclust(dend)
}

# Pairwise distance matrix heatmap
# jpeg("SNP_profile_distance_matrix_heatmap.jpg", width=10, height=8, unit="in", res=300)
heatmap_vcf = pheatmap(
    cor_dist_matrix,
    clustering_callback = callback,
    annotation_col = group,
    annotation_colors = custom_color,
    fontsize = 14
)
ggsave("/mnt/hdd/dminh/Cauris/SNPs/vcf_QC/SNP_profile_distance_matrix_heatmap.p", plot = heatmap_vcf, units = "in", width = 12, height = 8, dpi = 300, device = "pdf")
# dev.off()
