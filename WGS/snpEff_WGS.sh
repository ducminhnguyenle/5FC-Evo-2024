#!/usr/bin/env bash

# 80% alt reads in AD
vcf="WGS/finalfiltered_alt80.vcf.gz"
norm="WGS/finalfiltered.norm.vcf.gz"
ann="WGS/finalfiltered.norm.scaffold.ann.vcf"

# Decomposing and Normalizing vcf files
bcftools norm -m- "${vcf}" | awk '{if(/^#/) {print $0} else if($5 !~ /*/) {print $0} }' | bgzip > "${norm}"
tabix -p vcf "${norm}"

# Rename contigs to match annotation database
bcftools annotate --rename-chrs \
    "WGS/Cauris_contigs_rev_rename.txt" \
    "${norm}" | bgzip > "${norm/vcf.gz/scaffold.vcf.gz}"
tabix -p vcf "WGS/finalfiltered.norm.scaffold.vcf.gz"

# Annotation using snpEff, using this annotation database "_candida_auris_gca_002759435" to match the reference genome
## Check available annotation databases for C.auris
snpEff databases | grep -i "candida_auris"
## Running snpEff
snpEff \
    -v \
    -stats "${ann/vcf/html}" \
    "_candida_auris_gca_002759435" \
    "WGS/finalfiltered.norm.scaffold.vcf.gz" \
    > "${ann}"

# Remap contigs back to original to match reference genome
bcftools annotate --rename-chrs \
    "WGS/Cauris_contigs_ori_rename.txt" \
    "${ann}" | bgzip > "${ann/scaffold/PEKT}.gz"
tabix -p vcf "WGS/finalfiltered.norm.PEKT.ann.vcf.gz"

# Filtering variants that PASS hard filtering and overall DP > 30
bcftools view -H "WGS/finalfiltered.norm.PEKT.ann.vcf.gz" | wc -l        # 1067 variants
bcftools filter -i 'INFO/DP >= 30 & FILTER == "PASS"' "WGS/finalfiltered.norm.PEKT.ann.vcf.gz" | grep -v -c "^#"    # 1046

# For loop for all pairwise comparisons against strain_48 (multi-samples vcf is filtered with "PASS" hard filtering, overall DP >= 30 and masked genotype filtering)
strains=("S185" "S186" "S187" "S188" "S189" "S190" "S191" "S192")
for strain in "${strains[@]}"; do
    echo "${strain}"
    bcftools view -s "${strain},S48" "WGS/finalfiltered.norm.PEKT.ann.vcf.gz" | bcftools filter -i 'GT="0"' - | bcftools query -f 'GT=[%GT]\n' - | sort | uniq -c
    bcftools view -s "${strain},S48" "WGS/finalfiltered.norm.PEKT.ann.vcf.gz" | \
    bcftools filter -i 'GT="0"' - | \
    awk '{
        if (/^#/) {
            print $0
        }
        else if ((substr($10,1,1) != substr($11,1,1)) && (substr($10,1,1) != "\.") && (substr($11,1,1) != "\.")) {
            print $0
        }
    }' > "WGS/pairs/cauris_wgs_${strain}_vs_S48_filtered_PEKT.ann.vcf"
done

# Subset regions of interest
for file in WGS/pairs/*_S48*.ann.vcf; do
    bcftools view -H "${file}" | awk 'BEGIN{FS=OFS="\t"}; {print $1,$2}' >> "WGS/pairs/regions.tsv"
done

# Selected SNPs
bcftools view -T "WGS/alt_AD_80/pairs/regions.tsv" "WGS/alt_AD_80/norm_ann/finalfiltered_alt80.norm.PEKT.ann.vcf.gz" > "WGS/alt_AD_80/results/cauris.wgs.snps.of.interest.ann.vcf"
bcftools view -T "WGS/alt_AD_80/pairs/regions.tsv" "WGS/alt_AD_80/norm_ann/finalfiltered_alt80.norm.PEKT.ann.vcf.gz" | grep -v "^##" | awk 'BEGIN{FS=OFS="\t"};{sub(/#/, "", $1); print $1,$2,$4,$5,$6,$7,$10,$11,$12,$13,$14,$15,$16,$17,$18}' > "WGS/alt_AD_80/results/cauris.wgs.snps.of.interest.tsv"

# Parsing snpeff output
SnpSift extractFields -s "," -e "." \
    "WGS/alt_AD_80/results/cauris.wgs.snps.of.interest.ann.vcf" \
    CHROM POS REF ALT "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" "GEN[*].GT" | \
    awk 'BEGIN{FS=OFS="\t"}; {sub("\\.1", "", $1); print $1":"$2,$3,$4,$5,$6,$7,$8,$9,$10}' \
    > "WGS/alt_AD_80/results/cauris.wgs.snps.of.interest.single-ann.tsv"

# S187
bcftools view -s "S187,S48" "WGS/finalfiltered.norm.PEKT.ann.vcf.gz" | \
    bcftools filter -i 'GT="0"' - | \
    awk '{
        if (/^#/) {
            print $0
        }
        else if (substr($10,1,1) != substr($11,1,1)) {
            print $0
        }
    }' | less -S

# S48 (Wild type) vs B8441
B9J08_000964
B9J08_001270

bcftools view -s S48 WGS/alt_AD_80/norm_ann/finalfiltered_alt80.norm.PEKT.ann.vcf.gz | bcftools filter -i 'GT="1"' - | bcftools view -i 'FILTER="PASS"' - > "WGS/alt_AD_80/results/S48.snps.ann.vcf"

SnpSift extractFields -s "," -e "." \
    "WGS/alt_AD_80/results/S48.snps.ann.vcf" \
    "CHROM" "POS" "REF" "ALT" "FILTER" "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" "GEN[0].GT" "GEN[0].AD" | \
    awk 'BEGIN{FS=OFS="\t"}; {sub("\\.1", "", $1); print $1":"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
    > "WGS/alt_AD_80/results/S48.snps.ann.tsv"

# S48-YPD vs B8441
bcftools view -s S48-YPD WGS/alt_AD_80/norm_ann/finalfiltered_alt80.norm.PEKT.ann.vcf.gz |\
    bcftools filter -i 'GT="1" & FILTER == "PASS"' - > "WGS/alt_AD_80/S48-YPD_vs_B8441/S48-YPD.WGS.snps.ann.vcf"

SnpSift extractFields -s "," -e "." \
    "WGS/alt_AD_80/S48-YPD_vs_B8441/S48-YPD.WGS.snps.ann.vcf" \
    "CHROM" "POS" "REF" "ALT" "FILTER" "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" "GEN[0].GT" "GEN[0].AD" | \
    awk 'BEGIN{FS=OFS="\t"}; {sub("\\.1", "", $1); print $1":"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
    > "WGS/alt_AD_80/S48-YPD_vs_B8441/S48-YPD.WGS.snps.ann.tsv"

# bedtools intersect \
#     -a "WGS/alt_AD_80/results/S48.snps.ann.vcf" \
#     -b "WGS/alt_AD_80/S48-YPD_vs_B8441/S48-YPD.WGS.snps.ann.vcf" | wc -l    # 967
