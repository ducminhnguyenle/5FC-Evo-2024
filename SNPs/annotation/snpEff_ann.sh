#!/usr/bin/env bash
#----------------------------------------------------------#
vcf="Cauris_merged_output/HaplotypeCaller/selectedsnpsfiltered"
ann="Cauris_merged_output/Annotation"
norm="Cauris_merged_output/HaplotypeCaller/normalizing"

# Clean sample's name in vcf files
bcftools reheader \
    -s "Cauris_merged_output/HaplotypeCaller/selectedsnpsfiltered/samples_rename.txt" \
    "${vcf}/cauris.combined.genotype.snps.filtered.vcf.gz" \
    > "${vcf}/cauris.combined.genotype.snps.filtered.rename.vcf.gz"
tabix -p vcf "${vcf}/cauris.combined.genotype.snps.filtered.rename.vcf.gz"

# Decomposing and Normalizing vcf files
bcftools norm -m- "${vcf}/cauris.combined.genotype.snps.filtered.rename.vcf.gz" | awk '{if(/^#/) {print $0} else if($5 !~ /*/) {print $0} }' | bgzip > "${norm}/cauris.combined.genotype.snps.filtered.norm.vcf.gz"
tabix -p vcf "${norm}/cauris.combined.genotype.snps.filtered.norm.vcf.gz"

# Rename contigs to match annotation database
bcftools annotate --rename-chrs \
    "Cauris_merged_output/Annotation/Cauris_contigs_rev_rename.txt" \
    "${norm}/cauris.combined.genotype.snps.filtered.norm.vcf.gz" | bgzip > "${ann}/cauris.snp.filtered.scaffold.vcf.gz"
tabix -p vcf "${ann}/cauris.snp.filtered.scaffold.vcf.gz"

# Annotate using snpEff annotation database "_candida_auris_gca_002759435" to match the reference genome
## Check available annotation databases for C.auris
snpEff databases | grep -i "candida_auris"
## Running snpEff
snpEff \
    -v \
    -stats "${ann}/cauris.snp.filtered.scaffold.ann.html" \
    "_candida_auris_gca_002759435" \
    "${ann}/cauris.snp.filtered.scaffold.vcf.gz" \
    > "${ann}/cauris.snp.filtered.scaffold.ann.vcf"

# Remap contigs back to original to match reference genome
bcftools annotate --rename-chrs \
    "${ann}/Cauris_contigs_ori_rename.txt" \
    "${ann}/cauris.snp.filtered.scaffold.ann.vcf" | bgzip > "${ann}/cauris.snp.filtered.PEKT.ann.vcf.gz"
tabix -p vcf "${ann}/cauris.snp.filtered.PEKT.ann.vcf.gz"

# Filtering variants that PASS hard filtering and overall DP > 30
bcftools view -H cauris.snp.filtered.PEKT.ann.vcf.gz | wc -l
bcftools filter -i 'INFO/DP >= 30 & FILTER == "PASS"' "./SNPs/annotation/marked_filtered/cauris.snp.filtered.PEKT.ann.vcf.gz" > "./SNPs/annotation/marked_filtered/filtered/cauris.snps.PASS.DP30.PEKT.ann.vcf"

# Genotype filtering (based on FORMAT column - masking the GT)
# gzip -c -d "Cauris_merged_output/Annotation/cauris.snp.filtered.PEKT.ann.vcf.gz" > "Cauris_merged_output/Annotation/cauris.snp.filtered.PEKT.ann.vcf"

filterGatkGenotypes="./bin/filterGatkGenotypes.py"
python3 ${filterGatkGenotypes} "./SNPs/annotation/marked_filtered/filtered/cauris.snps.PASS.DP30.PEKT.ann.vcf" \
                                --min_GQ "50" \
                                --keep_GQ_0_refs \
                                --min_percent_alt_in_AD "0.8" \
                                --min_total_DP "30" \
                                --keep_all_ref \
                                > "./SNPs/annotation/marked_filtered/filtered/cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.vcf"

# For loop for all pairwise comparisons against strain_48 (multi-samples vcf is filtered with "PASS" hard filtering, overall DP >= 30 and masked genotype filtering)
strains=("strain_185" "strain_188" "strain_189" "strain_191")
for strain in "${strains[@]}"; do
    echo "${strain}"
    type="${strain/strain_/}"
    bcftools view -s "${strain},strain_48" "filtered/cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.vcf" | bcftools filter -i 'GT="0"' - | bcftools query -f 'GT=[%GT]\n' - | sort | uniq -c
    bcftools view -s "${strain},strain_48" "filtered/cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.vcf" | \
    bcftools filter -i 'GT="0"' - | \
    awk '{
        if (/^#/) {
            print $0
        }
        else if ((substr($10,1,1) != substr($11,1,1)) && (substr($10,1,1) != "\.") && (substr($11,1,1) != "\.")) {
            print $0
        }
    }' > "./filtered/cauris_${type}_vs_48_filtered_PEKT.ann.vcf"
done

# Subset regions of interest
for file in ./*_48*.ann.vcf; do
    bcftools view -H "${file}" | awk 'BEGIN{FS=OFS="\t"}; {print $1,$2}' >> regions.tsv
done

# Selected SNPs
bcftools view -T regions.tsv "cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.vcf" | grep -v "^##" | awk 'BEGIN{FS=OFS="\t"};{sub(/#/, "", $1); print $1,$2,$4,$5,$6,$7,$10,$11,$12,$13,$14}' > "cauris.snps.of.interest.tsv"