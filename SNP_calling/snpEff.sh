vcf="/mnt/rdisk/dminh/Cauris/Cauris_merged_output/HaplotypeCaller/selectedsnpsfiltered"
ann="/mnt/rdisk/dminh/Cauris/Cauris_merged_output/Annotation"
# Clean sample's name in vcf files
bcftools reheader \
    -s "/mnt/rdisk/dminh/Cauris/Cauris_merged_output/HaplotypeCaller/selectedsnpsfiltered/samples_rename.txt" \
    "${vcf}/cauris.combined.genotype.snps.filtered.vcf.gz" \
    > "${vcf}/cauris.combined.genotype.snps.filtered.rename.vcf.gz"
tabix -p vcf "${vcf}/cauris.combined.genotype.snps.filtered.rename.vcf.gz"

# Decomposing and Normalizing vcf files
norm="/mnt/rdisk/dminh/Cauris/Cauris_merged_output/HaplotypeCaller/normalizing"
bcftools norm -m- "${vcf}/cauris.combined.genotype.snps.filtered.rename.vcf.gz" | awk '{if(/^#/) {print $0} else if($5 !~ /*/) {print $0} }' | bgzip > "${norm}/cauris.combined.genotype.snps.filtered.norm.vcf.gz"
tabix -p vcf "${norm}/cauris.combined.genotype.snps.filtered.norm.vcf.gz"

# Rename contigs to match annotation database
bcftools annotate --rename-chrs \
    "/mnt/rdisk/dminh/Cauris/Cauris_merged_output/Annotation/Cauris_contigs_rev_rename.txt" \
    "${norm}/cauris.combined.genotype.snps.filtered.norm.vcf.gz" | bgzip > "${ann}/cauris.snp.filtered.scaffold.vcf.gz"
tabix -p vcf "${ann}/cauris.snp.filtered.scaffold.vcf.gz"

# Annotation using snpEff, using this annotation database "_candida_auris_gca_002759435" to match the reference genome
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
    "/mnt/rdisk/dminh/Cauris/Cauris_merged_output/Annotation/Cauris_contigs_ori_rename.txt" \
    "${ann}/cauris.snp.filtered.scaffold.ann.vcf" | bgzip > "${ann}/cauris.snp.filtered.PEKT.ann.vcf.gz"
tabix -p vcf "${ann}/cauris.snp.filtered.PEKT.ann.vcf.gz"

# Filtering variants that PASS hard filtering and overall DP > 30
bcftools view -H cauris.snp.filtered.PEKT.ann.vcf.gz | wc -l        # 1198 variants
bcftools filter -i 'INFO/DP >= 30 & FILTER == "PASS"' "cauris.snp.filtered.PEKT.ann.vcf.gz" | grep -v "^#" | wc -l    # 926
bcftools filter -i 'INFO/DP >= 30 & FILTER == "PASS"' "cauris.snp.filtered.PEKT.ann.vcf.gz" > "filtered/cauris.snps.PASS.DP30.PEKT.ann.vcf"

# FILTER_GATK_GENOTYPES (based on FORMAT) (by masking the GT)
filterGatkGenotypes="/mnt/rdisk/dminh/Cauris/bin/filterGatkGenotypes.py"
python3 ${filterGatkGenotypes}  "filtered/cauris.snps.PASS.DP30.PEKT.ann.vcf" \
                                --min_GQ "50" \
                                --keep_GQ_0_refs \
                                --min_percent_alt_in_AD "0.8" \
                                --min_total_DP "30" \
                                --keep_all_ref \
                                > "filtered/cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.vcf"

# All pairwise comparisons against strain_48
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
    }' > "./filtered/cauris_${type}_vs_48_finalfiltered_PEKT.ann.vcf"
done

# VCF_TO_FASTA (Multi-samples fasta)
MAX_PERC_AMB_SAMPLES=10
NUM_SAMPLES=$(cat "filtered/samplelist.txt" | wc -l)
if [ ${MAX_PERC_AMB_SAMPLES} > 0 ]; then
    MAX_AMB_SAMPLES=$(echo "${NUM_SAMPLES} 10" | awk '{x=$1*($2/100); y=int(x); x=(y<1?1:y)} END {print x}')
else
    MAX_AMB_SAMPLES=10000000
fi
## Note: filter with min_depth=30 is not good for consensus distance matrix
vcfSnpsToFasta="/mnt/rdisk/dminh/Cauris/bin/vcfSnpsToFasta.py"
gzip -c -d "filtered/cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.vcf.gz" > "filtered/cauris.snps.PASS.DP30.genotypefiltered.PEKT.ann.vcf"
python3 ${vcfSnpsToFasta} \
    --max_amb_samples $MAX_AMB_SAMPLES \
    --min_depth 0 \
    "filtered/cauris.snps.PASS.DP30.PEKT.ann.vcf" > "filtered/vcf_to_fasta/passvcf_DP30_to_fasta.fasta"

# SEQKIT_REPLACE
seqkit \
    replace \
    -s -p '\*' -r '-' \
    --threads 2 \
    -i "filtered/vcf_to_fasta/passvcf_DP30_to_fasta.fasta" \
    -o "filtered/vcf_to_fasta/multi_vcfs_to_fasta.fasta"

# SNPDISTS
snp-dists \
    -b \
    "filtered/vcf_to_fasta/multi_vcfs_to_fasta.fasta" > "filtered/vcf_to_fasta/pairwise_snp_dist.tsv"
