# --java-options "-Xmx4G"
# Prepare reference genome
REF_DIR="/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
samtools faidx "${REF_DIR}/genome.fa"
gatk CreateSequenceDictionary -R "${REF_DIR}/genome.fa"

# Handling Splicing Events in RNASeq Data using SplitNCigarReads: Split reads into exon and intron segments
mrkdup="/mnt/rdisk/dminh/Cauris/Cauris_output/MarkDuplicates"
split="/mnt/rdisk/dminh/Cauris/Cauris_output/SplitReads"
for bam in `ls ${mrkdup} | grep -vE "metrics|bai"`; do 
    echo "${bam}"
    gatk SplitNCigarReads \
        --reference "${REF_DIR}/genome.fa" \
        --input "${mrkdup}/${bam}" \
        --output "${split}/${bam/markDups/splitReads}"
done

# Adding read group (RG) information (Remember adding RG before SplitNCigarReads step, after alignment step)
split="/mnt/rdisk/dminh/Cauris/Cauris_output/SplitReads"
RG="/mnt/rdisk/dminh/Cauris/Cauris_output/SplitReads/RG"
for file in `ls ${split}/*.bam`; do
    bam=$(basename "${file}")
    sample=$(basename "${file}" .splitReads.bam)
    echo "${sample}"
    gatk --java-options "-Xmx128G" AddOrReplaceReadGroups \
        -I "${split}/${bam}" \
        -O "${RG}/${bam/.splitReads.bam/.splitReads.RG.bam}" \
        --RGLB "lib1" \
        --RGPL "Illumina" \
        --RGPU "unit1" \
        --RGSM "${sample}"
done

## Index splitreads with added RG bam file
RG="/mnt/rdisk/dminh/Cauris/Cauris_output/SplitReads/RG"
for bam in `ls ${RG}`; do
    echo "${RG}/${bam}"
    samtools index -@ 60 "${RG}/${bam}"
done

# Preparing snp known-sites for C.aruis
## Rename contigs in truth set and consensus snp vcf files from this paper "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10210944/"
cat genome.fa | grep "^>" | awk 'BEGIN{FS=OFS=" "}; {gsub(",|>", "", $0); print $6,$1}' > "Cauris_contigs_rename.txt"
bcftools annotate --rename-chrs Cauris_contigs_rename.txt B11221_CA05.vcf.gz | bgzip > "B11221_CA05_contigs_rename.vcf.gz"
bcftools annotate --rename-chrs Cauris_contigs_rename.txt B11245_CA06.vcf.gz | bgzip > "B11245_CA06_contigs_rename.vcf.gz"
bcftools annotate --rename-chrs Cauris_contigs_rename.txt cauris.consensus.vcf.gz | bgzip > "cauris.consensus_contigs_rename.vcf.gz"
## Decomposing, normalising multi-allelic variants and removing invalid variants
# bcftools norm -m- cauris.consensus_contigs_rename.vcf.gz | bgzip > "cauris.consensus_contigs_rename_split.vcf.gz"
# zcat cauris.consensus_contigs_rename_split.vcf.gz | awk '{if(/^#/) {print $0} else if($5 !~ /-/) {print $0} }' | bgzip > "cauris.consensus_contigs_rename_split_edited.vcf.gz"
bcftools norm -m- cauris.consensus_contigs_rename.vcf.gz | awk '{if(/^#/) {print $0} else if($5 !~ /-/) {print $0} }' | bgzip > "cauris.consensus_contigs_rename_split_edited.vcf.gz"
## Index all vcf files
tabix -p vcf "known-sites"

# BQSR using adjusted SNPs vcf files
## 1. BaseRecalibrator (First pass)
split_RG="/mnt/rdisk/dminh/Cauris/Cauris_output/SplitReads/RG"
BQSR="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recal_table"
for file in `ls ${split_RG}/*.bam`; do
    bam=$(basename "${file}")
    gatk --java-options "-Xmx128G" BaseRecalibrator \
        --reference "${REF_DIR}/genome.fa" \
        --input "${split_RG}/${bam}" \
        --output "${BQSR}/${bam/splitReads.RG.bam/recal.table}" \
        --known-sites "/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11221_CA05_contigs_rename.vcf.gz" \
        --known-sites "/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11245_CA06_contigs_rename.vcf.gz" \
        --known-sites "/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/cauris.consensus_contigs_rename_split_edited.vcf.gz"
done
## 2. ApplyBQSR
REF_DIR="/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
split_RG="/mnt/rdisk/dminh/Cauris/Cauris_output/SplitReads/RG"
BQSR="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recal_table"
recal_bam="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated"
for file in `ls ${split_RG}/*.bam`; do
    bam=$(basename "${file}")
    gatk --java-options "-Xmx128G" ApplyBQSR \
        --reference "${REF_DIR}/genome.fa" \
        --input "${split_RG}/${bam}" \
        --bqsr-recal-file "${BQSR}/${bam/splitReads.RG.bam/recal.table}" \
        --output "${recal_bam}/${bam/RG/recal}"
done
## 3. Second pass BQSR on the recalibrated bam files
BQSR2="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recal_table2"
recal_bam="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated"
for file in `ls ${recal_bam}/*.bam`; do
    bam=$(basename "${file}")
    gatk --java-options "-Xmx128G" BaseRecalibrator \
        --reference "${REF_DIR}/genome.fa" \
        --input "${recal_bam}/${bam}" \
        --output "${BQSR2}/${bam/splitReads.recal.bam/recal.table2}" \
        --known-sites "/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11221_CA05_contigs_rename.vcf.gz" \
        --known-sites "/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11245_CA06_contigs_rename.vcf.gz" \
        --known-sites "/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/cauris.consensus_contigs_rename_split_edited.vcf.gz"
done
## 4. Evaluate and compare base quality score recalibration tables
BQSR="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recal_table"
BQSR2="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recal_table2"
covariate="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/covariate"
for table in `ls ${BQSR}`; do
    gatk --java-options "-Xmx128G" AnalyzeCovariates \
        -before "${BQSR}/${table}" \
        -after "${BQSR2}/${table/table/table2}" \
        -plots "${covariate}/${table/recal.table/covariates.pdf}" \
        -csv "${covariate}/${table/recal.table/covariates.csv}"
done

# Germline variant calling with HaplotypeCaller
REF_DIR="/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
recal_bam="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated"
variant_calling="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller"
for file in `ls ${recal_bam}/*.bam`; do
    bam=$(basename "${file}")
    gatk --java-options "-Xmx128G" HaplotypeCaller \
        --reference "${REF_DIR}/genome.fa" \
        --input "${recal_bam}/${bam}" \
        --output "${variant_calling}/${bam/splitReads.recal.bam/haplotypecaller.vcf.gz}" \
        --dont-use-soft-clipped-bases \
        -stand-call-conf 20
done

# Extract only SNP variants
raw_vcf="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/raw_vcf"
snp="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/SNPs"
for vcf in `ls ${raw_vcf} | grep -v "tbi"`; do
    gatk SelectVariants \
        -R "${REF_DIR}/genome.fa" \
        -V "${raw_vcf}/${vcf}" \
        --select-type-to-include SNP \
        -O "${snp}/${vcf/haplotypecaller/haplotypecaller.snps}"
done

# SNPs filtering
## Hard filtering (labeling)
snp="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/SNPs"
for vcf in `ls ${snp} | grep -v "tbi"`; do
    gatk VariantFiltration \
        -R "${REF_DIR}/genome.fa" \
        -V "${snp}/${vcf}" \
        -O "${snp}/${vcf/snps/filtered.snps}" \
        --filter-name "QD_filter" -filter "QD < 2.0" \
        --filter-name "FS_filter" -filter "FS > 60.0" \
        --filter-name "MQ_filter" -filter "MQ < 40.0" \
        --filter-name "SOR_filter" -filter "SOR > 3.0" \
        --filter-name "MQRankSum_filter" -filter "MQRankSum <-12.5" \
        --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum <-8.0"
done
## Select Variants that "PASS" filters
snp="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/SNPs"
filtered_snp="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/SNPs/filtered"
for file in `ls ${snp}/*filtered.snps* | grep -v "tbi"`; do
    labeled_vcf=$(basename "${file}")
    gatk SelectVariants \
        --exclude-filtered \
        -V "${snp}/${labeled_vcf}" \
        -O "${filtered_snp}/${labeled_vcf/filtered/filteredPASS}"
done
## Decomposing(splitting) multi-allelic snp variants and removal of invalid variant sites
filtered_snp="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/SNPs/filtered"
norm_filtered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/SNPs/filtered/norm_filtered"
for file in `ls ${filtered_snp}/*.vcf.gz`; do
    filtered_vcf=$(basename "${file}")
    bcftools norm -m- "${filtered_snp}/${filtered_vcf}" | awk '{if(/^#/) {print $0} else if($5 !~ /*/) {print $0} }' | bgzip > "${norm_filtered}/${filtered_vcf/filteredPASS/norm.filteredPASS}"
done

### Indexing
for file in `ls ${norm_filtered}`; do
    echo $file
    tabix -p vcf $file
done

# Merge all individual vcf files into a single vcf file
bcftools merge -Oz -o merged.vcf.gz *.vcf.gz
bcftools norm -m- -Oz -o merged.norm.vcf.gz merged.vcf.gz   # Splitting multi-allelic sites

### Testing
MAX_AMB_SAMPLES=$(echo "15 10" | awk '{x=$1*($2/100); y=int(x); x=(y<1?1:y)} END {print x}')
/mnt/rdisk/dminh/mycosnp-nf/bin/vcfSnpsToFasta.py --max_amb_samples $MAX_AMB_SAMPLES --min_depth 10 "merged.norm.vcf" > combined_vcf-to-fasta.fasta


### Consensus testing
for file in *.vcf.gz; do 
    sample=${file/.haplotypecaller.norm.filteredPASS.snps.vcf.gz/}
    echo $sample
    echo ">${sample}" > "consensus_test/${sample}.fasta"
    bcftools consensus -s "${sample}" -f "${REF_DIR}/genome.fa" "${sample}.haplotypecaller.norm.filteredPASS.snps.vcf.gz" | grep -E -v "^>" | grep -E -v "^$" >> "consensus_test/${sample}.fasta"
done
cat *.fasta > combined.fasta     ### Concat all consensus fasta files including reference genome
snp-dists combined.fasta        ### SNPs distance matrix

fasttree \
    -gtr -gamma -fastest \
    -log fasttree_phylogeny.tre.log \
    -nt combined.fasta \
    > fasttree_phylogeny.tre