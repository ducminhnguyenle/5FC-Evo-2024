# Merges multiple bam files with the same Cauris strain into one a single bam file
bam="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated"
merged_bam="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated/merged_bam"
input=("${bam}/F1_185_1.splitReads.recal.bam" "${bam}/F4_188_1.splitReads.recal.bam" "${bam}/F7_189_1.splitReads.recal.bam" "${bam}/F10_191_1.splitReads.recal.bam" "${bam}/S4_48_1.splitReads.recal.bam")
for file in "${input[@]}"; do
    file=$(basename $file)
    ## Strain 185
    if [[ ${file} =~ "185" ]]; then
        gatk --java-options "-Xmx128G" MergeSamFiles \
        -I "${bam}/F1_185_1.splitReads.recal.bam" -I "${bam}/F2_185_2.splitReads.recal.bam" -I "${bam}/F3_185_3.splitReads.recal.bam" \
        -O "${merged_bam}/strain_185_merged.bam"
    ## Strain 188
    elif [[ ${file} =~ "188" ]]; then
        gatk --java-options "-Xmx128G" MergeSamFiles \
        -I "${bam}/F4_188_1.splitReads.recal.bam" -I "${bam}/F5_188_2.splitReads.recal.bam" -I "${bam}/F6_188_3.splitReads.recal.bam" \
        -O "${merged_bam}/strain_188_merged.bam"
    ## Strain 189
    elif [[ ${file} =~ "189" ]]; then
        gatk --java-options "-Xmx128G" MergeSamFiles \
        -I "${bam}/F7_189_1.splitReads.recal.bam" -I "${bam}/F8_189_2.splitReads.recal.bam" -I "${bam}/F9_189_3.splitReads.recal.bam" \
        -O "${merged_bam}/strain_189_merged.bam"
    ## Strain 191
    elif [[ ${file} =~ "191" ]]; then
        gatk --java-options "-Xmx128G" MergeSamFiles \
        -I "${bam}/F10_191_1.splitReads.recal.bam" -I "${bam}/F11_191_2.splitReads.recal.bam" -I "${bam}/F12_191_3.splitReads.recal.bam" \
        -O "${merged_bam}/strain_191_merged.bam"
    ## Strain 48 (WT)
    elif [[ ${file} =~ "48" ]]; then
        gatk --java-options "-Xmx128G" MergeSamFiles \
        -I "${bam}/S4_48_1.splitReads.recal.bam" -I "${bam}/S5_48_2.splitReads.recal.bam" -I "${bam}/S6_48_3.splitReads.recal.bam" \
        -O "${merged_bam}/strain_48_merged.bam"
    fi
done

# Germline variant calling with HaplotypeCaller in gVCF (Phred score confident = 30)
REF_DIR="/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
MERGED_BAM="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated/merged_bam"
MERGED_BAM_GVCF="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/gvcf/merged_bam_gvcf"
for file in `ls ${MERGED_BAM}/*.bam`; do
    BAM=$(basename "${file}")
    gatk --java-options "-Xmx128G" HaplotypeCaller \
        -R "${REF_DIR}/genome.fa" \
        -I "${MERGED_BAM}/${BAM}" \
        -O "${MERGED_BAM_GVCF}/${BAM/_merged.bam/.haplotypecaller.g.vcf.gz}" \
        --dont-use-soft-clipped-bases \
        -ERC GVCF --sample-ploidy "1" \
        -stand-call-conf 30
done

# Check single-sample or multi-sample bam file
samtools samples "/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated/merged_bam/strain_185_merged.bam" |\
cut -f1 | sort | uniq | cat -n      # multi-sample bam

----------------------------------------------------------------------------------------------------------------------------