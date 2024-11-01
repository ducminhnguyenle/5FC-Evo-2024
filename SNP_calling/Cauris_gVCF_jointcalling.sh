# Germline variant calling with HaplotypeCaller in gVCF (Phred score confident = 30)
REF_DIR="/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
recal_bam="/mnt/rdisk/dminh/Cauris/Cauris_output/BaseRecalibrator/recalibrated"
gvcf="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/gvcf"
for file in `ls ${recal_bam}/*.bam`; do
    bam=$(basename "${file}")
    gatk --java-options "-Xmx256G" HaplotypeCaller \
        -R "${REF_DIR}/genome.fa" \
        -I "${recal_bam}/${bam}" \
        -O "${gvcf}/${bam/splitReads.recal.bam/haplotypecaller.g.vcf.gz}" \
        --dont-use-soft-clipped-bases \
        -ERC GVCF --sample-ploidy "1"
done

# GATK4_LOCALCOMBINEGVCFS: Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
REF_DIR="/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
combined="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/combinegvcfs"
gvcf="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/gvcf"
gatk --java-options "-Xmx128g" CombineGVCFs \
      -R "${REF_DIR}/genome.fa" \
      -O "${combined}/cauris.combined.g.vcf.gz" \
      -V "${gvcf}/F10_191_1.haplotypecaller.g.vcf.gz" -V "${gvcf}/F11_191_2.haplotypecaller.g.vcf.gz" -V "${gvcf}/F1_185_1.haplotypecaller.g.vcf.gz" -V "${gvcf}/F12_191_3.haplotypecaller.g.vcf.gz" -V "${gvcf}/F2_185_2.haplotypecaller.g.vcf.gz" -V "${gvcf}/F3_185_3.haplotypecaller.g.vcf.gz" -V "${gvcf}/F4_188_1.haplotypecaller.g.vcf.gz" -V "${gvcf}/F5_188_2.haplotypecaller.g.vcf.gz" -V "${gvcf}/F6_188_3.haplotypecaller.g.vcf.gz" -V "${gvcf}/F7_189_1.haplotypecaller.g.vcf.gz" -V "${gvcf}/F8_189_2.haplotypecaller.g.vcf.gz" -V "${gvcf}/F9_189_3.haplotypecaller.g.vcf.gz" -V "${gvcf}/S4_48_1.haplotypecaller.g.vcf.gz" -V "${gvcf}/S5_48_2.haplotypecaller.g.vcf.gz" -V "${gvcf}/S6_48_3.haplotypecaller.g.vcf.gz"

# GATK4_GENOTYPEGVCFS: Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
combined="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/combinegvcfs"
genotype="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/genotypegvcfs"
gatk --java-options "-Xmx128g" GenotypeGVCFs \
    -R "${REF_DIR}/genome.fa" \
    -V "${combined}/cauris.combined.g.vcf.gz" \
    -O "${genotype}/cauris_combined_genotype.vcf.gz"

# Extract only SNP variants
genotype="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/genotypegvcfs"
snp="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/selectedsnps"
gatk --java-options "-Xmx64G" SelectVariants \
    -R "${REF_DIR}/genome.fa" \
    -V "${genotype}/cauris_combined_genotype.vcf.gz" \
    -O "${snp}/cauris_combined_genotype_snps_selectvariants.vcf.gz" \
    --select-type-to-include "SNP"

# SNPs hard filtering
snp="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/selectedsnps"
snpfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/selectedsnpsfiltered"
gatk --java-options "-Xmx64G" VariantFiltration \
        -R "${REF_DIR}/genome.fa" \
        -V "${snp}/cauris_combined_genotype_snps_selectvariants.vcf.gz" \
        -O "${snpfiltered}/cauris_combined_genotype_snps_filtered.vcf.gz" \
        --filter-name "QD_filter" -filter "QD < 2.0" \
        --filter-name "FS_filter" -filter "FS > 60.0" \
        --filter-name "MQ_filter" -filter "MQ < 40.0" \
        --filter-name "SOR_filter" -filter "SOR > 3.0" \
        --filter-name "MQRankSum_filter" -filter "MQRankSum <-12.5" \
        --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum <-8.0"

# Select SNP variants that "PASS" filters (based on INFO)
snpfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/selectedsnpsfiltered"
gatk SelectVariants \
    --exclude-filtered \
    -V "${snpfiltered}/cauris_combined_genotype_snps_filtered.vcf.gz" \
    -O "${snpfiltered}/cauris_combined_genotype_snps_filteredPASS.vcf.gz"

# FILTER_GATK_GENOTYPES (based on FORMAT) (by masking the GT)
snpfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/selectedsnpsfiltered"
genofiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/genotypesfiltered"
filterGatkGenotypes="/mnt/rdisk/dminh/Cauris/bin/filterGatkGenotypes.py"

gzip -c -d "${snpfiltered}/cauris_combined_genotype_snps_filteredPASS.vcf.gz" > "${snpfiltered}/cauris_combined_genotype_snps_filteredPASS.vcf"

python3 ${filterGatkGenotypes}  "${snpfiltered}/cauris_combined_genotype_snps_filteredPASS.vcf" \
                                --min_GQ "50" \
                                --keep_GQ_0_refs \
                                --min_percent_alt_in_AD "0.8" \
                                --min_total_DP "10" \
                                --keep_all_ref \
                                > "${genofiltered}/cauris_combined_genotype_filtered_snps_filtered.vcf"

bgzip "${genofiltered}/cauris_combined_genotype_filtered_snps_filtered.vcf"
tabix -p vcf "${genofiltered}/cauris_combined_genotype_filtered_snps_filtered.vcf.gz"

# BCFTOOLS_VIEW_CONVERT
finalfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/finalfiltered"
bcftools view \
    -o finalfiltered/finalfiltered.vcf.gz \
    -Oz \
    genotypesfiltered/cauris_combined_genotype_filtered_snps_filtered.vcf.gz

# BCFTOOLS_INDEX
finalfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/finalfiltered"
bcftools index \
    --threads 2 \
    finalfiltered/finalfiltered.vcf.gz
tabix -p vcf finalfiltered/finalfiltered.vcf.gz

# DECOMPOSING and NORMALIZING FINAL VCF FILES
finalfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/finalfiltered"
bcftools norm -m- "${finalfiltered}/finalfiltered.vcf.gz" | awk '{if(/^#/) {print $0} else if($5 !~ /*/) {print $0} }' | bgzip > "${finalfiltered}/finalfiltered_norm.vcf.gz"
tabix -p vcf "${finalfiltered}/finalfiltered_norm.vcf.gz"

# SPLIT_VCF
finalfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/finalfiltered"
splitvcf="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/splitvcf"
bcftools query \
    -l \
    "${finalfiltered}/finalfiltered_norm.vcf.gz" > "${splitvcf}/samplelist.txt"
for SAMPLE in $(cat "${splitvcf}/samplelist.txt"); do
    echo ${SAMPLE}
    bcftools view -Oz -s $SAMPLE -o "${splitvcf}/${SAMPLE}.vcf.gz" "${finalfiltered}/finalfiltered_norm.vcf.gz"
done

# VCF_CONSENSUS
splitvcf="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/splitvcf"
finalfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/finalfiltered"
REF_DIR="/mnt/rdisk/dminh/Cauris/refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
consensus="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/consensus"
for SAMPLE in $(cat "${splitvcf}/samplelist.txt"); do
    echo ">${SAMPLE}" > "${consensus}/${SAMPLE}.fasta"
    bcftools consensus -s ${SAMPLE} -f "${REF_DIR}/genome.fa" ${finalfiltered}/finalfiltered_norm.vcf.gz |  grep -E -v "^>" | grep -E -v "^$" >> "${consensus}/${SAMPLE}.fasta"
    gzip "${consensus}/${SAMPLE}.fasta"
done

# VCF_TO_FASTA (Multi-fasta)
splitvcf="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/splitvcf"
MAX_PERC_AMB_SAMPLES=10
NUM_SAMPLES=$(cat "${splitvcf}/samplelist.txt" | wc -l)
if [ ${MAX_PERC_AMB_SAMPLES} > 0 ]; then
    MAX_AMB_SAMPLES=$(echo "${NUM_SAMPLES} 10" | awk '{x=$1*($2/100); y=int(x); x=(y<1?1:y)} END {print x}')
else
    MAX_AMB_SAMPLES=10000000
fi

finalfiltered="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/finalfiltered"
vcfSnpsToFasta="/mnt/rdisk/dminh/Cauris/bin/vcfSnpsToFasta.py"
vcf_to_fasta="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/vcf_to_fasta"
gzip -c -d "${finalfiltered}/finalfiltered_norm.vcf.gz" > "${finalfiltered}/finalfiltered_norm.vcf"
python3 ${vcfSnpsToFasta} \
    --max_amb_samples $MAX_AMB_SAMPLES \
    --min_depth 10 \
    "${finalfiltered}/finalfiltered_norm.vcf" > "${vcf_to_fasta}/combined_vcf-to-fasta.fasta"

# VCF_QC
vcf_qc="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/vcf_qc_report"
printf "Sample Name\tLength\tNumber-N\n" > "${vcf_qc}/vcf-qc-report.txt"
awk '$0 ~ ">" {if (NR > 1) {print c "\t" d;} c=0;d=0;printf substr($0,2,200) "\t"; } $0 !~ ">" {c+=length($0);d+=gsub(/N/, "");d+=gsub(/n/, "")} END { print c "\t" d; }' "${vcf_to_fasta}/combined_vcf-to-fasta.fasta" >> "${vcf_qc}/vcf-qc-report.txt"

# SEQKIT_REPLACE
vcf_to_fasta="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/vcf_to_fasta"
seqkit \
    replace \
    -s -p '\*' -r '-' \
    --threads 2 \
    -i combined_vcf-to-fasta.fasta \
    -o vcf-to-fasta.fasta

# SNPDISTS
snpdists="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/snpdists"
vcf_to_fasta="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/vcf_to_fasta"
snp-dists \
    -b \
    "${vcf_to_fasta}/vcf-to-fasta.fasta" > "${snpdists}/combined.tsv"

#----------------------------------------CREATE_PHYLOGENY------------------------------------#
## RAPIDNJ: neighbour-joining
rapidnj="/mnt/rdisk/dminh/Cauris/Cauris_output/phylogeny/rapidnj"
python \
    -c 'from Bio import SeqIO; SeqIO.convert("/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/vcf_to_fasta/vcf-to-fasta.fasta", "fasta", "alignment.sth", "stockholm")'

rapidnj \
    alignment.sth \
    -t d -b 1000 -n \
    -i sth \
    -c 2 \
    -x rapidnj_phylogeny.nwk

## FASTTREE: approximately-maximum-likelihood phylogenetic trees
vcf_to_fasta="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/vcf_to_fasta"
fasttree="/mnt/rdisk/dminh/Cauris/Cauris_output/phylogeny/fasttree"
fasttree \
    -gtr -gamma -fastest \
    -log fasttree_phylogeny.tre.log \
    -nt "${vcf_to_fasta}/vcf-to-fasta.fasta" \
    > fasttree_phylogeny.nh

## IQTREE: maximum likelihood
vcf_to_fasta="/mnt/rdisk/dminh/Cauris/Cauris_output/HaplotypeCaller/vcf_to_fasta"
iqtree="/mnt/rdisk/dminh/Cauris/Cauris_output/phylogeny/iqtree"
iqtree \
    -alrt 1000 -B 1000 -m MFP -czb \
    -s "${vcf_to_fasta}/vcf-to-fasta.fasta" \
    -nt AUTO \
    -ntmax 2 \
    -mem 6GB

## RAXMLNG: maximum likelihood
raxml-ng \
    --all --model GTR+G --bs-trees 1000 \
    --msa "${vcf_to_fasta}/vcf-to-fasta.fasta" \
    --threads 2 \
    --prefix output

## QUICKSNP: builds a NJ tree from a SNP distance matrix (Not good)
QuickSNP.py \
    --dm combined.tsv \
    --outtree quicksnp_tree.nwk