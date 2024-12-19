#!/usr/bin/env bash

# Calculate the average coverage in WT bam sample
avg_coverage=$(samtools depth -a "WGS/output/samples/S48/finalbam/S48.bam" | awk '{sum+=$3; count++} END {print sum/count}')
# Binning Cauris's genome into 1kb bin size and calculate mean coverage per bin
bedtools makewindows -g "Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta/genome.fa.fai" -w 1000 > "WGS/CNV/genome_1kb.bed"
bedtools coverage -a "WGS/CNV/genome_1kb.bed" -b "WGS/output/samples/S48/finalbam/S48.bam" -mean > "WGS/CNV/S48_window1kb_coverage.bed"
# Remove any bins with empty coverage and log2 scale
awk '($4 != 0){print $0}' "WGS/CNV/S48/S48_window1kb_coverage.bed" | sed '1i scaffold\tstart\tend\tmean_cov' > "WGS/CNV/S48/S48_window1kb_coverage.tsv"
awk -v avg="${avg_coverage}" '{print $1, log($4/avg)/log(2)}' "WGS/CNV/S48/S48_window1kb_coverage.tsv" | less -S

# Create a bait interval
awk 'BEGIN{FS=OFS="\t"}; {gsub("gnl\\|WGS:PEKT\\|mrna_","",$0); print $1,$2,$3,$4}' "Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Annotation/Genes/genes.bed" | sed 's/T0$//'  > "WGS/CNV/ref/baits.bed"

# CNV calling: S48 (S48-YPD used as reference)
cnvkit.py batch \
    "WGS/output/samples/S48/finalbam/S48.bam" \
    --normal "WGS/output/samples/S48-YPD/finalbam/S48-YPD.bam" \
    -m wgs \
    -f "Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta/genome.fa" \
    -t "WGS/CNV/ref/baits.bed" \
    --output-reference "WGS/CNV/ref/S48-YPD_reference.cnn" --output-dir "WGS/CNV/S48/" \
    --scatter

# CNV calling: S187 (S48 and S48-YPD used as normal samples for reference)
cnvkit.py batch \
    "WGS/output/samples/S187/finalbam/S187.bam" \
    --normal "WGS/output/samples/S48/finalbam/S48.bam" "WGS/output/samples/S48-YPD/finalbam/S48-YPD.bam" \
    -m wgs \
    -f "refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta/genome.fa" \
    -t "WGS/CNV/ref/baits.bed" \
    --output-reference "WGS/CNV/ref/reference.cnn" --output-dir "WGS/CNV/S187/"

# CNV calling: S191 (Reusing reference.cnn from S48 and S48-YPD)
cnvkit.py batch \
    "WGS/output/samples/S191/finalbam/S191.bam" \
    -r "WGS/CNV/ref/reference.cnn" \
    -d "WGS/CNV/S191/"
