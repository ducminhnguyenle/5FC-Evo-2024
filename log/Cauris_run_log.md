# _Candida auris_

## 2 groups only (WT vs Treated with 5FC)

```bash
nextflow run main.nf \
    -work-dir "/home/minhnguyen/work/" \
    --design "/mnt/rdisk/minhnguyen/Cauris/input/design.csv" \
    --comparisons "/mnt/rdisk/minhnguyen/Cauris/input/comparisons.csv" \
    --outdir "/home/minhnguyen/Cauris_output/" \
    --genome "Cand_auris_B8441_V2" \
    --protocol "illumina" \
    --skip_trimming true \
    --star_twopass \
    --name "Cauris" \
    -profile docker \
    -resume \
    2>&1 | tee Cauris_headnode.log
```

## Rerun with comparisons for each cell line

```bash
nextflow run main.nf \
    -work-dir "/mnt/rdisk/minhnguyen/Cauris/work/" \
    --design "/mnt/rdisk/minhnguyen/Cauris/input/design.csv" \
    --comparisons "/mnt/rdisk/minhnguyen/Cauris/input/comparisons.csv" \
    --outdir "/mnt/rdisk/minhnguyen/Cauris/output_comparisons/" \
    --genome "Cand_auris_B8441_V2" \
    --merged_counts "/home/minhnguyen/Cauris_output/featureCounts/merged_gene_counts.txt" \
    --name "Cauris" \
    -profile docker \
    2>&1 | tee Cauris_headnode2.log
```
