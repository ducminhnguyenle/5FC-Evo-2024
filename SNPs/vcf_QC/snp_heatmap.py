#!/usr/bin/env python3
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from natsort import natsort_keygen

# Func to extract variant depth based on genotype
def get_variant_depth(cell):
    elements = cell.split(":")
    
    if elements[0] == "0":
        return 0
    elif elements[0] == ".":
        return pd.NA
    elif elements[0] == "1":
        if len(elements) >= 3:
            try:
                return int(elements[2]) - int(elements[1].split(",")[0])
            except ValueError:
                return pd.NA
        else:
            return pd.NA
    else:
        return(":".join(elements[:3]))

def genotype_label(cell):
    if pd.isna(cell):
        return "."
    elif cell > 0:
        return 1
    elif cell == 0:
        return 0
    else:
        return cell

# Prepare snps of interest data
snp_df = pd.read_csv("cauris.snps.of.interest.tsv", sep = "\t", header = 0)

# Rename scaffold 
scaffold = snp_df["CHROM"].str.split(".", n = 1, expand = True)
snp_df["CHROM"] = scaffold[0]

# Apply func to each cell
snp_df.iloc[:, 6:] = snp_df.iloc[:, 6:].applymap(get_variant_depth)

# Concat "CHROM" and "POS" columns together and drop unnecssary columns
snp_df["CHR_POS"] = snp_df["CHROM"].str.cat(snp_df["POS"].astype(str), sep=":")
snp_strain = snp_df.drop(labels = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"], axis = 1)

# Set "CHR_POS" column to be the new index
snp_strain.set_index("CHR_POS", inplace=True)

# Covert all values to numeric for compatible with heatmap plot
snp_of_strain_df = snp_strain.apply(pd.to_numeric)
columns = {
    "strain_185": "185-R",
    "strain_188": "188-S",
    "strain_189": "189-R",
    "strain_191": "191-I",
    "strain_48": "WT-S"
}
# Rename sample inlucding phenotype and rerrange based on matching phenotype
snp_of_strain_df.rename(columns = columns, inplace = True)
snp_of_strain_df = snp_of_strain_df.reindex(columns = ["185-R", "189-R", "191-I", "188-S", "WT-S"])
# Reorder scaffold
snp_of_strain_df.sort_values(
    by = "CHR_POS",
    key = natsort_keygen(),
    inplace = True
)
# print(snp_of_strain_df)

# Labelling appropriate genotype
ann_snp_of_strain_df = snp_of_strain_df.applymap(genotype_label)
# print(ann_snp_of_strain_df)

# Heatmap plot
plt.figure(figsize = (10, 8))
my_cmap = sns.color_palette("flare", as_cmap = True).copy()
my_cmap.set_under("#468B97")
ax = sns.heatmap(snp_of_strain_df, annot = ann_snp_of_strain_df, fmt = "", cmap = my_cmap, square = True,
                 linewidth = 0.5, linecolor = "white", mask = snp_of_strain_df.isnull(),
                 cbar_kws = {"shrink": 0.6, "extend": "min", "extendrect": True},
                 vmin = 0.1
                 )
ax.set(xlabel = "", ylabel = "")
ax.xaxis.tick_top()
plt.xticks(weight = "bold")
plt.tick_params(
    axis = "x",
    which = "both",
    top = False
)
plt.tick_params(
    axis = "y",
    which = "both",
    left = False
)
plt.show()