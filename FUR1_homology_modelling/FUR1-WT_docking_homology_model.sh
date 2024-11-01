#!/usr/bin/env bash
# Prepare the receptor for AutoDock
## Clean up the raw homology pdb file
## removing atoms from ligands, water molecules and other non-protein/non-nucleic residues
for receptor in FUR1_homology_modelling/*.pdb; do
    prepare_receptor4.py \
        -r "${receptor}" \
        -o "${receptor/cleaned.pdb/prep.pdbqt}" \
        -A hydrogens \
        -U nphs_lps \
        -v
done

# Prepare 5FU and PRPP ligands
wget https://go.drugbank.com/structures/small_molecule_drugs/DB00544.smiles
wget https://go.drugbank.com/structures/small_molecule_drugs/DB01632.smiles
## Convert SMILES to SDF format
obabel \
    -i smi "DB00544.smiles" \
    -o sdf \
    -O "5FU.sdf" \
    -p 7.4 \
    --gen3d --best --canonical
obabel \
    -i smi "DB01632.smiles" \
    -o sdf \
    -O "PRPP.sdf" \
    -p 7.4 \
    --gen3d --best --canonical
## Convert ligand into AutoDock structure
### 5FU
mk_prepare_ligand.py \
    -i "5FU.sdf" \
    -o "5FU.pdbqt"
### PRPP
mk_prepare_ligand.py \
    -i "PRPP.sdf" \
    -o "PRPP.pdbqt"

# Using Vina forcefield (predefined space search box)
## 5FU ligand and PRPP ligand separately
ligands=("5FU.pdbqt" "PRPP.pdbqt")
for ligand in "${ligands[@]}"; do
    ligand_name=${ligand%.*}
    vina \
        --receptor "Model3-FUR1-WT-monomer-prep.pdbqt" \
        --ligand "${ligand}" \
        --config "config/FUR1-WT.config" \
        --out "Model3-FUR1-WT_${ligand_name}_vinaFF_out.pdbqt" \
        2>&1 | tee "Model3-FUR1-WT_${ligand_name}_vinaFF_out.log"
done
## 5FU + PRPP
vina \
    --receptor "Model3-FUR1-WT-monomer-prep.pdbqt" \
    --ligand "PRPP.pdbqt" "5FU.pdbqt" \
    --config "config/FUR1-WT.config" \
    --out "FUR1-WT/model3-FUR1-WT_PRPP-5FU_vinaFF_out.pdbqt" \
    2>&1 | tee "FUR1-WT/output_WT_PRPP-5FU.log"
