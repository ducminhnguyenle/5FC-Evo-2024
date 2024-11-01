#!/usr/bin/env bash
# Prepare the receptor for AutoDock
## Clean up the raw homology pdb file
## removing atoms from ligands, water molecules and other non-protein/non-nucleic residues
prepare_receptor4.py \
    -r "Model3-FUR1-R214T-monomer-cleaned.pdb" \
    -o "Model3-FUR1-R214T-monomer-prep.pdbqt" \
    -A hydrogens \
    -U nphs_lps \
    -v

# Using Vina forcefield (predefined space search box)
## 5FU ligand
vina \
    --receptor "Model3-FUR1-R214T-monomer-prep.pdbqt" \
    --ligand "5FU.pdbqt" \
    --config "config/FUR1-R214T.config" \
    --out "FUR1-R214T/model3-FUR1-R214T_5FU_vinaFF_out.pdbqt" \
    2>&1 | tee "FUR1-R214T/output_R214T_5FU.log"

## PRPP ligand
vina \
    --receptor "Model3-FUR1-R214T-monomer-prep.pdbqt" \
    --ligand "PRPP.pdbqt" \
    --config "config/FUR1-R214T.config" \
    --out "FUR1-R214T/model3-FUR1-R214T_PRPP_vinaFF_out.pdbqt" \
    2>&1 | tee "FUR1-R214T/output_R214T_PRPP.log"

## 5FU + PRPP
vina \
    --receptor "Model3-FUR1-R214T-monomer-prep.pdbqt" \
    --ligand "PRPP.pdbqt" "5FU.pdbqt" \
    --config "config/FUR1-R214T.config" \
    --out "FUR1-R214T/model3-FUR1-R214T_PRPP-5FU_vinaFF_out.pdbqt" \
    2>&1 | tee "FUR1-R214T/output_R214T_PRPP-5FU.log"
