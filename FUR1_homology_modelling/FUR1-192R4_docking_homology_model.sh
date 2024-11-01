#!/usr/bin/env bash
# Prepare the receptor for AutoDock
## Clean up the raw homology pdb file
## removing atoms from ligands, water molecules and other non-protein/non-nucleic residues
prepare_receptor4.py \
    -r "Model3-FUR1-192R4-monomer-cleaned.pdb" \
    -o "Model3-FUR1-192R4-monomer-prep.pdbqt" \
    -A hydrogens \
    -U nphs_lps \
    -v

# Using Vina forcefield (predefined space search box)
## 5FU ligand
vina \
    --receptor "Model3-FUR1-192R4-monomer-prep.pdbqt" \
    --ligand "5FU.pdbqt" \
    --config "config/FUR1-192R4.config" \
    --out "FUR1-192R4/model3-FUR1-192R4_5FU_vinaFF_out.pdbqt" \
    2>&1 | tee "FUR1-192R4/output_192R4_5FU.log"

## PRPP ligand
vina \
    --receptor "Model3-FUR1-192R4-monomer-prep.pdbqt" \
    --ligand "PRPP.pdbqt" \
    --config "config/FUR1-192R4.config" \
    --out "FUR1-192R4/model3-FUR1-192R4_PRPP_vinaFF_out.pdbqt" \
    2>&1 | tee "FUR1-192R4/output_192R4_PRPP.log"

## 5FU + PRPP
vina \
    --receptor "Model3-FUR1-192R4-monomer-prep.pdbqt" \
    --ligand "PRPP.pdbqt" "5FU.pdbqt" \
    --config "config/FUR1-192R4.config" \
    --out "FUR1-192R4/model3-FUR1-192R4_PRPP-5FU_vinaFF_out.pdbqt" \
    2>&1 | tee "FUR1-192R4/output_192R4_PRPP-5FU.log"
