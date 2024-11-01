#!/usr/bin/env bash
# Molecular docking on AlphaFold structures

# Using AutoDock Vina

## Convert SMILES to SDF format
obabel \
    -i smi 5FU.smiles \
    -o sdf \
    -O 5FU.sdf \
    -p 7.4 \
    --gen3d --best --canonical
## Convert protein into AutoDock structure:
## PDB > prepare_receptor > protein PDBQT
### FUR1 WT
prepare_receptor4.py \
    -r "AF-A0A2H0ZP91-F1-model_v4.pdb" \
    -o "AF-A0A2H0ZP91-F1-model_v4.pdbqt"
### FUR1 R214T mutation
prepare_receptor4.py \
    -r "FUR1_R214T.pdb" \
    -o "FUR1_R214T.pdbqt"

## Convert ligand into AutoDock structure:
## SDF > mk_prepare_ligand.py > ligand PDBQT
### 5FU ligand
mk_prepare_ligand.py \
    -i "5FU.sdf" \
    -o "5FU.pdbqt"
### PRPP ligand
mk_prepare_ligand.py \
    -i "PRPP.sdf" \
    -o "PRPP.pdbqt" 

## Specify AutoGrid4 calculation for defining grid space:
## Protein PDBQT, ligand PDBQT > prepare_gdf4 > ligand GPF
# ./autodock_scripts/prepare_gpf.py \
#     -l "5FU.pdbqt" \
#     -r "AF-A0A2H0ZP91-F1-model_v4.pdbqt" \
#     -y

## Create affinity map using AutoGrid4:
## Ligand GPF > autogrid4 > affinity map XYZ
# /home/dminh/mgltools/autodocksuite_4.2.6/autogrid4 \
#     -p "AF-A0A2H0ZP91-F1-model_v4.gpf" \
#     -l "AF-A0A2H0ZP91-F1-model_v4.glg"

## Perform the docking using AutoDock Vina:
## Protein PDBQT, ligand PDBQT, affinity map XYZ > Vina().dock
### Using AutoDock4 forcefield (not good)
# vina \
#     --ligand "5FU.pdbqt" \
#     --maps "AF-A0A2H0ZP91-F1-model_v4" \
#     --scoring ad4 \
#     --exhaustiveness 32 \
#     --out "FUR1_5FU_ad4_out.pdbqt"
### Using Vina forcefield (predefined space search box) - (much better)
#### FUR1 WT
vina \
    --receptor "AF-A0A2H0ZP91-F1-model_v4.pdbqt" \
    --ligand "PRPP.pdbqt" \
    --config "autodock_config/config_singledock" \
    --out "FUR1-WT_PRPP_vinaFF_out.pdbqt" \
    2>&1 | tee output_WT_PRPP.log
#### FUR1 R214T mutation
vina \
    --receptor "FUR1_R214T.pdbqt" \
    --ligand "PRPP.pdbqt" \
    --config "autodock_config/config_singledock" \
    --out "FUR1-R214T_PRPP_vinaFF_out.pdbqt" \
    2>&1 | tee output_R214T_PRPP.log