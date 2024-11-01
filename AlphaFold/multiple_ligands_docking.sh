#!/usr/bin/env bash
# Multiple ligands docking

## FUR1 WT with PRPP and 5FU ligands
vina \
    --receptor "AF-A0A2H0ZP91-F1-model_v4.pdbqt" \
    --ligand "PRPP.pdbqt" "5FU.pdbqt" \
    --config "autodock_config/config_singledock" \
    --out "FUR1-WT_PRPP-5FU_vinaFF_out.pdbqt" \
    2>&1 | tee output_WT_PRPP_5FU.log

## FUR1 R214T mutation with PRPP and 5FU ligands
vina \
    --receptor "FUR1_R214T.pdbqt" \
    --ligand "PRPP.pdbqt" "5FU.pdbqt" \
    --config "autodock_config/config_singledock" \
    --out "FUR1-R214T_PRPP-5FU_vinaFF_out.pdbqt" \
    2>&1 | tee output_R214T_PRPP_5FU.log