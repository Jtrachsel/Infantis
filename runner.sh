#!/bin/bash

rm metadata/PDG*

Rscript scripts/01_setup.R
Rscript scripts/02_ppangg_setup.R
Rscript scripts/03_overview_tree_setup.R
./scripts/03_core_genome_overview_tree.sh
Rscript scripts/04_SNP_tree_download.R
Rscript scripts/05_phydelity_setup.R
Rscript scripts/06_pESI_ref_extract.R
# SCRIPT TO RUN MINIMAP2 and BLASTN
Rscript scripts/07_parse_phydelity.R
Rscript scripts/08_plots_n_stats.R
