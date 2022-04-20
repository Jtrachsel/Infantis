#!/bin/bash

ppanggolin annotate --fasta overview_tree_ppanggolin.tsv --cpu 60 -o overview_tree -f
ppanggolin cluster -p overview_tree/pangenome.h5 --cpu 60
ppanggolin graph -p overview_tree/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p overview_tree/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p overview_tree/pangenome.h5 -o overview_tree/MSA --cpu 60
#ppanggolin rgp -p overview_tree/pangenome.h5 -c 16
#ppanggolin module -p overview_tree/pangenome.h5 -c 16
#ppanggolin spot -p overview_tree/pangenome.h5 -c 16
#ppanggolin write -o overview_tree/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p overview_tree/pangenome.h5
#ppanggolin fasta -p ./all_pan/pangenome.h5 --output REP_PROTS --prot_families all






raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n overview_tree -s ./overview_tree/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7

