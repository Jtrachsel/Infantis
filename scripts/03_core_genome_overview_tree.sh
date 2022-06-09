#!/bin/bash

# For overview tree showing relatedness of SNP clusters

conda activate ppanggolin

ppanggolin annotate --fasta overview_tree_ppanggolin.tsv --cpu 50 -o overview_tree -f
ppanggolin cluster -p overview_tree/pangenome.h5 --cpu 50
ppanggolin graph -p overview_tree/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p overview_tree/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p overview_tree/pangenome.h5 -o overview_tree/MSA --cpu 50
# ppanggolin rgp -p overview_tree/pangenome.h5 -c 16
# ppanggolin module -p overview_tree/pangenome.h5 -c 16
# ppanggolin spot -p overview_tree/pangenome.h5 -c 16
# ppanggolin write -o overview_tree/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p overview_tree/pangenome.h5
# ppanggolin fasta -p ./all_pan/pangenome.h5 --output REP_PROTS --prot_families all





trimal -in core_genome_alignment.aln -out TRIMAL_core.aln -fasta -noallgaps 

raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n overview_tree -s ./overview_tree/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7
raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n overview_tree_reduced -s ./overview_tree/MSA/core_genome_alignment.aln.reduced -T 35 -x 7 -N autoMRE -p 7
raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n overview_TRIMAL_core -s ./overview_tree/MSA/TRIMAL_core.aln -T 30 -x 7 -N autoMRE -p 7



