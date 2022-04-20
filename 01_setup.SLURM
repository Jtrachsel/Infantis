#!/bin/bash

#SBATCH --job-name=setup                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 4                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=35G   
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

# ENTER COMMANDS HERE:
module load miniconda
module load R
source activate /project/fsepru/conda_envs/ppanggolin/

ppanggolin annotate --fasta all_ppangg.tsv --cpu 70 -o all_pan -f
ppanggolin cluster -p all_pan/pangenome.h5 --cpu 70
ppanggolin graph -p all_pan/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p all_pan/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p all_pan/pangenome.h5 -o all_pan/MSA --cpu 70
ppanggolin rgp -p all_pan/pangenome.h5 -c 16
ppanggolin module -p all_pan/pangenome.h5 -c 16
ppanggolin spot -p all_pan/pangenome.h5 -c 16
ppanggolin write -o all_pan/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p all_pan/pangenome.h5
ppanggolin fasta -p ./all_pan/pangenome.h5 --output REP_PROTS --prot_families all




