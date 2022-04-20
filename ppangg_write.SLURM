#!/bin/bash

#SBATCH --job-name=ppangg                              # name of the job submitted
#SBATCH -p scavenger                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 4                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 8:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=500G   
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

# ENTER COMMANDS HERE:
module load miniconda
source activate /project/fsepru/conda_envs/ppanggolin_dev

ppanggolin write -o all_pan/WRITE -f --projection --regions --modules --families_tsv -p ./all_pan/pangenome.h5

#End of file
