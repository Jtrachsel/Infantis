#!/bin/bash

source activate /home/Julian.Trachsel/miniconda3/envs/gifrop2/


mkdir -p plasmid_extraction/REPLACE

# move in draft genomes
while read line
do
	cp $line plasmid_extraction/REPLACE/

done < REPLACE

# move in reference genomes
cp reference_genomes/LT2.fna plasmid_extraction/REPLACE/
cp reference_genomes/FSIS1502916.fna plasmid_extraction/REPLACE/


cd plasmid_extraction/REPLACE/

pan_pipe --roary_args '-p 70 -s' --gifrop_args '--threads 40'
#roary -p 40 -s *gff
#gifrop -r LT2 --threads 40 --get_islands -m 10

cd -

