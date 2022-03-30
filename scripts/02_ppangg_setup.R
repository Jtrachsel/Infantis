library(pdtools)
library(tidyverse)

## add check for fna existence here ##


Infantis_metadata <- read_tsv('./metadata/Infantis_metadata.tsv') %>% 
  filter(fna_exists) %>% 
  mutate(fna_path=sub('.gz', '',fna_dest))


complete_genomes <- 
  Infantis_metadata %>%
  filter(asm_level == 'Complete Genome') 

draft_genomes <- 
  Infantis_metadata %>%
  filter(asm_level != 'Complete Genome') 



ppangg <- 
  build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes$fna_path, 
                               incomplete_genome_paths = draft_genomes$fna_path) %>% 
  write_tsv('all_ppangg.tsv', col_names = FALSE)





Infantis_metadata %>% 
  # check_files_exist() %>% 
  filter(fna_exists) %>% 
  count(ag_match) %>%
  arrange(desc(n))

check <- Infantis_metadata %>%
  filter(ag_match == 'Bovine_Goat') %>% 
  pull(isolation_source)

Infantis_metadata %>%
  count(SERO) %>%
  arrange(desc(n))

Infantis_metadata %>% 
  count(SERO, PDS_acc) %>% 
  # arrange(desc(n)) %>% 
  group_by(PDS_acc) %>% 
  summarise(num_infantis=n[which(SERO == 'Infantis')], 
            num_others)

