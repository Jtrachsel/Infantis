library(tidyverse)
library(pdtools)

infmeta <- read_tsv('metadata/Infantis_metadata.tsv')



set.seed(7)
tree_dat <- 
  infmeta %>%
  group_by(PDS_acc) %>%  
  summarise(num_genomes=n(), 
            rep_genome=sample(asm_acc, size = 1 )) %>% 
  arrange(desc(num_genomes)) %>% 
  mutate(fasta_path=map_chr(rep_genome, ~list.files(path='assemblies', pattern=paste0(.x, '.fna'), full.names = T) )) %>% 
  write_tsv('output/tree_rep_genomes.tsv')




build_ppanggolin_file_fastas(incomplete_genome_paths = tree_dat$fasta_path) %>% 
  write_tsv(col_names = FALSE, 'overview_tree_ppanggolin.tsv')


### this doesnt work
# 
# tree_dat <- tree_dat %>% 
#   mutate(fasta_path=map_chr(rep_genome, ~list.files(path='assemblies', pattern=paste0(.x, '.fna'), full.names = T) )) 
# 

infmeta %>% group_by(SERO) %>% tally() %>% arrange(desc(n))
infmeta %>% group_by(country) %>% tally() %>% arrange(desc(n))
infmeta %>% group_by(ag_match, country) %>% tally() %>% arrange(desc(n))

LOOK <- infmeta %>% filter(grepl('_', ag_match))
