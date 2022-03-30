library(tidyverse)

infmeta <- read_tsv('metadata/Infantis_metadata.tsv')


infmeta %>% count(Year, PDS_acc) %>% arrange(desc(n))

infmeta %>% count(country, PDS_acc) %>% arrange(desc(n))

# UK strains?
# one isolate from pork in denmark 2018, the rest are from 'food' or 'environmental' isolations
# in the UK
look <- 
  infmeta %>%
  filter(PDS_acc == 'PDS000073471.105') %>% 
  filter(ag_match != 'Human') #%>% select(Year,country,host, ag_match)


### This one!
infmeta %>% 
  mutate(PDS_lump=fct_lump_n(f = PDS_acc, n = 8)) %>% 
  count(Year, PDS_lump) %>% 
  filter(Year < 2022) %>% 
  filter(Year > 2015) %>% 
  ggplot(aes(x=Year, y=n, color=PDS_lump)) +
  geom_line(size=1.5, color='black', aes(group=PDS_lump))+
  geom_line(size=1.2) +
  ylab('number of genomes') + 
  guides(color=guide_legend(title="SNP cluster")) +
  scale_color_brewer(palette = 'Set1') + 
  theme_bw()

set.seed(7)
tree_dat <- 
  infmeta %>%
  group_by(PDS_acc) %>%  
  summarise(num_genomes=n(), 
            rep_genome=sample(asm_acc, size = 1 )) %>% 
  arrange(desc(num_genomes)) %>% 
  mutate(fasta_path=map_chr(rep_genome, ~list.files(path='assemblies', pattern=paste0(.x, '.fna'), full.names = T) )) %>% 
  write_tsv('output/tree_rep_genomes.tsv')

library(pdtools)



build_ppanggolin_file_fastas(incomplete_genome_paths = tree_dat$fasta_path) %>% 
  write_tsv(col_names = FALSE, 'overview_tree_ppanggolin.tsv')


### this doesnt work

tree_dat <- tree_dat %>% 
  mutate(fasta_path=map_chr(rep_genome, ~list.files(path='assemblies', pattern=paste0(.x, '.fna'), full.names = T) )) 


infmeta %>% group_by(SERO) %>% tally() %>% arrange(desc(n))
infmeta %>% group_by(country) %>% tally() %>% arrange(desc(n))
infmeta %>% group_by(ag_match, country) %>% tally() %>% arrange(desc(n))

LOOK <- infmeta %>% filter(grepl('_', ag_match))
