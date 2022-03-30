library(tidyverse)
library(Biostrings)
library(pafr)
library(tictoc)
library(furrr)
library(cowplot)

plan(multisession, workers = 16)  
# pafr::plot_coverage()

# extract only the pESI plasmid sequence from a well known Infantis assembly
FSIS <- readDNAStringSet('reference_genomes/FSIS1502916.fna')
FSIS['FSIS1502916_2'] %>% writeXStringSet('./reference_genomes/pESI.fna')

## ran these commands to map all contigs to the pESI plasmid sequence:
# mkdir pESI_screen
#parallel -j 35 'minimap2 -x asm5 -t 1 -o ./pESI_screen/{/.}.PAF ./reference_genomes/pESI.fna {}' ::: ./assemblies/*fna


return_percent_pESI <- function(path){
  
  PAF <- read_paf(path)
  
  prim_alignment <- 
    filter_secondary_alignments(PAF) %>%
    filter(alen > 1000 & mapq > 40) %>%
    unique()
  
  tibble(asm_acc=sub('pESI_screen/(.*).PAF','\\1',path), 
         perc_pESI=sum(prim_alignment$nmatch) / 322518 * 100)
  
}
# tst_paf <- read_paf('pESI_screen/GCA_000230875.1.PAF')

# plot_coverage(tst_paf)

# tst_paf


# test it out...
return_percent_pESI('pESI_screen/GCA_000230875.1.PAF')

# run it on all
# about a minute
TIC <- tic()

pESI_screen_files <- list.files('pESI_screen', pattern = 'PAF', full.names = T)
pESI_presence <- future_map(.x =pESI_screen_files, .f = ~return_percent_pESI(.x) ) %>%
  bind_rows()

TOC <- toc()


pESI_presence %>% 
  filter(perc_pESI > .5) %>%
  ggplot(aes(x=perc_pESI)) +
  geom_histogram(bins=50) + 
  geom_vline(xintercept = 1) + 
  # annotate(geom='text', label='1%', x=5, y=5000) + 
  # geom_vline(xintercept = 25) +
  # geom_vline(xintercept = 75) +
  annotate(geom='text', label='Partial', x=45, y=1000, size=5)+
  geom_vline(xintercept = 75) +
  annotate(geom='text', label='Full', x=100, y=1000, size=5) + 
  theme_half_open() + 
  cowplot::background_grid()





pESI_presence <- 
  pESI_presence %>%
  mutate(pESI_presence=case_when(
    perc_pESI < 1 ~ 'absent',
    perc_pESI > 1 & perc_pESI < 75 ~ 'partial',
    perc_pESI > 75 ~ 'full'
  )) %>% 
  write_tsv('./output/pESI_presence.tsv') 


pESI_asm_accs <- 
  pESI_presence %>% 
  filter(pESI_presence != 'absent') %>% 
  pull(asm_acc)

### function to ID the contig names of pESI containing matches

collect_pESI_contigs <- function(paf_path, genome_fasta_path){
  PAF <- read_paf(paf_path)
  
  fasta <- Biostrings::readDNAStringSet(genome_fasta_path)
  
  matching_contig_IDs <- 
    filter_secondary_alignments(PAF) %>%
    filter(alen > 1000 & mapq > 40) %>%
    unique() %>%
    pull(qname)
  
  result_fasta <- fasta[matching_contig_IDs]
  
  
  return(result_fasta)
  
}


tic()
pESI_contigs <- 
  tibble(fasta_path=list.files('assemblies', pattern = '.fna', full.names = T)) %>% 
  mutate(asm_acc=sub('assemblies/(.*).fna','\\1',fasta_path), 
         paf_path=paste0('pESI_screen/', asm_acc,'.PAF')) %>% 
  filter(asm_acc %in% pESI_asm_accs) %>% 
  mutate(pESI_contigs=future_map2(paf_path, fasta_path, collect_pESI_contigs)) %>% 
  mutate(contigs_path=paste0('pESI_contigs/pESI_', asm_acc, '.fna'))

TOC <- toc()
# 200 seconds

library(usethis)

use_directory('pESI_contigs')

pESI_contigs$pESI_contigs


# write contigs, one per file

# these contigs are probably contaminated with chromosomal contigs
map2(.x=pESI_contigs$pESI_contigs,
     .y=pESI_contigs$contigs_path, .f = ~writeXStringSet(.x, .y))


  


