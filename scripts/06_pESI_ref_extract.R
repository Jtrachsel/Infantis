library(tidyverse)
library(Biostrings)
library(pafr)
library(tictoc)
library(furrr)
library(cowplot)

plan(multisession, workers = 32)  
# pafr::plot_coverage()

# extract only the pESI plasmid sequence from a well known Infantis assembly
FSIS <- readDNAStringSet('reference_genomes/FSIS1502916.fna')
FSIS['FSIS1502916_2'] %>% writeXStringSet('./reference_genomes/pESI.fna')

## ran these commands to map all contigs to the pESI plasmid sequence:
# mkdir -p pESI_screen
# parallel -j 35 'minimap2 -x asm5 -t 1 -o ./pESI_screen/{/.}.PAF ./reference_genomes/pESI.fna {}' ::: ./assemblies/*fna
# mkdir -p pESI_screen2
# system("blastn -query reference_genomes/pESI_CDC_replicon.fna -subject reference_genomes/pESI.fna -outfmt 6")
# parallel -j 35 'blastn -query reference_genomes/pESI_CDC_replicon.fna -subject {} -outfmt 6 > ./pESI_screen2/{/.}.blast' ::: ./assemblies/*fna

return_percent_pESI <- function(path){
  
  PAF <- read_paf(path)
  
  prim_alignment <- 
    filter_secondary_alignments(PAF) %>%
    filter(alen > 1000 & mapq > 40) %>%
    unique()
  
  tibble(asm_acc=sub('pESI_screen/(.*).PAF','\\1',path), 
         perc_pESI=sum(prim_alignment$nmatch) / 322518 * 100)
  
}

blast_col_names <- str_split('qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', pattern = ' ') %>% 
  unlist()


return_pESI_replicon_blast <- function(path){
  blast_col_names <- str_split('qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', pattern = ' ') %>% 
    unlist()
  # browser()
  BLAST <- read_tsv(path,
                    col_names = blast_col_names, 
                    col_types=c('ccdddddddddd'))
  BLAST <- 
    BLAST %>% filter(pident > 90 & length > 500)
  
  res <- tibble(asm_acc=sub('pESI_screen2/(.*).blast','\\1',path), 
                pESI_replicon=ifelse(nrow(BLAST) > 0, 'present', 'absent'))
  return(list(res, BLAST))
  
}

# return_pESI_replicon_blast('pESI_screen2/GCA_006289335.1.blast')

list.files('pESI_screen2', pattern = '*.blast', full.names = TRUE)
TIC <- tic()

pESI_replicon_files <- list.files('pESI_screen2', pattern = 'blast', full.names = T)
pESI_replicon_presence <- future_map(.x =pESI_replicon_files, .f = ~return_pESI_replicon_blast(.x) ) #%>%
  # bind_rows()
replicon_blasts <-
  pESI_replicon_presence %>% 
  map(2) %>% 
  bind_rows() %>% 
  mutate(genome=sub('(GCA_.*)_[0-9]+','\\1',sseqid))



pESI_replicon_presence <- 
  pESI_replicon_presence %>% 
  map(1) %>%
  bind_rows()

TOC <- toc()



# tst_paf <- read_paf('pESI_screen/GCA_000230875.1.PAF')

# plot_coverage(tst_paf)

# tst_paf


# test it out...
# return_percent_pESI('pESI_screen/GCA_000230875.1.PAF')

# run it on all
# about a minute
TIC <- tic()

pESI_screen_files <- list.files('pESI_screen', pattern = 'PAF', full.names = T)
pESI_presence <- future_map(.x =pESI_screen_files, .f = ~return_percent_pESI(.x) ) %>%
  bind_rows()

TOC <- toc()

pESI_replicon_presence <- pESI_replicon_presence %>% left_join(pESI_presence)


pESI_replicon_presence %>% 
  # filter(perc_pESI > .5) %>%
  ggplot(aes(x=perc_pESI, fill=pESI_replicon)) +
  geom_histogram(bins=50) + 
  geom_vline(xintercept = 1) + 
  # annotate(geom='text', label='1%', x=5, y=5000) + 
  # geom_vline(xintercept = 25) +
  # geom_vline(xintercept = 75) +
  annotate(geom='text', label='Partial', x=45, y=1000, size=5)+
  geom_vline(xintercept = 75) +
  annotate(geom='text', label='Full', x=100, y=1000, size=5) + 
  theme_half_open() + 
  cowplot::background_grid() + 
  ylab('number of genomes') + 
  xlab('percent of reference pESI sequence present')

ggsave('pESI_presence_histogram_with_replicon_presence.jpeg', bg='white')
pESI_replicon_presence %>% 
  # filter(perc_pESI > .5) %>%
  ggplot(aes(x=perc_pESI, fill=pESI_replicon)) +
  geom_histogram(bins=50) + 
  # geom_vline(xintercept = 1) + 
  # annotate(geom='text', label='1%', x=5, y=5000) + 
  # geom_vline(xintercept = 25) +
  # geom_vline(xintercept = 75) +
  # annotate(geom='text', label='Partial', x=45, y=1000, size=5)+
  # geom_vline(xintercept = 75) +
  # annotate(geom='text', label='Full', x=100, y=1000, size=5) + 
  theme_half_open() + 
  cowplot::background_grid() + 
  ylab('number of genomes') + 
  xlab('percent of reference pESI sequence present') +
  xlim(0,45) +
  ylim(0,250)

ggsave('pESI_presence_zoom.jpeg', bg='white')

# two problems 
# 1) pESI associated sequences present but replicon absent
# 2) pESI replicon present but very little other pESI associated sequence


no_replicon_but_lots_of_pESI_sequence <- 
  pESI_replicon_presence %>%
  filter(pESI_replicon == 'absent' & perc_pESI > 1) %>% 
  arrange(desc(perc_pESI)) %>% 
  # slice_head(n=10)
  write_tsv('pESI_replicon_absent_lots_of_sequence.tsv')



replicon_but_little_sequence <- 
  pESI_replicon_presence %>%
  filter(pESI_replicon == 'present') %>% 
  filter(perc_pESI < 35) %>% 
  arrange(perc_pESI) %>%
  # slice_head(n=10)
  write_tsv('pESI_replicon_present_little_sequence.tsv')

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

### this garbage down here needs to be run after the 08 plots... script
# Shawn wanted pESI negative isolates, these 10 isoaltes are the only turkey
# FSIS Infantis pESI negative isolates.  
# each has small streches of DNA that match the pESI plasmid sequece,
# this is supposed to extract those sequences
  
LOOK <- read_paf('pESI_screen/GCA_010600745.1.PAF')

LOOK <- LOOK %>%
  filter_secondary_alignments() %>%
  filter(alen > 1000 & mapq > 40) %>%
  unique()

LOOK$qstart
LOOK <- 
  turkey_NOpESI_FSIS %>%
  select(asm_acc) %>% 
  mutate(paf_path=paste0('./pESI_screen/', asm_acc, '.PAF'), 
         PAF=map(.x=paf_path,
                 .f = ~read_paf(.x) %>% filter_secondary_alignments() %>%
                   filter(alen > 1000 & mapq > 40) %>%
                   unique()), 
         qname=map(.x=PAF, .f=~.x %>% pull(qname)), 
         qstart=map(.x=PAF, .f=~.x %>% pull(qstart)), 
         qend=map(.x=PAF, .f=~.x %>% pull(qend))) %>% 
  mutate(genome_seq_path=paste0('assemblies/', asm_acc, '.fna'), 
         genome_seq=map(.x=genome_seq_path, readDNAStringSet),
         pesi_contigs=map2(.x=genome_seq, .y=qname,.f=~.x[.y]), 
         pesi_seqs=pmap(.l = list(pesi_contigs, qstart, qend), 
                        ~subseq(x = ..1, start = ..2, end = ..3)))

LOOK %>%
  pull(pesi_seqs) %>%
  DNAStringSetList() %>% 
  unlist() %>% 
  writeXStringSet('FSIS_infantis_turkey_NOpESI_matches.fna')
