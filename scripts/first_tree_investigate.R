library(tidyverse)
library(ggtree)
library(treeio)
library(DECIPHER)
library(Biostrings)

CORE_ALN <- readDNAStringSet('overview_tree/MSA/core_genome_alignment.aln', 'fasta')

CORE_DIST2 <- DistanceMatrix(CORE_ALN, 
                            penalizeGapLetterMatches = T,
                            correction = 'Jukes-Cantor', processors = 8)

cumulative_core_dist2 <- 
  CORE_DIST2 %>%
  as.data.frame() %>% 
  rownames_to_column(var='asm_acc') %>% 
  pivot_longer(cols=-asm_acc, names_to = 'to', values_to='phydist') %>% 
  group_by(asm_acc) %>% 
  summarise(total_distance=sum(phydist)) %>% 
  arrange(desc(total_distance)) 

hist(cumulative_core_dist$total_distance, breaks = 100)

cumulative_core_dist %>%
  mutate(BAD_REP=ifelse(asm_acc %in% BAD_REPS, TRUE, FALSE)) %>% 
  ggplot(aes(x=total_distance, fill=BAD_REP)) + geom_histogram()


cumulative_core_dist2 %>%
  mutate(BAD_REP=ifelse(asm_acc %in% BAD_REPS, TRUE, FALSE)) %>% 
  ggplot(aes(x=total_distance, fill=BAD_REP)) + geom_histogram()


# hist(phylo_dist$phydist)

cumulative_dists <- 
  phylo_dist %>%
  group_by(asm_acc) %>% 
  summarise(total_distance=sum(phydist)) %>% 
  arrange(desc(total_distance)) 



### read raxml tree ###
tr <- read.raxml('./RAxML_bipartitionsBranchLabels.overview_tree')

ggtree(tr) #+# geom_nodelab(aes(label=node))+
  # geom_hilight(node=525) + geom_highlight(node=864, fill='green') + geom_highlight(node=686, fill='red')

library(ape)
# coph_dists <- ape::cophenetic.phylo(tr@phylo)
# 
# tst <- 
#   coph_dists %>%
#   as.data.frame() %>%
#   rownames_to_column(var='asm_acc') %>% 
#   pivot_longer(cols = -asm_acc, names_to = 'to', values_to = 'phydist') %>% 
#   group_by(asm_acc) %>% 
#   summarise(TOT=sum(phydist)) %>% 
#   arrange(desc(TOT))

phylo_dist <- 
  ape::cophenetic.phylo(tr@phylo) %>%
  as.data.frame() %>% 
  rownames_to_column(var='asm_acc') %>% 
  pivot_longer(cols=-asm_acc, names_to = 'to', values_to='phydist')

# hist(phylo_dist$phydist)

cumulative_dists <- 
  phylo_dist %>%
  group_by(asm_acc) %>% 
  summarise(total_distance=sum(phydist)) %>% 
  arrange(desc(total_distance)) 

hist(cumulative_dists$total_distance, breaks = 100)



BAD_REPS <- cumulative_dists %>% filter(total_distance > 1) %>% pull(asm_acc)


tree_rep_dat <- read_tsv('output/tree_rep_genomes.tsv')


tree_rep_dat %>% filter(rep_genome %in% BAD_REPS)


sum(coph_dists > .001)

FSIS_USDA <- 
  meta %>%
  filter(grepl('FSIS|USDA', collection_agency))

FSIS_PDSs <- FSIS_USDA %>% pull(PDS_acc) %>% unique()


PDS_summary <- 
  meta %>% 
  filter(!is.na(PDS_acc)) %>% 
  filter(Year > 2017) %>% 
  group_by(PDS_acc, Year, country, ag_match) %>%
  tally() %>%
  arrange(desc(n)) %>%
  ungroup() %>% 
  group_by(PDS_acc) %>% 
  # nest()
  summarise(num_years=length(unique(Year)),
            total_isolates=sum(n),
            most_recent_year=unique(Year[which.max(Year)]), 
            most_recent_year_num=sum(n[which.max(Year)]), 
            previous_year_num = sum(n[which(Year == Year[which.max(Year)] - 1)]),
            all_ag_matches=paste(unique(sort(ag_match)), collapse = '_'),
            .groups='drop') %>%
  mutate(hosts_except_human=sub('Human', '', all_ag_matches), 
         hosts_except_human=sub('__','_',hosts_except_human), 
         hosts_except_human=sub('_$','',hosts_except_human), 
         hosts_except_human=sub('^_','',hosts_except_human)) %>% 
  arrange(desc(num_years)) %>% 
  mutate(USDA=ifelse(PDS_acc %in% FSIS_PDSs, 'TRUE', 'FALSE'))


tree_data <- 
  tibble(asm_acc=tr@phylo$tip.label) %>%
  left_join(meta) %>% 
  left_join(PDS_summary) %>% 
  mutate(LABEL=asm_acc)

tree_data$lineage

ggtr <- ggtree(tr) %<+% tree_data
# tree_data$LABEL

P_O157_1 <- 
  ggtr +
  geom_hilight(node=525, alpha=.25) +
  geom_highlight(node=864, fill='green', alpha=.25) +
  geom_highlight(node=686, fill='red', alpha=.25)+
  geom_tippoint(aes(fill=lineage, alpha=Year, size=total_isolates), shape=21) + 
  # geom_point2(aes(subset= FSIS == 'TRUE'),position = position_nudge(x = .00009, y = 0), color='red', size=1) + 
  geom_text2(aes(subset= !grepl('GCA',LABEL), label=LABEL), nudge_x = .00005, size=2) + 
  geom_cladelabel(node=525, label='clade_one') + 
  geom_cladelabel(node=864, label='clade_two') + 
  geom_cladelabel(node=686, label='clade_three')

P_O157_1
