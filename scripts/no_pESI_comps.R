library(tidyverse)
library(pdtools)
library(usethis)


infmeta <- read_tsv('output/08_metadata.tsv')

turkey_NOpESI_FSIS <- 
  infmeta %>%
  filter(grepl('FSIS', collection_agency)) %>%
  filter(pESI_presence == 'absent') %>% 
  filter(ag_match == 'Turkey') 

turkey_NOpESI_FSIS %>% 
  select(asm_acc, strain, geo_loc_name, Year,collection_date, isolation_source, pESI_presence, perc_pESI, stress_genotypes, virulence_genotypes, AMR_genotypes) %>%
  write_tsv('Infantis_FSIS_Turkey_NOpESI.tsv')



turkey_NOpESI_FSIS

# reference strains with pESI
# FSIS12029062 
# FSIS22028788




pesi_refs <- infmeta %>% filter(strain %in% c('FSIS22028788', 'FSIS12029062'))


tmp <- 
  bind_rows(turkey_NOpESI_FSIS, pesi_refs) %>% 
  mutate(fasta_path=map_chr(asm_acc, ~list.files(path='assemblies', pattern=paste0(.x, '.fna'), full.names = T) ))



LOOK <- 
  build_ppanggolin_file_fastas(incomplete_genome_paths = tmp$fasta_path) %>%
  write_tsv('no_pESI_comp.tsv', col_names = FALSE)

# ran the ppanggolin script

use_directory('no_pESI_comp/')

tmp2 <- 
  tmp %>% 
  select(asm_acc, strain, fasta_path) %>%
  mutate(fasta_dest=sub('assemblies/','no_pESI_comp/',fasta_path))


file.copy(tmp2$fasta_path, tmp2$fasta_dest)


# run gifrop

###

library(ggtree)


library(treeio)

tmp$pESI_presence

tmp <- tmp %>% select(asm_acc, everything())

tr <- read.tree('RAxML_bestTree.no_pESI_comp') 
ggtr <- ggtree(tr) 
ggtr <- ggtr %<+% tmp


p1 <- 
  ggtr +
  # geom_tippoint(aes(size=perc_pESI, fill=perc_pESI), shape=21) +
  geom_tiplab(aes(label=strain, color=pESI_presence), hjust = -.1) +
  coord_cartesian(clip='off', xlim = c(0,.00015), ylim=c(-2,17)) 

p1

# p2 <- 
#   ggtr +
#   # geom_tippoint(aes(size=perc_pESI, fill=perc_pESI), shape=21) +
#   geom_tiplab(aes(label=strain, color=pESI_presence), hjust = -.1)



genome_pESI_status <- tmp %>% transmute(genome_name=asm_acc,strain, genome_pESI=pESI_presence)

cii <- read_csv('no_pESI_comp/pan/gifrop_out/clustered_island_info.csv') %>% 
  left_join(genome_pESI_status)

plas_islands <- cii %>% filter(only_island)
chrome_islands <- cii %>% filter(!only_island)


cii$genome_name

cii %>% 
  ggplot(aes(x=genome_pESI, y=quat_cluster, color=island_type)) +
  geom_point()

non_pesi_island_quats <- 
  cii %>% 
  filter(genome_pESI == 'absent') %>%
  pull(quat_cluster) %>%
  unique()


non_pesi_islands <- cii %>% filter(quat_cluster %in% non_pesi_island_quats)

cii %>% filter(quat_cluster == 15) %>% pull(strain)

cii %>% filter(quat_cluster == 15) %>% pull(megares_type)



tmp$stress_genotypes

tmp$strain


LOOK <- tmp %>% filter(strain == 'FSIS22028824')# %>% pull(asm_acc)

LOOK <- cii %>% filter(genome_name == 'GCA_014368505.1')


tmp %>%
  filter(asm_acc == 'GCA_014368505.1') %>% 
  select(asm_acc, fna_download) %>% 
  make_dest_paths('fna', 'no_pESI_comp/GOD_DAMNIT_BRAD') %>% 
  download_genomes(type = 'fna')

LOOK <- read_tsv('no_pESI_comp/GOD_DAMNIT_BRAD/LOOK.tsv')


######## from 08_plots...




AMR_reference <- read_tsv('./metadata/ReferenceGeneCatalog.txt')

# 
# AMR_reference %>% 
#   pivot_longer(cols = c(allele, gene_family), names_to='NAMES_TO', values_to="VALUES_TO")



GENE_2_CLASS_DICT <- 
  AMR_reference %>%
  dplyr::select(allele, gene_family, scope:subclass) %>% 
  pivot_longer(cols = c(allele, gene_family),
               names_to='ORIGIN',
               values_to="GENE_ID") %>% 
  filter(!is.na(GENE_ID)) %>%
  unique()



AMR_finder <- 
  tmp %>% 
  dplyr::select(asm_acc, AMR_genotypes, stress_genotypes, virulence_genotypes) %>% 
  transmute(asm_acc=asm_acc, 
            AMR_ID=gsub('"','',AMR_genotypes), 
            stress_ID=gsub('"','',stress_genotypes),
            virulence_ID=gsub('"','',virulence_genotypes))


# separate_rows(ends_with('ID'), sep = ',') %>% 
# filter(!grepl('PARTIAL', AMR_ID)) %>% 
# filter(!grepl('PARTIAL', virulence_ID)) %>% 
# filter(!grepl('PARTIAL', stress_ID))

# COUNTS OF AMR CLASSES HERE
AMR <- 
  AMR_finder %>%
  transmute(asm_acc=asm_acc, 
            GENE_ID=AMR_ID) %>% 
  separate_rows(GENE_ID, sep = ',') %>% 
  filter(!grepl('PARTIAL', GENE_ID)) %>% 
  mutate(GENE_ID=sub('(.*)=.*', '\\1', GENE_ID)) %>% 
  left_join(GENE_2_CLASS_DICT) 

AMR <- tmp %>%
  dplyr::select(asm_acc,ag_match_lump, Year, country, pESI_presence) %>%
  right_join(AMR)  %>% 
  separate_rows(class, sep = '/') %>%
  unique()

AMR_sum <- 
  AMR %>%
  separate_rows(class, sep = '/') %>% 
  group_by(asm_acc) %>%
  arrange(class) %>%
  summarise(all_AMR_genes=paste(unique(GENE_ID), collapse = '_'), 
            AMR_classes=paste(unique(class[!is.na(class)]), collapse='_'), 
            num_AMR_classes=length(unique(class[!is.na(class)])))


AMR_df <- 
  AMR_sum %>%
  select(asm_acc, AMR_classes) %>% 
  separate_rows(AMR_classes) %>% 
  mutate(presence='TRUE') %>% 
  pivot_wider(names_from = AMR_classes, values_from = presence)%>% 
  column_to_rownames(var='asm_acc') 



#####
STRESS <-  
  AMR_finder %>%
  transmute(asm_acc=asm_acc, 
            GENE_ID=stress_ID) %>% 
  separate_rows(GENE_ID, sep = ',') %>% 
  filter(!grepl('PARTIAL', GENE_ID)) %>% 
  mutate(GENE_ID=sub('(.*)=.*', '\\1', GENE_ID)) %>% 
  left_join(GENE_2_CLASS_DICT) 

STRESS <- tmp %>%
  dplyr::select(asm_acc,ag_match_lump, Year, country, pESI_presence) %>%
  right_join(STRESS)  %>% 
  separate_rows(class, sep = '/') %>%
  unique()


STRESS_sum <- 
  STRESS %>%
  separate_rows(class, sep = '/') %>% 
  group_by(asm_acc) %>%
  arrange(class) %>%
  summarise(all_stress_genes=paste(unique(GENE_ID), collapse = '_'), 
            stress_classes=paste(unique(class[!is.na(class)]), collapse='_'), 
            num_stress_classes=length(unique(class[!is.na(class)])))

stress_df <- 
  STRESS_sum %>%
  select(asm_acc, stress_classes) %>% 
  mutate(stress_classes=sub('QUATERNARY AMMONIUM', 'QACs',stress_classes)) %>% 
  separate_rows(stress_classes, sep = '_') %>% 
  mutate(presence='TRUE') %>% 
  pivot_wider(names_from = stress_classes, values_from = presence) %>% 
  column_to_rownames(var='asm_acc') 



#### virulence

VIR <-  
  AMR_finder %>%
  transmute(asm_acc=asm_acc, 
            GENE_ID=virulence_ID) %>% 
  separate_rows(GENE_ID, sep = ',') %>% 
  filter(!grepl('PARTIAL', GENE_ID)) %>% 
  mutate(GENE_ID=sub('(.*)=.*', '\\1', GENE_ID)) %>% 
  left_join(GENE_2_CLASS_DICT) %>% filter(GENE_ID!='NULL')

VIR <- tmp %>%
  dplyr::select(asm_acc,ag_match_lump, Year, country, pESI_presence) %>%
  right_join(VIR)  %>% 
  separate_rows(class, sep = '/') %>%
  unique()


VIR_sum <- 
  VIR %>%
  separate_rows(class, sep = '/') %>% 
  group_by(asm_acc) %>%
  arrange(class) %>%
  summarise(all_vir_genes=paste(unique(GENE_ID), collapse = '_'), 
            vir_classes=paste(unique(class[!is.na(class)]), collapse='_'), 
            num_vir_classes=length(unique(class[!is.na(class)])))

vir_df <- 
  VIR_sum %>%
  select(asm_acc, all_vir_genes) %>% 
  # mutate(stress_classes=sub('QUATERNARY AMMONIUM', 'QACs',stress_classes)) %>% 
  separate_rows(all_vir_genes, sep = '_') %>% 
  mutate(presence='TRUE') %>% 
  pivot_wider(names_from = all_vir_genes, values_from = presence) %>% 
  column_to_rownames(var='asm_acc') 
vir_df
##### end vir




p2 <- gheatmap(p = p1, data = stress_df, offset=.000030, width=.3,
         colnames_angle=90,
         colnames_position = 'top',
         colnames_offset_y = .00015,
         font.size = 2.5,  hjust = 0, legend_title='presence' )
p2

p3 <- gheatmap(p = p2, data = AMR_df, offset=.000049, width=.55,
         colnames_angle=90,
         colnames_position = 'top',
         colnames_offset_y = .00015,
         font.size = 2.5, hjust = 0, legend_title='presence' )


p3

p4 <- gheatmap(p = p3, data = vir_df, offset=.000079, width=.18,
               colnames_angle=90,
               colnames_position = 'top',
               colnames_offset_y = .00015,
               font.size = 3.5, hjust = 0, legend_title='presence' )


p4


###


# unique(AMR$class)

# AMR %>%
#   select(genome, class) %>% 
#   separate_rows(class, sep = '/') %>% 
#   unique() %>% 
#   mutate(PRESENT = 1) %>% 
#   spread(key = class, value=PRESENT, fill = 0)



