library(tidyverse)
library(usethis)
library(pdtools)

# remotes::install_github('jtrachsel/pdtools', ref = 'dev', force = TRUE)
# setwd('/project/fsepru113/trachsel/Infantis')
# usethis::proj_sitrep()

use_directory('./metadata/')
use_directory('./assemblies/')
use_directory('./output/')
use_directory('./scripts/')


if (!file.exists('./metadata/assembly_summary.txt')){
  gbk <- download_gbk_assembly_summary('./metadata/assembly_summary.txt')
}

PDG_files <- list.files(path = './metadata', pattern = 'PDG', full.names = T)

if (length(PDG_files) == 0){
  download_most_recent_complete('Salmonella', folder_prefix = './metadata/')
  PDG_files <- list.files(path = './metadata', pattern = 'PDG', full.names = T)
}
PDG_files <- list.files(path = './metadata', pattern = 'PDG', full.names = T)
PDG_files <- sort(PDG_files, decreasing = TRUE)
PDG_files
all_sal <- 
  read_tsv(PDG_files[2]) %>% 
  left_join(read_tsv(PDG_files[1])) %>% 
  mutate(SERO=sub('serotype=(.*),antigen_formula=(.*)','\\1',computed_types))


Infantis_PDS <- 
  all_sal %>% 
  filter(SERO=='Infantis') %>% 
  pull(PDS_acc) %>% 
  unique() %>% 
  na.omit()

Infantis_related <- 
  all_sal %>% 
  filter(PDS_acc %in% Infantis_PDS) %>%
  filter(!is.na(asm_acc)) %>% 
  filter(!is.null(asm_acc)) %>% 
  filter(asm_acc != 'NULL')

year_data <- Infantis_related %>% extract_earliest_year()
host_data <- Infantis_related %>% extract_consensus_ag_species()
country_data <- Infantis_related %>% extract_country()
agency_data <- Infantis_related %>% extract_collection_agency()



# agency_data$collection_agency %>% table()
# year_data$Year %>% table()
# country_data$country %>% table()
# host_data$ag_match %>% table()
# these <- host_data %>% filter(ag_match == 'Bovine_Chicken') %>% pull(target_acc)
# LOOK <- Infantis_related %>% filter(target_acc %in% these)


Infantis_metadata <- 
  Infantis_related %>% 
  left_join(year_data) %>% 
  left_join(host_data) %>% 
  left_join(country_data) %>% 
  left_join(agency_data) %>% 
  make_ftp_paths(assembly_summary_path = 'metadata/assembly_summary.txt') %>% 
  make_download_urls('fna') %>% 
  make_dest_paths(type = 'fna', dest_dir = './assemblies/') %>% 
  make_download_urls(type = 'gff') %>% 
  make_dest_paths(type = 'gff', dest_dir = './assemblies/') %>% 
  write_tsv('./metadata/Infantis_metadata.tsv')

#Infantis_metadata %>% count(asm_acc) %>% arrange(desc(n))

library(future)
future::plan(multisession, workers=4)

Infantis_metadata <- 
  Infantis_metadata %>% 
  download_genomes(type = 'fna', PARALLEL = TRUE) %>% 
  write_tsv('./metadata/Infantis_metadata.tsv')


Infantis_metadata %>% count(fna_exists)

### The 'type' argument of download_genomes doesnt work right,
### checks for 'fna' files even though gff specified.
Infantis_metadata <- 
  Infantis_metadata %>% 
  download_genomes(type='gff', PARALLEL = TRUE) %>% 
  write_tsv('./metadata/Infantis_metadata.tsv')
  
Infantis_metadata %>% count(gff_exists)








### gunzip all files ###

# remove those that dont gunzip properly


# LOOK <- Infantis_metadata %>% filter(!gff_exists)

# 
# extract_computed_serotype <- function(data){
#   mutate(SERO=sub('serotype=(.*),antigen_formula=(.*)','\\1',computed_types), 
#          ant_form=sub('serotype=(.*),antigen_formula=(.*)','\\2',computed_types))
# }
# 
# extract_continent <- function(data){
#   
# }


### add referene genome download here



