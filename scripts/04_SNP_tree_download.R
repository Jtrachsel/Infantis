library(tidyverse)
library(pdtools)
library(usethis)
library(ape)

use_directory('/90daydata/fsepru113/jtrachsel/Infantis/SNP_trees')

get_PDG_version <- function(data_dir){
  sub('(PDG.*).amr.metadata.tsv','\\1',list.files(data_dir, '(PDG.*).amr.metadata.tsv')[1])
}



download_SNP_trees <- function(data){
  
  original_options <- base::options(timeout = 10000)
  base::on.exit(base::options(original_options))
  
  ### check if SNPS exist here
  data %>% 
    mutate(SNP_tree_dl=map2(.x = SNP_tree_url, .y=SNP_tree_dest, .f = ~download.file(.x, .y)))
}

make_SNPtree_urls(organism = 'Salmonella', 
                  data = infmeta, 
                  PDG = get_PDG_version('metadata'))

make_SNP_tree_dest <- function(data, data_dir){
  data %>% 
    mutate(SNP_tree_dest=
             paste0(data_dir, '/', sub('.*SNP_trees/(PDS[0-9]+.[0-9]+.tar.gz)','\\1',SNP_tree_url)))
}

get_PDG_version('metadata')

infmeta <- read_tsv('metadata/Infantis_metadata.tsv')

SNP_tree_download_info <- 
  infmeta %>% 
  group_by(PDS_acc) %>% 
  tally() %>%
  filter(n>10) %>% 
  make_SNPtree_urls(organism = 'Salmonella',
                    PDG = get_PDG_version('metadata')) %>% 
  enframe(name=NULL, value='SNP_tree_url') %>% 
  make_SNP_tree_dest(data_dir='./SNP_trees/') %>% 
  download_SNP_trees()


SNP_tree_download_info %>% 
  unnest(SNP_tree_dl) %>%
  write_tsv('./metadata/SNP_tree_info.tsv') %>% 
  filter(SNP_tree_dl != 0)


## untar trees here... ###



