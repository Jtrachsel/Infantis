library(tidyverse)

### Since each of these files was calculated for each SNP cluster in isolation
### the numbers describing the cluster identity are not unique
### need to rename the clusters with SNP clusters in mind
phydelity_files <- list.files('./output/phydelity_results/', full.names = T)


read_tsv(phydelity_files[1])

parse_phydelity <- 
  function(file){
  read_tsv(file) %>% 
    transmute(target_acc=TAXA, 
              phydel_clust=CLUSTER)
  }

phydel_dat <- lapply(phydelity_files, parse_phydelity) %>% bind_rows()


meta <- read_tsv('metadata/Infantis_metadata.tsv') %>% filter(!is.na(ftp_path))



clean_clusts <- 
  meta %>% 
  left_join(phydel_dat) %>%
  filter(!is.na(phydel_clust)) %>% 
  group_by(PDS_acc, phydel_clust) %>% 
  mutate(phy_clust=cur_group_id()) %>% 
  ungroup() %>% 
  select(asm_acc, phy_clust) %>%
  arrange(phy_clust) %>% 
  write_tsv('./output/phydelity_clusters.tsv')

meta <- meta %>% left_join(clean_clusts)


####### end
