library(tidyverse)

infmeta <- read_tsv('metadata/Infantis_metadata.tsv') %>% 
  left_join(read_tsv('output/pESI_presence.tsv')) %>% 
  left_join(read_tsv('output/phydelity_clusters.tsv')) %>% 
  mutate(pESI_presence = factor(pESI_presence, levels = c('absent', 'partial','full')))


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



### then follow up with others showing the breakdowns by various categories

# read in outputs and merge into metadata
infmeta %>% 
  mutate(PDS_lump=fct_lump_n(f = PDS_acc, n = 8), 
         ag_match_lump=fct_lump_n(f=ag_match, n=4)) %>% 
  count(Year,ag_match_lump, PDS_lump) %>% 
  filter(Year < 2022) %>% 
  filter(Year > 2015) %>% 
  # filter(country_lump )
  ggplot(aes(x=Year, y=n, color=PDS_lump)) +
  geom_line(size=1.5, color='black', aes(group=PDS_lump))+
  geom_line(size=1.2) +
  ylab('number of genomes') + 
  guides(color=guide_legend(title="SNP cluster")) +
  scale_color_brewer(palette = 'Set1') + 
  theme_bw() + facet_wrap(~ag_match_lump, scales = 'free_y') 

### split into US and non-US

infmeta %>%
  mutate(PDS_lump=fct_lump_n(f = PDS_acc, n = 8), 
         ag_match_lump=fct_lump_n(f=ag_match, n=6)) %>% 
  filter(country != 'USA') %>% 
  count(Year,ag_match_lump, PDS_lump) %>% 
  filter(Year < 2022) %>% 
  filter(Year > 2015) %>% 
  # filter(country_lump )
  ggplot(aes(x=Year, y=n, color=PDS_lump)) +
  geom_line(size=1.5, color='black', aes(group=PDS_lump))+
  geom_line(size=1.2) +
  ylab('number of genomes') + 
  guides(color=guide_legend(title="SNP cluster")) +
  scale_color_brewer(palette = 'Set1') + 
  theme_bw() + facet_wrap(~ag_match_lump, scales = 'free_y') 





infmeta %>% count(Year, PDS_acc) %>% arrange(desc(n))
infmeta %>% count(country, PDS_acc) %>% arrange(desc(n))

# UK strains?
# one isolate from pork in denmark 2018, the rest are from 'food' or 'environmental' isolations
# in the UK
look <- 
  infmeta %>%
  filter(PDS_acc == 'PDS000073471.105') %>% 
  filter(ag_match != 'Human') %>%
  select(Year,country,host, ag_match)

###

# pESI presence

infmeta %>%
  group_by(pESI_presence, ag_match) %>% 
  tally() %>% 
  arrange(desc(n))



infmeta %>%
  group_by(pESI_presence, country) %>% 
  tally() %>% 
  arrange(desc(n))

# by ag species
### THIS ONE!!!
infmeta %>%
  mutate(ag_match_lump=fct_lump_n(f=ag_match, n=6)) %>% 
  group_by(pESI_presence, Year, ag_match_lump) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(ag_match_lump, Year) %>% 
  mutate(tot_year=sum(n), 
         perc_year=n/tot_year) %>% 
  filter(Year < 2022) %>% 
  filter(Year > 2005) %>% 
  arrange(desc(n)) %>% 
  ggplot(aes(x=Year, y=perc_year, fill=pESI_presence)) +
  geom_col()+
  # geom_line(size=2.2, color='black') + 
  # geom_line(size=2) + 
  facet_wrap(~ag_match_lump, scales = 'free')

# by PDS accession #


infmeta %>%
  mutate(PDS_lump=fct_lump_n(f=PDS_acc, n=8)) %>% 
  group_by(pESI_presence, Year, PDS_lump) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(Year, PDS_lump) %>% 
  mutate(perc_year=n/sum(n)) %>% 
  filter(Year < 2022) %>% 
  filter(Year > 2010) %>% 
  arrange(desc(n)) %>% 
  ggplot(aes(x=Year, y=perc_year, fill=pESI_presence, group=pESI_presence)) +
  # geom_line(size=2.2, color='black') + 
  geom_col()+
  # geom_text(aes(label=(n)))+
  # geom_line(size=2) + 
  facet_wrap(~PDS_lump)

full_pESI_PDS <- 
  infmeta %>% 
  group_by(PDS_acc, pESI_presence) %>% 
  tally() %>% 
  filter(pESI_presence == 'full') %>% 
  # filter(n >5) %>% 
  pull(PDS_acc) %>% unique()

infmeta %>% 
  filter(PDS_acc %in% full_pESI_PDS) %>% 
  group_by(PDS_acc, pESI_presence) %>% 
  tally() %>%
  ungroup() %>% 
  group_by(PDS_acc) %>% 
  mutate(tot_PDS=sum(n), 
         perc_PDS=n/tot_PDS) %>% 
  # filter(pESI_presence != 'absent')
  ggplot(aes(y=PDS_acc, x=perc_PDS, fill=pESI_presence)) + 
  geom_col() + 
  ggtitle('percent pESI plasmid presence in SNP cluster with at least one "full" genome ')





absent_pESI_PDS <- 
  infmeta %>% 
  group_by(PDS_acc, pESI_presence) %>% 
  tally() %>% 
  filter(pESI_presence == 'absent') %>% 
  # filter(n >5) %>% 
  pull(PDS_acc) %>% unique()

infmeta %>% 
  filter(PDS_acc %in% absent_pESI_PDS) %>% 
  group_by(PDS_acc, pESI_presence) %>% 
  tally() %>%
  ungroup() %>% 
  group_by(PDS_acc) %>% 
  mutate(tot_PDS=sum(n), 
         perc_PDS=n/tot_PDS) %>% 
  # filter(pESI_presence != 'absent')
  ggplot(aes(y=PDS_acc, x=perc_PDS, fill=pESI_presence)) + 
  geom_col() + 
  ggtitle('percent pESI plasmid presence in SNP cluster with at least one "absent" genome ')










#########
# select representatives for overview tree
########



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




#### orphaned from transmission summary script ###




phydel_sum <- 
  meta %>% 
  filter(!is.na(phy_clust)) %>% 
  group_by(PDS_acc, phy_clust) %>% 
  summarise(num_isolates=n(),
            hosts=paste(sort(unique(ag_match)), collapse = '_'), 
            years=paste(sort(unique(Year)), collapse = '_'), 
            countries=paste(sort(unique(country)), collapse = '_')) %>% 
  arrange(desc(num_isolates)) %>% 
  mutate(PDS_size=sum(num_isolates), 
         prop_tot_PDS=num_isolates/PDS_size)
phydel_sum
hist(phydel_sum$prop_tot_PDS)


phydel_sum$hosts %>% unique()



#### pESI presence determination ####


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



### should run stats against PESI pres, 
glm(family='binomial', pesi_prez ~ Year + ag_host + country, data=meta)

library(tidyverse)

ex_dat <- tibble(genome=paste0('genome_', seq(75)),
                 year_category=sample(c('Before 2017', 'After 2017'),replace=T, size = 75),
                 phage_presence=sample(c(0,1), size = 75, replace = T))



ex_dat


### Fisher exact test
fish_test <-
  fisher.test(
    table(ex_dat$year_category, ex_dat$phage_presence)
  )

fish_test


### logistic regression
log_reg <- glm(data=ex_dat, formula=phage_presence ~ year_category, family = 'binomial')
summary(log_reg)

broom::tidy(log_reg, conf.int=T, exponentiate=T)



### fisher's exact gives very similar results to logistic regression
# logistic regression odds ratio:
broom::tidy(log_reg, conf.int=T, exponentiate=T)$estimate[2]

# fisher odds ratio:
fish_test$estimate












































## EXPERIMENTAL ZONE  ###

##### gene presence/absence variability within clusters
# given a distance matrix and a set of isolates
# calculate cumulative pairwise distances divided by the number of isolates
#   mean_pw_dist = clusters with higher mean_pw_dist have isolates that have gained/lost genes in the transmission event?




trans_meta <- meta %>% filter(!is.na(phy_clust))


gpa <- read_tsv('all_pan/WRITE/gene_presence_absence.Rtab')

gpa <- gpa %>% column_to_rownames(var='Gene') %>% as.matrix()

gpa <- gpa[,colnames(gpa) %in% trans_meta$asm_acc]
gpa <- gpa[rowSums(gpa) > 0,]



meta <- meta %>% filter(asm_acc %in% colnames(gpa))

trans_meta <- trans_meta %>% filter(asm_acc %in% colnames(gpa))
# library(parallelDist)
# library(usedist)

# trans_dists <- parallelDist(t(gpa), method = 'binary')

# write_rds(trans_dists,'output/transmission_cluster_dists.rds')

# attributes(trans_dists)$Labels

# calc_gene_variability <- function(dist_mat){
#   # 1) subset dist to only contain the asm_accs
#   # 2) move to long format
#   # 3) sum(dists) / n()
#   # browser()
#   from_centroid <- dist_to_centroids(dist_mat, rep(1, length(attributes(dist_mat)$Labels)))
#   return(from_centroid)
#   
#   
# }

# sum(trans_dists) / attributes(trans_dists)$Size
# attributes(trans_dists)$Size
# length(trans_dists)

extract_PA <- function(asm_accs, pan_PA){
  # browser()
  pan_PA <- pan_PA[,colnames(pan_PA) %in% asm_accs]
  pan_PA <- pan_PA[rowSums(pan_PA) > 0,]
  return(pan_PA)
}


extract_variable_genes <- function(gpa){
  # browser()
  loggpa <- gpa > 0
  res <- apply(X = loggpa, MARGIN=1, all)
  variable_genes <- names(res[!res])
  var_gene_PA <- gpa[rownames(gpa) %in% variable_genes,]
  
}

test <- 
  trans_meta %>% 
  group_by(phy_clust) %>%
  nest() %>% 
  left_join(phydel_sum) %>% 
  filter(num_isolates > 10) %>% 
  ungroup()

test1 <- 
  test %>% 
  mutate(clust_PA=map(.x=data, ~extract_PA(pan_PA = gpa, asm_accs = .x$asm_acc)), 
         clust_dist=map(clust_PA, ~dist(.x, method = 'binary'))) 


# test2 <- test1 %>% mutate(dists_to_centroids=map(clust_dist, calc_gene_variability))

test2 <- test1 %>%
  mutate(var_genes=map(clust_PA, extract_variable_genes), 
         num_var_genes=map_int(var_genes, ~dim(.x)[1])) %>% 
  arrange(desc(num_var_genes))

test2 %>% 
  mutate(num_years=map_int(.x=strsplit(years, '_'), .f = length)) %>% 
  ggplot(aes(x=num_years, y=num_var_genes)) + geom_point() + 
  geom_smooth(method = 'lm')

test2 %>% filter(phy_clust == 1248) 

test %>% filter(phy_clust == 919) %>% 
  mutate(clust_PA=map(.x=data, ~extract_PA(pan_PA = gpa, asm_accs = .x$asm_acc)))



LOOK <- meta %>% filter(phy_clust == 919)
usedist::dist_to_centroids()


which.min(test$num_isolates)
test %>% filter(num_isolates == 1)
gpa
.GlobalEnv$gpa
trans_meta %>% filter(phydel)

dist_subset(trans_dists, trans_meta[1:10]$asm_acc)
all(attributes(trans_dists)$Labels %in% trans_meta$asm_acc)


all(trans_meta$asm_acc %in% attributes(trans_dists)$Labels)

trans_meta <- trans_meta %>% filter(asm_acc %in% attributes(trans_dists)$Labels)

m4 <- matrix(1:16, nrow=4, dimnames=list(LETTERS[1:4]))
dm4 <- dist(m4)
dist_subset(dm4, c("A", "B", "C"))
dist_subset(dm4, c("D", "C", "B", "A"))

attributes(dm4)

TESTTEST <- c('A', 'B', 'C', 'D')
dist_subset(dm4, TESTTEST)




#######



rgps <- read_tsv('./all_pan/WRITE/plastic_regions.tsv')

modules <- read_tsv('./all_pan/WRITE/modules_summary.tsv')


