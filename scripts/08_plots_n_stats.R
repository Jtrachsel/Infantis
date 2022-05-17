library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(countrycode)
library(pdtools)

infmeta <-
  read_tsv('metadata/Infantis_metadata.tsv') %>% 
  left_join(read_tsv('output/pESI_presence.tsv')) %>% 
  left_join(read_tsv('output/phydelity_clusters.tsv')) %>% 
  mutate(pESI_presence = factor(pESI_presence, levels = c('absent', 'partial','full')))


infmeta %>%
  count(phy_clust, pESI_presence) %>%
  group_by(phy_clust) %>% tally() %>% 
  filter(n !=1) %>% 
  arrange(desc(n))






country_vector

infmeta$continent <- 
  countrycode(sourcevar = infmeta$country,
              origin = "country.name",
              destination = "continent")



###

# lump SNP cluster levels
# lump ag_match levels

infmeta <- 
  infmeta %>% 
  mutate(PDS_lump=fct_lump_n(f = PDS_acc, n = 8), 
         ag_match_lump=fct_lump_n(f=ag_match, n=6)) 

infmeta %>% write_tsv('output/08_metadata.tsv')
# set SNP cluster levels
PDS_colors <- brewer.pal(9, 'Set1')
names(PDS_colors) <- infmeta$PDS_lump %>% levels()


### This one!
# Overall SNP cluster abundance
P_inf_1 <- 
  infmeta %>% 
  # mutate(PDS_lump=fct_lump_n(f = PDS_acc, n = 8)) %>% 
  count(Year, PDS_lump) %>% 
  filter(Year < 2022) %>% 
  filter(Year > 2015) %>% 
  ggplot(aes(x=Year, y=n, color=PDS_lump)) +
  geom_line(size=1.5, color='black', aes(group=PDS_lump))+
  geom_line(size=1.2) +
  ylab('number of genomes') + 
  guides(color=guide_legend(title="SNP cluster")) +
  scale_color_manual(values = PDS_colors) + 
  theme_bw() + 
  guides(fill=guide_legend(title='SNP cluster')) +
  ggtitle('Major Infantis SNP clusters over time') + 
  theme_half_open() + 
  background_grid()
P_inf_1


### then follow up with others showing the breakdowns by various categories

# read in outputs and merge into metadata
P_inf_2 <- infmeta %>% 
  # mutate(PDS_lump=fct_lump_n(f = PDS_acc, n = 8), 
  #        ag_match_lump=fct_lump_n(f=ag_match, n=6)) %>% 
  count(ag_match_lump, PDS_lump) %>% 
  # filter(Year < 2022) %>% 
  # filter(Year > 2015) %>% 
  # filter(country_lump )
  ggplot(aes(x=ag_match_lump, y=n, fill=PDS_lump)) +
  # geom_line(size=1.5, color='black', aes(group=PDS_lump))+
  # geom_line(size=1.2) +
  ylab('number of genomes') +
  geom_col(color='black')+
  guides(color=guide_legend(title="SNP cluster")) +
  scale_fill_manual(values = PDS_colors) + 
  theme_bw()  + 
  xlab('Host') + 
  guides(fill=guide_legend(title='SNP cluster')) +
  ggtitle('Major Infantis SNP clusters by host')+ 
  theme_half_open()+
  background_grid()

P_inf_2
### split into US and non-US



# US only
P_inf_3 <- 
  infmeta %>%
  filter(country == 'USA') %>% 
  count(ag_match_lump, PDS_lump) %>% 
  ggplot(aes(x=ag_match_lump, y=n, fill=PDS_lump)) +
  geom_col(color='black')+
  ylab('number of genomes') + 
  guides(fill=guide_legend(title="SNP cluster")) +
  scale_fill_manual(values = PDS_colors) + 
  theme_half_open()+
  background_grid()+
  ggtitle('US isolates only')

P_inf_3



P_inf_4 <- 
  infmeta %>%
  filter(country != 'USA') %>% 
  count(ag_match_lump, PDS_lump) %>% 
  ggplot(aes(x=ag_match_lump, y=n, fill=PDS_lump)) +
  geom_col(color='black')+
  ylab('number of genomes') + 
  guides(fill=guide_legend(title="SNP cluster")) +
  scale_fill_manual(values = PDS_colors) + 
  theme_half_open()+
  background_grid()+
  # facet_wrap(~ag_match_lump, scales = 'free_y') +
  ggtitle('Non-US')

P_inf_4

# UK strains?
# one isolate from pork in denmark 2018, the rest are from 'food' or 'environmental' isolations
# in the UK
# what's the deal with the yellow SNP cluster?
# UK human isolates mostly.  One swine isolate from Denmark.

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

### NEED HISTOGRAM

# raw histogram
P_inf_5 <- 
  infmeta %>% 
  # filter(perc_pESI > .5) %>%
  ggplot(aes(x=perc_pESI)) +
  geom_histogram(bins=50) + 
  # geom_vline(xintercept = 1) + 
  # annotate(geom='text', label='1%', x=5, y=5000) + 
  # geom_vline(xintercept = 25) +
  # geom_vline(xintercept = 75) +
  # annotate(geom='text', label='Partial', x=45, y=1000, size=5)+
  # geom_vline(xintercept = 75) +
  # annotate(geom='text', label='Full', x=100, y=1000, size=5) + 
  theme_half_open() + 
  xlab('Percent reference pESI plasmid coverage')+
  cowplot::background_grid()
P_inf_5


P_inf_6 <- 
  infmeta %>% 
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
  xlab('Percent reference pESI plasmid coverage')+
  ylab('Count') +
  cowplot::background_grid()

P_inf_6





# by ag species
### THIS ONE!!!
P_inf_7 <- 
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
  geom_col(color='white')+
  # geom_line(size=2.2, color='black') + 
  # geom_line(size=2) + 
  facet_wrap(~ag_match_lump, scales = 'free')

P_inf_7

# by PDS accession #


# infmeta %>%
#   mutate(PDS_lump=fct_lump_n(f=PDS_acc, n=8)) %>% 
#   group_by(pESI_presence, Year, PDS_lump) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   group_by(Year, PDS_lump) %>% 
#   mutate(perc_year=n/sum(n)) %>% 
#   filter(Year < 2022) %>% 
#   filter(Year > 2010) %>% 
#   arrange(desc(n)) %>% 
#   ggplot(aes(x=Year, y=perc_year, fill=pESI_presence, group=pESI_presence)) +
#   # geom_line(size=2.2, color='black') + 
#   geom_col()+
#   # geom_text(aes(label=(n)))+
#   # geom_line(size=2) + 
#   facet_wrap(~PDS_lump)

full_pESI_PDS <- 
  infmeta %>% 
  group_by(PDS_acc, pESI_presence) %>% 
  tally() %>% 
  filter(pESI_presence == 'full') %>% 
  # filter(n >5) %>% 
  pull(PDS_acc) %>% unique()

P_inf_8 <- 
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
  theme(axis.text.y = element_blank())+
  ggtitle('percent pESI plasmid presence in SNP cluster with at least one "full" genome ')

P_inf_8




absent_pESI_PDS <- 
  infmeta %>% 
  group_by(PDS_acc, pESI_presence) %>% 
  tally() %>% 
  filter(pESI_presence == 'absent') %>% 
  # filter(n >5) %>% 
  pull(PDS_acc) %>% unique()

P_inf_9 <- 
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
  theme(axis.text.y = element_blank())+
  ggtitle('percent pESI plasmid presence in SNP cluster with at least one "absent" genome ')

P_inf_9

### SNP clusters basically either have the full plasmid or dont


#### orphaned from transmission summary script ###


### REMOVED SOME AG MATCHES BECAUSE OF UNDERSCORES

phydel_sum <- 
  infmeta %>% 
  filter(!grepl('_', ag_match)) %>% 
  filter(!is.na(phy_clust)) %>% 
  group_by(PDS_acc, phy_clust) %>% 
  summarise(num_isolates=n(),
            hosts=paste(sort(unique(ag_match)), collapse = '_'), 
            years=paste(sort(unique(Year)), collapse = '_'), 
            continents=paste(sort(unique(continent)), collapse = '_'), 
            pESI=paste(sort(unique(pESI_presence)))) %>% 
  arrange(desc(num_isolates)) %>% 
  ungroup() %>% 
  group_by(PDS_acc) %>% 
  mutate(PDS_size=sum(num_isolates), 
         prop_tot_PDS=num_isolates/PDS_size) 



major_hosts <- 
  phydel_sum %>%
  group_by(hosts) %>%
  summarise(tot_isos=sum(num_isolates)) %>% 
  arrange(desc(tot_isos)) %>% 
  filter(tot_isos >100) %>% 
  pull(hosts)





P_inf_10 <- 
  phydel_sum %>%
  filter(hosts %in% major_hosts) %>% 
  filter(num_isolates >3) %>% 
  # filter(num_isolates > 7 ) %>% 
  # mutate(hosts2=fct_count(hosts, sort = T)) %>% 
  ggplot(aes(y=fct_reorder(hosts, num_isolates, .fun = sum ), x=num_isolates, fill=continents)) +
  geom_col(color='white', size=.1) +
  theme_half_open() + 
  background_grid() + 
  ylab('Hosts within transmission cluster') + 
  xlab('Number of isolates') + 
  ggtitle('Predicted transmission clusters by host')
P_inf_10



#### pESI presence determination ####

# 
# pESI_presence %>% 
#   filter(perc_pESI > .5) %>%
#   ggplot(aes(x=perc_pESI)) +
#   geom_histogram(bins=50) + 
#   geom_vline(xintercept = 1) + 
#   # annotate(geom='text', label='1%', x=5, y=5000) + 
#   # geom_vline(xintercept = 25) +
#   # geom_vline(xintercept = 75) +
#   annotate(geom='text', label='Partial', x=45, y=1000, size=5)+
#   geom_vline(xintercept = 75) +
#   annotate(geom='text', label='Full', x=100, y=1000, size=5) + 
#   theme_half_open() + 
#   cowplot::background_grid()
# 
# 



##### AMR #######


tree_clade_year_iso_totals <-
  infmeta %>%
  filter(Year < 2022) %>% 
  count(pESI_presence, ag_match_lump, Year)

tree_clade_iso_totals <- 
  infmeta %>%
  count(pESI_presence, ag_match_lump, .drop=F,name="tot_n_genomes") 



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
  infmeta %>% 
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

AMR <- infmeta %>%
  dplyr::select(asm_acc,ag_match_lump, Year, country, pESI_presence) %>%
  right_join(AMR)  %>% 
  separate_rows(class, sep = '/') %>%
  unique()



# unique(AMR$class)

# AMR %>%
#   select(genome, class) %>% 
#   separate_rows(class, sep = '/') %>% 
#   unique() %>% 
#   mutate(PRESENT = 1) %>% 
#   spread(key = class, value=PRESENT, fill = 0)



AMR_sum <- 
  AMR %>%
  separate_rows(class, sep = '/') %>% 
  group_by(asm_acc) %>%
  arrange(class) %>%
  summarise(all_AMR_genes=paste(unique(GENE_ID), collapse = '_'), 
            AMR_classes=paste(unique(class[!is.na(class)]), collapse='_'), 
            num_AMR_classes=length(unique(class[!is.na(class)])))


Percent_resist_ISO_clade <- 
  AMR %>%
  dplyr::select(asm_acc, ag_match_lump, class, pESI_presence) %>%
  unique() %>%
  # mutate(class=factor(class), 
  #        tree_clade=factor(tree_clade)) %>% 
  count(class, ag_match_lump, pESI_presence, .drop=F) %>%
  left_join(tree_clade_iso_totals) %>% 
  filter(class != 'EFFLUX') %>%
  mutate(prop_genomes=n/tot_n_genomes, 
         perc_genomes=prop_genomes *100) 



mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(16)


names(mycolors) <- unique(Percent_resist_ISO_clade$class)



###### THIS ONE
P_inf_11 <-
  Percent_resist_ISO_clade %>% 
  filter(pESI_presence == 'full') %>%
  ggplot(aes(x=class, y=perc_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=50, label=round(perc_genomes, 1))) +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~ag_match_lump+ tot_n_genomes, nrow=1) +
  scale_y_continuous(limits = c(0,100),breaks = c(0,50,100))+  coord_flip()  +
  theme(axis.text.x=element_text(size=10, angle = -45, hjust = 0), 
        legend.position = 'none') + 
  ylab('Percent of Genomes') + 
  xlab('') + 
  theme_half_open()+
  background_grid()+
  ggtitle('AMR resistance for genomes with full pESI plasmid')

P_inf_11

###### THIS ONE
P_inf_12 <-
  Percent_resist_ISO_clade %>% 
  filter(pESI_presence == 'partial') %>%
  ggplot(aes(x=class, y=perc_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=50, label=round(perc_genomes, 1))) +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~ag_match_lump+ tot_n_genomes, nrow=1) +
  scale_y_continuous(limits = c(0,100),breaks = c(0,50,100))+  coord_flip()  +
  theme(axis.text.x=element_text(size=10, angle = -45, hjust = 0), 
        legend.position = 'none') + 
  ylab('Percent of Genomes') + 
  xlab('') + 
  theme_half_open()+
  background_grid()+
  ggtitle('AMR resistance for genomes with partial pESI plasmid')

P_inf_12

#######

###### THIS ONE
P_inf_13 <-
  Percent_resist_ISO_clade %>% 
  filter(pESI_presence == 'absent') %>%
  ggplot(aes(x=class, y=perc_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=50, label=round(perc_genomes, 1))) +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~ag_match_lump+ tot_n_genomes, nrow=1) +
  scale_y_continuous(limits = c(0,100),breaks = c(0,50,100))+  coord_flip()  +
  theme(axis.text.x=element_text(size=10, angle = -45, hjust = 0), 
        legend.position = 'none') + 
  ylab('Percent of Genomes') + 
  xlab('') + 
  theme_half_open()+
  background_grid()+
  ggtitle('AMR resistance for genomes with no pESI plasmid')

P_inf_13



infantis_plots <- 
list(P_inf_1, 
     P_inf_2,
     P_inf_3,
     P_inf_4, 
     P_inf_5, 
     P_inf_6, 
     P_inf_7, 
     P_inf_8, 
     P_inf_9, 
     P_inf_10, 
     P_inf_11, 
     P_inf_12, 
     P_inf_13)


infantis_plots %>% write_rds('output/infantis_plots.rds')



######


# 
# ### should run stats against PESI pres, 
# glm(family='binomial', pesi_prez ~ Year + ag_host + country, data=meta)
# 
# library(tidyverse)
# 
# ex_dat <- tibble(genome=paste0('genome_', seq(75)),
#                  year_category=sample(c('Before 2017', 'After 2017'),replace=T, size = 75),
#                  phage_presence=sample(c(0,1), size = 75, replace = T))
# 
# 
# 
# ex_dat
# 
# 
# ### Fisher exact test
# fish_test <-
#   fisher.test(
#     table(ex_dat$year_category, ex_dat$phage_presence)
#   )
# 
# fish_test
# 
# 
# ### logistic regression
# log_reg <- glm(data=ex_dat, formula=phage_presence ~ year_category, family = 'binomial')
# summary(log_reg)
# 
# broom::tidy(log_reg, conf.int=T, exponentiate=T)
# 
# 
# 
# ### fisher's exact gives very similar results to logistic regression
# # logistic regression odds ratio:
# broom::tidy(log_reg, conf.int=T, exponentiate=T)$estimate[2]
# 
# # fisher odds ratio:
# fish_test$estimate
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## EXPERIMENTAL ZONE  ###
# 
# ##### gene presence/absence variability within clusters
# # given a distance matrix and a set of isolates
# # calculate cumulative pairwise distances divided by the number of isolates
# #   mean_pw_dist = clusters with higher mean_pw_dist have isolates that have gained/lost genes in the transmission event?
# 
# 
# 
# 
# trans_meta <- meta %>% filter(!is.na(phy_clust))
# 
# 
# gpa <- read_tsv('all_pan/WRITE/gene_presence_absence.Rtab')
# 
# gpa <- gpa %>% column_to_rownames(var='Gene') %>% as.matrix()
# 
# gpa <- gpa[,colnames(gpa) %in% trans_meta$asm_acc]
# gpa <- gpa[rowSums(gpa) > 0,]
# 
# 
# 
# meta <- meta %>% filter(asm_acc %in% colnames(gpa))
# 
# trans_meta <- trans_meta %>% filter(asm_acc %in% colnames(gpa))
# # library(parallelDist)
# # library(usedist)
# 
# # trans_dists <- parallelDist(t(gpa), method = 'binary')
# 
# # write_rds(trans_dists,'output/transmission_cluster_dists.rds')
# 
# # attributes(trans_dists)$Labels
# 
# # calc_gene_variability <- function(dist_mat){
# #   # 1) subset dist to only contain the asm_accs
# #   # 2) move to long format
# #   # 3) sum(dists) / n()
# #   # browser()
# #   from_centroid <- dist_to_centroids(dist_mat, rep(1, length(attributes(dist_mat)$Labels)))
# #   return(from_centroid)
# #   
# #   
# # }
# 
# # sum(trans_dists) / attributes(trans_dists)$Size
# # attributes(trans_dists)$Size
# # length(trans_dists)
# 
# extract_PA <- function(asm_accs, pan_PA){
#   # browser()
#   pan_PA <- pan_PA[,colnames(pan_PA) %in% asm_accs]
#   pan_PA <- pan_PA[rowSums(pan_PA) > 0,]
#   return(pan_PA)
# }
# 
# 
# extract_variable_genes <- function(gpa){
#   # browser()
#   loggpa <- gpa > 0
#   res <- apply(X = loggpa, MARGIN=1, all)
#   variable_genes <- names(res[!res])
#   var_gene_PA <- gpa[rownames(gpa) %in% variable_genes,]
#   
# }
# 
# test <- 
#   trans_meta %>% 
#   group_by(phy_clust) %>%
#   nest() %>% 
#   left_join(phydel_sum) %>% 
#   filter(num_isolates > 10) %>% 
#   ungroup()
# 
# test1 <- 
#   test %>% 
#   mutate(clust_PA=map(.x=data, ~extract_PA(pan_PA = gpa, asm_accs = .x$asm_acc)), 
#          clust_dist=map(clust_PA, ~dist(.x, method = 'binary'))) 
# 
# 
# # test2 <- test1 %>% mutate(dists_to_centroids=map(clust_dist, calc_gene_variability))
# 
# test2 <- test1 %>%
#   mutate(var_genes=map(clust_PA, extract_variable_genes), 
#          num_var_genes=map_int(var_genes, ~dim(.x)[1])) %>% 
#   arrange(desc(num_var_genes))
# 
# test2 %>% 
#   mutate(num_years=map_int(.x=strsplit(years, '_'), .f = length)) %>% 
#   ggplot(aes(x=num_years, y=num_var_genes)) + geom_point() + 
#   geom_smooth(method = 'lm')
# 
# test2 %>% filter(phy_clust == 1248) 
# 
# test %>% filter(phy_clust == 919) %>% 
#   mutate(clust_PA=map(.x=data, ~extract_PA(pan_PA = gpa, asm_accs = .x$asm_acc)))
# 
# 
# 
# LOOK <- meta %>% filter(phy_clust == 919)
# usedist::dist_to_centroids()
# 
# 
# which.min(test$num_isolates)
# test %>% filter(num_isolates == 1)
# gpa
# .GlobalEnv$gpa
# trans_meta %>% filter(phydel)
# 
# dist_subset(trans_dists, trans_meta[1:10]$asm_acc)
# all(attributes(trans_dists)$Labels %in% trans_meta$asm_acc)
# 
# 
# all(trans_meta$asm_acc %in% attributes(trans_dists)$Labels)
# 
# trans_meta <- trans_meta %>% filter(asm_acc %in% attributes(trans_dists)$Labels)
# 
# m4 <- matrix(1:16, nrow=4, dimnames=list(LETTERS[1:4]))
# dm4 <- dist(m4)
# dist_subset(dm4, c("A", "B", "C"))
# dist_subset(dm4, c("D", "C", "B", "A"))
# 
# attributes(dm4)
# 
# TESTTEST <- c('A', 'B', 'C', 'D')
# dist_subset(dm4, TESTTEST)
# 
# 
# 
# 
# #######
# 
# 
# 
# rgps <- read_tsv('./all_pan/WRITE/plastic_regions.tsv')
# 
# modules <- read_tsv('./all_pan/WRITE/modules_summary.tsv')
# 

