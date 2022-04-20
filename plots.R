library(tidyverse)

infmeta <- read_tsv('metadata/Infantis_metadata.tsv')


infmeta %>% count(Year, PDS_acc) %>% arrange(desc(n))

infmeta %>% count(country, PDS_acc) %>% arrange(desc(n))

look <- 
  infmeta %>%
  filter(PDS_acc == 'PDS000073471.99') %>% 
  filter(ag_match != 'Human')

infmeta %>% 
  mutate(PDS_lump=fct_lump_n(f = PDS_acc, n = 8)) %>% 
  count(Year, PDS_lump) %>% 
  filter(Year < 2022) %>% 
  filter(Year > 2015) %>% 
  ggplot(aes(x=Year, y=n, color=PDS_lump)) +
  geom_line() 
