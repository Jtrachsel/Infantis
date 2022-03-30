# fix trees set up for phydelity

# I can't get gurobi to work very well on SCINet
# gurobi needs a license and this is linked to a 'Host_ID'
# it seems like this Host_ID changes for each compute node...
# so it will work the first time you set it up 
# but the next time you try and call gurobi, you're probably on a different 'host'
# I set up phydelity and gurobi on fsep11 instead of fighting with it



library(tidyverse)
library(ape)

fix_trees <- function(SNP_tree_path, dest_path){
  
  tr <- read.tree(SNP_tree_path)
  write.tree(tr, dest_path)
  return(TRUE)
  
}

fixed_for_phydelity <- 
  tibble(OG_path=list.files(path = '/90daydata/fsepru113/jtrachsel/Infantis/SNP_trees', pattern = 'newick', full.names = TRUE), 
         new_path=paste0(sub('(.*).newick_tree.newick','\\1',OG_path), '_fix.newick')) %>% 
  mutate(fixed=map2_lgl(.x = OG_path, .y=new_path, .f = ~fix_trees(SNP_tree_path = .x, dest_path = .y)))



list.files(path='/90daydata/fsepru113/jtrachsel/Infantis/SNP_trees/', pattern = 'fix')



       