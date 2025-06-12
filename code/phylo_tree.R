# tree plot
library(tidyverse)
library(tidytree)
library(treeio)
library(seqinr)
library(ggtree)
library(openxlsx)
library(treeio)
library(tidytree)
library(here)
library(readxl)
library(colorspace)
library(phytools)
library(ggtreeExtra)


### Load the phenotype data

## Bacteria
bac_scores = read_xlsx('supplementary_tables/Table S1.xlsx', 
                   sheet = 'screen_bact_growth', skip = 1) %>% 
  rename(Drug = `5FU`)

bac_sum = bac_scores %>% 
  group_by(Drug, Media, Bacteria) %>% 
  summarise(mean_AUC = mean(AUC)) %>% 
  ungroup


bac_sum_controls = bac_sum %>%
  filter(Drug == 0) %>% 
  rename(AUC_control = mean_AUC) %>% 
  ungroup %>% 
  select(-Drug)

bac_scores = bac_sum %>% 
  left_join(bac_sum_controls) %>% 
  mutate(AUC_norm = mean_AUC/AUC_control) %>%
  filter(Drug != 0) %>% 
  group_by(Media, Bacteria) %>% 
  summarise(AUC_norm_score = sum(AUC_norm)) %>% 
  ungroup



## Worm
worm_scores = read_xlsx('supplementary_tables/Table S1.xlsx', 
          sheet = 'screen_worm_scores', skip = 1) %>% 
  rename(Drug = `FU`)


worm_sum = worm_scores %>%
  group_by(Drug, Media, Bacteria) %>% 
  summarise(median_Score = median(Score)) %>% 
  ungroup

worm_scores = worm_sum %>% 
  group_by(Media, Bacteria) %>% 
  summarise(worm_score = sum(median_Score)) %>% 
  ungroup


# theme_set(theme_light())


# phylo tree --------------------------------------------------------------


tree = read.tree('other_data/core_tree_mod.treefile')


tree_metadata = read_csv('other_data/core_tree_metadata.csv')


new_tip_names = tree$tip.label %>% 
  as_tibble() %>% 
  rename(tree_tip = value) %>% 
  left_join(tree_metadata) %>% 
  pull(Bacteria)

# quick check that names are in the right order
# tibble(new_tip_names, tree$tip.label) %>% view

tree$tip.label = new_tip_names


# quick representation of the tree to see it's ok

ggtree(tree, layout = 'fan', branch.length = 'none') + 
  geom_tiplab() + theme(legend.position = 'none')


# complete tree

scores_auc_df = worm_scores %>%
  left_join(tree_metadata, by = c('Bacteria' = 'Bacteria')) %>% 
  select(Bacteria = strain_name, Media, worm_score) %>%
  # filter(Media == "NGM") %>% 
  as.data.frame


bac_scores_df = 
  bac_scores %>% 
  left_join(tree_metadata) %>% 
  select(Bacteria = strain_name , Media, bac_score = AUC_norm_score) %>% 
  as.data.frame

# MAKE A COPY 
tree_full = tree

tree_full = as_tibble(tree_full) %>% 
  full_join(tree_metadata %>% 
              rename(label = Bacteria)) %>% 
  select(-tree_tip) %>% 
  mutate(label = strain_name) %>% 
  select(-strain_name)

tree_full = as.phylo(tree_full)

tree_metadata_phyla = tree_metadata %>% 
  select(label = strain_name,
         phylum)

tree_full = full_join(tree_full , tree_metadata_phyla, by = 'label')


# base tree
p = ggtree(tree_full, layout = 'circular', 
           branch.length = 'none') + 
  geom_tiplab(aes(color = phylum),
              size = 2) +
  scale_color_manual(values =
                       c("Actinobacteria" = "#E0CB33",
                         "Firmicutes" = "darkred",
                         "Proteobacteria" = "#42A8E0"))
# rotate tips to make it look better
p 



# adding the worm scores
p = p + geom_fruit(data = scores_auc_df, 
                   geom = geom_bar, 
                   aes(x = worm_score, 
                       y = Bacteria,
                       fill = Media), 
                   width = 0.7,
                   stat = 'identity',
                   position = position_dodgex(hexpand = 69,
                                              vexpand = .5),
                   pwidth=.58,
                   orientation="y", 
                   axis.params=list(
                     axis       = "x",
                     text.size  = 1.8,
                     hjust      = 1,
                     vjust      = 0.5,
                     nbreak     = 3,
                   )
) +
  scale_fill_manual(values = c("#256EC2", "#F08622")) +
  labs(fill = "Worm development score")

# adding the bact scores
p + 
  ggnewscale::new_scale_fill() +
  geom_fruit(data = bac_scores_df, 
             geom = geom_bar, 
             aes(x = bac_score, 
                 y = Bacteria,
                 fill = Media), 
             width = 0.7,
             stat = 'identity',
             position = position_dodgex(hexpand = 90,
                                        vexpand = .5),
             pwidth=.8,
             orientation="y", 
             axis.params=list(
               axis       = "x",
               text.size  = 1.8,
               hjust      = 1,
               vjust      = 0.5,
               nbreak     = 3,
             )
  ) +
  scale_fill_manual(values = c("black", "grey")) +
  labs(fill = 'Bacterial growth score')

