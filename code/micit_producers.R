# script to generate the plots for the micit prediction part

library(tidyverse)
library(cowplot)
library(readxl)
library(glue)
library(broom)
library(ggrepel)
library(ggtern)
library(rstatix)
library(ggpubr)

theme_set(theme_cowplot(15))



# Cohorts -----------------------------------------------------------------


cohorts = read_excel("data/Table S4.xlsx", 
                     sheet = "cohort_production", skip = 1)


# stats
cohorts_stats = cohorts %>% 
  group_by(cohort) %>% 
  pairwise_wilcox_test(MicitFlux ~ Group, p.adjust.method = 'fdr') %>% 
  add_xy_position(x = "Group")



cohorts_stats$custom.label <- ifelse(cohorts_stats$p.adj <= 0.05, cohorts_stats$p.adj, "ns")


ggboxplot(cohorts, x = "Group", y = "MicitFlux", fill = "Group") +
  geom_jitter() +
  facet_wrap(~cohort, scales = 'free_x') +
  stat_pvalue_manual(cohorts_stats,label = "custom.label") 




# producers in cohorts ----------------------------------------------------

cohorts_p = read_excel("data/Table S4.xlsx", 
                     sheet = "cohort_sp_producers", skip = 1)

sp_select = cohorts_p %>% 
  group_by(cohort, Species) %>% 
  summarise(MicitProd_total = sum(MicitProd)) %>% 
  group_by(cohort) %>% 
  slice_max(MicitProd_total, n = 8) %>% 
  distinct(Species) %>% pull(Species)


cohorts_p %>% 
  filter(Species %in% sp_select) %>% 
  mutate(
    Species = str_replace(Species, "_", " "),
    Species = str_wrap(Species, 20)
  ) %>% 
  ggplot(aes(y = fct_reorder(Species, MicitProd, .desc = F), 
             x = MicitProd, fill = Species)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge(jitter.width = 0.4),
             size = 0.5,
             show.legend = F) +
  facet_wrap(~cohort) +
  theme_cowplot(14, font_family = "Arial") +
  labs(y = NULL, 
       x = "Predicted MiCit production") 





