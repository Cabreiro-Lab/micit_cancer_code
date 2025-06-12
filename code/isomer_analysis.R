library(tidyverse)
library(readr)
library(readxl)
library(cowplot)
library(rstatix)
library(extrafont)
library(viridis)
library(randomcoloR)
library(MESS)
library(growthrates)

theme_set(theme_cowplot(15))


## read the data 
isomers = read_excel("supplementary_tables/Table S5.xlsx", 
                       sheet = "isomers_cell_growth")

glimpse(isomers)

cells_auc = isomers  %>% 
  group_by(condition, micit_version, micit_mM, replicate) %>% 
  summarise(AUC = auc(time, confluence)) %>% 
  ungroup



cells_auc %>%
  group_by(condition, micit_version, micit_mM) %>% 
  mutate(replicate = 1:n()) %>%
  mutate(genotype = "HCT116") %>% 
  rename(Cell = genotype) %>% 
  filter(!(str_detect(micit_version, "isomer"))) %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "(2R,3S) Micit"),
         micit_version = str_replace(micit_version, "Target-2", "TetraMicit"),
         micit_version = str_replace(micit_version, "Target-1", "(2R,3S/2S,3R) Micit"),
         micit_version = str_replace(micit_version, "Target-3", "Inactive TetraMicit")) %>% 
  view


cells_auc %>% 
  distinct(micit_version)






# helper function
easyfit = function(df, h){
  fit = fit_easylinear(df$time, df$confluence, h = h)
  return(coef(fit))
}


cells_mumax = isomers  %>% 
  mutate(micit_version = str_replace(micit_version, "Target-1", "Newcastle")) %>% 
  group_by(condition, micit_version, micit_mM, replicate) %>% 
  filter(time > 30) %>% 
  nest() %>%
  summarise(numax = map(data, easyfit, 30)) %>% 
  unnest_wider(numax) %>% 
  ungroup

# # calculate the AUC per biorep in the control, Sigma and Newcastle version
# cells_mumax_controls = cells_mumax %>% 
#   filter(micit_version %in% c("Water", "Sigma", "Newcastle")) %>%
#   group_by(condition, micit_version, micit_mM) %>% 
#   summarise(mumax = mean(mumax)) %>% 
#   ungroup 
# 
# 
# 
# cells_mumax = cells_mumax %>% 
#   filter(!(micit_version %in% c("none", "Sigma", "Newcastle"))) %>% 
#   bind_rows(cells_mumax_controls) 

cells_mumax %>% 
  distinct(micit_version)

cells_mumax %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "(2R,3S)2-MiCit"),
         micit_version = str_replace(micit_version, "Target-2", "Trymethyl-2-MiCit"),
         micit_version = str_replace(micit_version, "Target-3", "Synthesis precursor"),
         micit_version = str_replace(micit_version, "Sigma", "(2R,3S/2S,3R)2-MiCit")) %>%  
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  mutate(micit_version = factor(micit_version, 
                                levels = c("Water", "DMSO", "(2R,3S/2S,3R)2-MiCit", 
                                           "(2R,3S)2-MiCit", "Trymethyl-2-MiCit", 
                                           "Synthesis precursor"))) %>% 
  drop_na(micit_version) %>%
  mutate(micit_mM = factor(micit_mM)) %>%
  ggplot(aes(x = micit_mM, y = mumax, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  # facet_wrap(~genotype) +
  theme_cowplot(font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c(
    "Water" = "grey70",
    "DMSO" = "grey20",
    "(2R,3S/2S,3R)2-MiCit" = "#925C4B",
    "(2R,3S)2-MiCit" = "#79C143",
    "Trymethyl-2-MiCit" = "#BE1E2D",
    "Synthesis precursor" = "#F9ED32"
  )) +
  # scale_y_continuous(
  #                    expand = expansion(mult = c(0., 0.05))) +
  labs(x = "2-MiCit (mM)",
       y = expression(paste(mu^{max}, " (", time^{-1}, ")")),
       fill = "Isomer") 



cells_mumax %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "(2R,3S)2-MiCit"),
         micit_version = str_replace(micit_version, "Target-2", "Trymethyl-2-MiCit"),
         micit_version = str_replace(micit_version, "Target-3", "Synthesis precursor"),
         micit_version = str_replace(micit_version, "Sigma", "(2R,3S/2S,3R)2-MiCit")) %>%  
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  unite(sample, micit_version, micit_mM) %>% 
  t_test(mumax ~ sample, ref.group = "Water_0", detailed=T,
         p.adjust.method = 'BH')

