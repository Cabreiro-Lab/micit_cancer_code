library(tidyverse)
library(readr)
library(readxl)
library(cowplot)
library(rstatix)
library(extrafont)
library(viridis)
library(randomcoloR)
library(MESS)

theme_set(theme_cowplot(15))



# data loading ------------------------------------------------------------


cells =  read_excel("raw_data/isomer_data_long.xlsx", sheet = 'growth') %>% 
  mutate(genotype = factor(genotype, levels = c("WT", "p53"))) %>% 
  mutate(micit_version = case_when(condition == 'DMSO' ~ 'DMSO', 
                                   TRUE ~ micit_version)) %>% 
  mutate(micit_version = factor(micit_version,
                                levels = c("none",
                                           "DMSO",
                                           "Sigma",
                                           "Newcastle",
                                           "isomer-1",
                                           "isomer-2",
                                           "isomer-3",
                                           "Target-2",
                                           "Target-3")),
         condition = factor(condition),
         micit_mM = factor(micit_mM, 
                           levels = c(0, 1, 5, 10)))








# remove bad data ---------------------------------------------------------


# replicate 2 from WT, exp_49, 1 mM and Sigma seems to be wrong, I'm deleting that

removals_1 = cells %>% filter(genotype == "WT", 
                 experiment == 'exp_49', 
                 micit_version == "Sigma", 
                 micit_mM == 1,
                 replicate == 2)

# replicate 2 from WT, exp_54, 10 mM and Target-2 seems to be wrong, I'm deleting that

removals_2 = cells %>% filter(genotype == "WT", 
                              experiment == 'exp_54', 
                              micit_version == "Target-2", 
                              micit_mM == 10,
                              replicate == 2)

# replicate 2 from WT, exp_54, 1 mM and Target-2 seems to be wrong, I'm deleting that

removals_3 = cells %>% filter(genotype == "WT", 
                              experiment == 'exp_54', 
                              micit_version == "Target-2", 
                              micit_mM == 1,
                              replicate == 2)

# replicate 2 from WT, exp_54, 1 mM and Target-2 seems to be wrong, I'm deleting that

removals_4 = cells %>% filter(genotype == "WT", 
                              experiment == 'exp_51', 
                              micit_version == "Sigma", 
                              micit_mM == 1,
                              replicate == 2)

# replicate 2 from WT, exp_54, 1 mM and Target-2 seems to be wrong, I'm deleting that

removals_5 = cells %>% filter(genotype == "p53", 
                              experiment == 'exp_51', 
                              micit_version == "Sigma", 
                              micit_mM == 1,
                              replicate == 2)

# replicate 2 from WT, exp_54, 1 mM and Target-2 seems to be wrong, I'm deleting that

removals_6 = cells %>% filter(genotype == "p53", 
                              experiment == 'exp_52', 
                              micit_version == "Target-3", 
                              micit_mM == 1,
                              biorep == 1,
                              replicate == 2)


removals_7 = cells %>% filter(genotype == "p53", 
                 experiment == 'exp_52', 
                 micit_version == "Target-3",
                 micit_mM == 10,
                 biorep == 2, 
                 replicate == 2 ) 


removals_8 = cells %>% filter(genotype == "WT", 
                              experiment == 'exp_54', 
                              micit_version == "DMSO",
                              micit_mM == 0,
                              biorep == 1, 
                              replicate == 1 ) 



cells = cells %>% 
  anti_join(removals_1) %>% 
  anti_join(removals_2) %>% 
  anti_join(removals_3) %>% 
  anti_join(removals_4) %>% 
  anti_join(removals_5) %>% 
  anti_join(removals_6) %>% 
  anti_join(removals_7) %>% 
  anti_join(removals_8)



# summarise data ----------------------------------------------------------

cells_sum = cells %>% 
  group_by(experiment, time, genotype, condition, micit_version, micit_mM, biorep) %>% 
  summarise(Mean = mean(confluence),
            SD = sd(confluence)) %>% 
  ungroup


# plot curves -------------------------------------------------------------

## experiment 49 ----------------



cells_sum %>% 
  filter(experiment == "exp_49") %>% 
  filter(genotype == 'WT') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(~micit_version)

ggsave('plots/growth_curves_exp49_WT.pdf',
       height = 4.5, width = 9)




cells_sum %>% 
  filter(experiment == "exp_49") %>% 
  filter(genotype == 'p53') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(~micit_version)

ggsave('plots/growth_curves_exp49_p53.pdf',
       height = 4.5, width = 9)




## experiment 54 ----------------


cells_sum %>% 
  filter(experiment == "exp_54") %>% 
  filter(genotype == 'WT') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(vars(micit_version, biorep))

ggsave('plots/growth_curves_exp54_WT.pdf',
       height = 4.5, width = 9)




cells_sum %>% 
  drop_na() %>% 
  filter(experiment == "exp_54") %>% 
  filter(genotype == 'p53') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(vars(micit_version, biorep))

ggsave('plots/growth_curves_exp54_p53.pdf',
       height = 4.5, width = 9)


## experiment 51 ----------------


cells_sum %>% 
  filter(experiment == "exp_51") %>% 
  filter(genotype == 'WT') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(vars(micit_version, biorep))

ggsave('plots/growth_curves_exp51_WT.pdf',
       height = 4.5, width = 9)




cells_sum %>% 
  drop_na() %>% 
  filter(experiment == "exp_51") %>% 
  filter(genotype == 'p53') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(vars(micit_version, biorep))

ggsave('plots/growth_curves_exp51_p53.pdf',
       height = 4.5, width = 9)



## experiment 52 ----------------


cells_sum %>% 
  filter(experiment == "exp_52") %>% 
  filter(genotype == 'WT') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(vars(micit_version, biorep))

ggsave('plots/growth_curves_exp52_WT.pdf',
       height = 4.5, width = 9)




cells_sum %>% 
  drop_na() %>% 
  filter(experiment == "exp_52") %>% 
  filter(genotype == 'p53') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_color_viridis_d(begin = .2, end = 0.9) +
  scale_fill_viridis_d(begin = .2, end = 0.9) +
  facet_wrap(vars(micit_version, biorep))

ggsave('plots/growth_curves_exp52_p53.pdf',
       height = 4.5, width = 9)

# calculate AUCs ----------------------------------------------------------




cells_auc = cells  %>% 
  group_by(experiment, genotype, condition, micit_version, micit_mM, biorep, replicate) %>% 
  summarise(AUC = auc(time, confluence)) %>% 
  ungroup



## exp 49 --------

cells_auc %>% 
  filter(experiment == 'exp_49') %>% 
  filter(genotype == 'WT') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp49_WT.pdf',
       height = 4.5, width = 8)

cells_auc %>% 
  filter(experiment == 'exp_49') %>% 
  filter(genotype == 'p53') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp49_p53.pdf',
       height = 4.5, width = 8)




## exp 51 --------

cells_auc %>% 
  filter(experiment == 'exp_51') %>% 
  filter(genotype == 'WT') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp51_WT.pdf',
       height = 4.5, width = 8)

cells_auc %>% 
  filter(experiment == 'exp_51') %>% 
  filter(genotype == 'p53') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp51_p53.pdf',
       height = 4.5, width = 8)



## exp 52 --------

cells_auc %>% 
  filter(experiment == 'exp_52') %>% 
  filter(genotype == 'WT') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp52_WT.pdf',
       height = 4.5, width = 8)

cells_auc %>% 
  filter(experiment == 'exp_52') %>% 
  filter(genotype == 'p53') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp52_p53.pdf',
       height = 4.5, width = 8)



## exp 54 --------

cells_auc %>% 
  filter(experiment == 'exp_54') %>% 
  filter(genotype == 'WT') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp54_WT.pdf',
       height = 4.5, width = 8)

cells_auc %>% 
  filter(experiment == 'exp_54') %>% 
  filter(genotype == 'p53') %>% 
  # unite(condition, condition, micit_version) %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(vars( genotype), scales = 'free')

ggsave('plots/AUC/exp54_p53.pdf',
       height = 4.5, width = 8)




# homogenise data ---------------------------------------------------------


# I need to filter times above the minimum time of experiment
min_time = cells %>% 
  group_by(experiment, genotype, replicate, biorep) %>% 
  slice_max(time) %>% 
  pull(time) %>% min


cells_good = cells %>% 
  filter(time <= min_time)


cells_good %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) %>% 
  write_csv("tables/cell_growth.csv")



cells_good_WT = cells_good %>% 
  filter(genotype == "WT") %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) 
  

cells_good_WT_control = cells_good_WT %>%
  filter(micit_version %in% c("Water", "Sigma", "Target-1")) %>%
  group_by(genotype, condition, time, micit_version, micit_mM, biorep) %>%
  summarise(confluence = mean(confluence, na.rm = T)) %>% 
  ungroup


cells_good_WT_control %>% 
  bind_rows(cells_good_WT %>% 
              filter(!(micit_version %in% c("Water", "Sigma", "Target-1")))) %>% 
  select(-biorep, -experiment) %>% 
  group_by(genotype, condition, micit_version, micit_mM, time) %>% 
  arrange(genotype, time, condition) %>% 
  mutate(replicate = 1:n()) %>%
  rename(Cell = genotype) %>% 
  mutate(Cell = "HCT116") %>% 
  write_csv("tables/cell_growth_WT.csv")


# cell good summary -------------------------------------------------------


## REMOVE EXP_49, NEWCASTLE, BOTH GENOTYPES FROM THE ANALYSIS, IS AN OUTLIER
cells_good_sum = cells_good %>%
  filter(!(micit_version == "Newcastle" & experiment == "exp_49")) %>%
  group_by(time, genotype, condition, micit_version, micit_mM) %>% 
  summarise(Mean = mean(confluence),
            SD = sd(confluence),
            SEM = SD / sqrt(n())) %>% 
  ungroup 


# cells_good_sum %>% 
#   write_csv("tables/cells_summary_stats.csv")

# 
# control_sigma_newc_means = cells_good_sum %>%
#   filter(micit_version %in% c("none", "Sigma", "Newcastle")) %>%
#   group_by(time, genotype, condition, micit_version, micit_mM) %>%
#   summarise(uMean = mean(Mean),
#             uSD = sd(Mean)) %>%
#   ungroup %>%
#   rename(Mean = uMean,
#          SD =  uSD)
# 
# 
# cells_good_sum = cells_good_sum %>%
#   filter(!(micit_version %in% c("none", "Sigma", "Newcastle"))) %>%
#   bind_rows(control_sigma_newc_means)






cols_growth = randomcoloR::randomColor(7)
cells_good_sum %>% 
  drop_na(Mean) %>% 
  filter(micit_mM %in% c(0,10)) %>%
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1"))) %>% 
  filter(genotype == 'p53') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_version, 
             fill = micit_version)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.5,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Confluence (%)'
  ) +
  scale_fill_manual(values = cols_growth) +
  scale_color_manual(values = cols_growth)

  

## PAPER VERSION ------

### WT =====
cells_good_sum %>% 
  drop_na(Mean) %>% 
  filter(micit_mM %in% c(0,10)) %>%
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) %>% 
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  mutate(micit_version = factor(micit_version, 
                                levels = c("Water", "DMSO", "Sigma", 
                                           "Target-1", "Target-2", "Target-3"))) %>% 
  filter(genotype == 'WT') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_version, 
             fill = micit_version)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.5,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Confluence (%) +- SD',
    color = "Isomer",
    fill = "Isomer"
  ) +
  scale_fill_manual(values = c(
    "Water" = "grey70",
    "DMSO" = "grey20",
    "Sigma" = "#915B4A",
    "Target-2" = "#296AE6",
    "Target-1" = "#75C82C",
    "Target-3" = "#E65729"
  )) +
  scale_color_manual(values = c(
    "Water" = "grey70",
    "DMSO" = "grey20",
    "Sigma" = "#915B4A",
    "Target-2" = "#296AE6",
    "Target-1" = "#75C82C",
    "Target-3" = "#E65729"
  ))


ggsave("plots/paper_figures/WT_growth_SEM.pdf", 
       height = 6, width = 10)


### p53 =====
cells_good_sum %>% 
  drop_na(Mean) %>% 
  filter(micit_mM %in% c(0,10)) %>%
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) %>% 
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  mutate(micit_version = factor(micit_version, 
                                levels = c("Water", "DMSO", "Sigma", 
                                           "Target-1", "Target-2", "Target-3"))) %>% 
  filter(genotype == 'p53') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_version, 
             fill = micit_version)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.5,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Confluence (%) +- SD',
    color = "Isomer",
    fill = "Isomer"
  ) +
  scale_fill_manual(values = c(
    "Water" = "grey70",
    "DMSO" = "grey20",
    "Sigma" = "#915B4A",
    "Target-2" = "#296AE6",
    "Target-1" = "#75C82C",
    "Target-3" = "#E65729"
  )) +
  scale_color_manual(values = c(
    "Water" = "grey70",
    "DMSO" = "grey20",
    "Sigma" = "#915B4A",
    "Target-2" = "#296AE6",
    "Target-1" = "#75C82C",
    "Target-3" = "#E65729"
  ))


ggsave("plots/paper_figures/p53_growth_SEM.pdf", 
       height = 6, width = 10)





### isomer growth -----------------------------------------------------------

cells_good_sum %>% 
  drop_na(Mean) %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) %>% 
  filter((micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2",
                               "Water"))) %>% 
  filter(genotype == 'WT',
         time > 48) %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_mM, 
             fill = micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.5,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Confluence (%) +- SD',
    color = "Isomer",
    fill = "Isomer"
  )  +
  facet_wrap(~micit_version)


### DMSO vs control plot ------------------------------

cells_good_sum %>% 
  drop_na(Mean) %>% 
  filter(micit_version %in% c("none", "DMSO")) %>% 
  # filter(genotype == 'p53') %>% 
  ggplot(aes(x = time, y = Mean, 
             color = micit_version, 
             fill = micit_version)) +
  geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM), alpha = 0.5,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Confluence (%)'
  ) +
  facet_wrap(~genotype, scale = 'free_y')

ggsave("plots/DMSO_control.pdf", 
       height = 6, width = 10)





cells_good %>% 
  filter(micit_version == "DMSO") %>% 
  unite(sample, experiment, condition, replicate, biorep, micit_mM) %>% 
  ggplot(aes(x = time, y = confluence, linetype = sample)) +
  geom_line() +
  facet_wrap(~genotype )



# AUC good ----------------------------------------------------------------



## REMOVE EXP_49, NEWCASTLE, BOTH GENOTYPES FROM THE ANALYSIS, IS AN OUTLIER

cells_auc = cells_good  %>% 
  filter(!(micit_version == "Newcastle" & experiment == "exp_49")) %>%
  group_by(experiment, genotype, condition, micit_version, micit_mM, biorep, replicate) %>% 
  summarise(AUC = auc(time, confluence)) %>% 
  ungroup

# calculate the AUC per biorep in the control, Sigma and Newcastle version
cells_auc_controls = cells_auc %>% 
  filter(micit_version %in% c("none", "Sigma", "Newcastle")) %>%
  group_by(experiment, genotype, condition, micit_version, micit_mM) %>% 
  summarise(Mean = mean(AUC)) %>% 
  rename(AUC = Mean) %>% 
  ungroup 


  
cells_auc = cells_auc %>% 
  filter(!(micit_version %in% c("none", "Sigma", "Newcastle"))) %>% 
  bind_rows(cells_auc_controls) %>% 
  select(-experiment, 
         -biorep)



cells_auc %>%
  group_by(genotype, condition, micit_version, micit_mM) %>% 
  mutate(replicate = 1:n()) %>% 
  write_csv("tables/cell_AUC_data.csv")

cells_auc %>%
  group_by(genotype, condition, micit_version, micit_mM) %>% 
  mutate(replicate = 1:n()) %>%
  mutate(genotype = "HCT116") %>% 
  rename(Cell = genotype) %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) %>% 
  write_csv("tables/cell_AUC_WT_data.csv")




# stats for AUC -----------------------------------------------------------


cells_auc %>%
  group_by(genotype, condition, micit_version, micit_mM) %>% 
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


### special case for isomers ---------

isomers_auc = cells_good  %>% 
  filter(micit_version %in% c("isomer-1", "isomer-2", "isomer-3"),
         time > 48) %>%
  group_by(experiment, genotype, condition, micit_version, micit_mM, biorep, replicate) %>% 
  summarise(AUC = auc(time, confluence)) %>% 
  ungroup

  
isomers_auc = cells_auc_controls %>% 
  group_by(genotype, condition, micit_version, micit_mM) %>% 
  mutate(replicate = 1:n()) %>% 
  bind_rows(isomers_auc %>% select(-biorep))
  





## controls plot --------

cells_auc %>% 
  filter(condition == "Control") %>% 
  ggplot(aes(x = genotype, y = AUC, fill = genotype)) +
  geom_boxplot() +
  geom_jitter() 

cells_auc %>% 
  drop_na(AUC) %>% 
  filter(micit_mM %in% c(0,10)) %>%
  # filter(!(micit_version %in% c("isomer-3",
  #                               "isomer-1"))) %>% 
  ggplot(aes(x = micit_version, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_jitter()+
  facet_wrap(~genotype) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("plots/AUC/AUC_boxplot_0vs10.pdf", 
       height = 6, width = 10)


cells_auc %>% 
  drop_na(AUC) %>% 
  filter(micit_mM %in% c(0,5)) %>%
  # filter(!(micit_version %in% c("isomer-3",
  #                               "isomer-1"))) %>% 
  ggplot(aes(x = micit_version, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_jitter()+
  facet_wrap(~genotype) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("plots/AUC/AUC_boxplot_0vs5.pdf", 
       height = 6, width = 10)




## BOXPLOTS PAPER ---------

cells_auc %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) %>% 
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  mutate(micit_version = factor(micit_version, 
                                levels = c("Water", "DMSO", "Sigma", 
                                           "Target-1", "Target-2", "Target-3"))) %>% 
  filter(genotype == 'WT') %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  facet_wrap(~genotype) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c(
    "Water" = "grey70",
    "DMSO" = "grey20",
    "Sigma" = "#915B4A",
    "Target-2" = "#296AE6",
    "Target-1" = "#75C82C",
    "Target-3" = "#E65729"
  )) +
  scale_y_continuous(limits = c(0, 15000), breaks = seq(0, 15000, 2500),
                     expand = expansion(mult = c(0., 0.02))) +
  labs(x = "2-MiCit (mM)",
       y = "AUC (a.u.)",
       fill = "Isomer") 

ggsave("plots/paper_figures/AUC_WT.pdf", 
       height = 6, width = 10)





## AUC stats paper ----------------------------------------------------------------



cells_auc %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "(2R,3S) Micit"),
         micit_version = str_replace(micit_version, "Target-2", "TetraMicit"),
         micit_version = str_replace(micit_version, "Target-3", "Inactive TetraMicit"),
         micit_version = str_replace(micit_version, "Sigma", "(2R,3S/2S,3R) Micit")) %>% 
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  filter(genotype == 'WT') %>% 
  unite(contrast, micit_version, micit_mM, remove = F) %>% 
  t_test(AUC ~ contrast, detailed = TRUE, p.adjust.method = 'fdr') %>% 
  rename(effect_size = estimate) %>% 
  select(-`.y.`) %>% 
  write_csv("tables/AUC_pairwise_stats.csv")

cells_auc %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "(2R,3S) Micit"),
         micit_version = str_replace(micit_version, "Target-2", "TetraMicit"),
         micit_version = str_replace(micit_version, "Target-3", "Inactive TetraMicit"),
         micit_version = str_replace(micit_version, "Sigma", "(2R,3S/2S,3R) Micit")) %>% 
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  filter(genotype == 'WT') %>% 
  unite(contrast, micit_version, micit_mM, remove = F) %>% 
  t_test(AUC ~ contrast, detailed = TRUE, p.adjust.method = 'fdr', 
         ref.group = "Water_0") %>% 
  rename(effect_size = estimate) %>% 
  mutate(effect_size = effect_size * -1) %>% 
  select(-`.y.`) %>% 
  write_csv("tables/AUC_pairwise_stats_Water.csv")



## isomers plot ------------------------------------------------------------

isomers_auc %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water")) %>% 
  filter((micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2",
                               "Water"))) %>% 
  mutate(micit_version = factor(micit_version, 
                                levels = c("Water", "isomer-1", "isomer-2", "isomer-3"))) %>% 
  filter(genotype == 'WT') %>% 
  ggplot(aes(x = micit_mM, y = AUC, fill = micit_version)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,
                                             jitter.width = 0.1)) +
  # facet_wrap(~genotype) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c(
    "Water" = "grey70",
    "isomer-1" = "#81A64C",
    "isomer-2" = "#336DA8",
    "isomer-3" = "#A856A1"
  )) +
  scale_y_continuous(limits = c(0, 13000), breaks = seq(0, 15000, 2500),
                     expand = expansion(mult = c(0., 0.02))) +
  labs(x = "2-MiCit (mM)",
       y = "AUC (a.u.)",
       fill = "Isomer") +
  theme_cowplot(font_family = "Arial")

ggsave("plots/paper_figures/AUC_isomers_WT.pdf", 
       height = 6, width = 10)





# AUC stats ---------------------------------------------------------------







cells_auc %>%
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "Target-1")) %>% 
  unite(condition, micit_version, micit_mM) %>% 
  group_by(genotype) %>% 
  pairwise_t_test(AUC ~ condition, detailed = T,
                  p.adjust.method = "fdr") %>% 
  mutate(group1 = str_replace(group1, "_", " "),
         group1 = str_c(group1, "mM"),
         group2 = str_replace(group2, "_", " "),
         group2 = str_c(group2, "mM")) %>% 
  arrange(genotype, group1, group2) %>%
  write_csv("tables/AUC_pairwise_stats.csv")





# max growth rate ---------------------------------------------------------



library(growthrates)
# helper function
easyfit = function(df, h){
  fit = fit_easylinear(df$time, df$confluence, h = h)
  return(coef(fit))
}



test = cells_good %>% 
  filter(micit_version == 'isomer-2', time > 10, 
         micit_mM == 10, biorep ==1, replicate ==1,
         genotype == "WT")

fit_easylinear(test$time, test$confluence, h = 30)

easyfit(test, h = 10)
fit_names = names(easyfit(test, h = 10))



cells_mumax = cells_good  %>% 
  filter(!(micit_version == "Newcastle" & experiment == "exp_49")) %>%
  group_by(experiment, genotype, condition, micit_version, micit_mM, biorep, replicate) %>% 
  filter(time > 30) %>% 
  filter(genotype == 'WT') %>%
  nest() %>%
  summarise(numax = map(data, easyfit, 30)) %>% 
  unnest_wider(numax) %>% 
  ungroup

# calculate the AUC per biorep in the control, Sigma and Newcastle version
cells_mumax_controls = cells_mumax %>% 
  filter(micit_version %in% c("none", "Sigma", "Newcastle")) %>%
  group_by(experiment, genotype, condition, micit_version, micit_mM) %>% 
  summarise(mumax = mean(mumax)) %>% 
  ungroup 



cells_mumax = cells_mumax %>% 
  filter(!(micit_version %in% c("none", "Sigma", "Newcastle"))) %>% 
  bind_rows(cells_mumax_controls) %>% 
  select(-experiment, 
         -biorep)



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
  filter(genotype == 'WT') %>% 
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

ggsave("plots/paper_figures/mumax_WT.pdf",
       height = 6, width = 10)




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
         p.adjust.method = 'BH') %>% 
  write_csv("tables/mumax_stats.csv")

cells_mumax %>% 
  mutate(micit_version = str_replace(micit_version, "none", "Water"),
         micit_version = str_replace(micit_version, "Newcastle", "(2R,3S)2-MiCit"),
         micit_version = str_replace(micit_version, "Target-2", "Trymethyl-2-MiCit"),
         micit_version = str_replace(micit_version, "Target-3", "Synthesis precursor"),
         micit_version = str_replace(micit_version, "Sigma", "(2R,3S/2S,3R)2-MiCit")) %>%  
  filter(!(micit_version %in% c("isomer-3",
                                "isomer-1",
                                "isomer-2"))) %>% 
  arrange() %>% 
  write_csv("tables/mumax_data.csv")
  