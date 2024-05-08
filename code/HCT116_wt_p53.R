library(tidyverse)
library(readr)
library(readxl)
library(cowplot)
library(rstatix)
library(extrafont)

theme_set(theme_cowplot(15))

# data loading ------------------------------------------------------------


wt_cells =  read_excel("data_for_R.xlsx", sheet = 'WT') %>% 
  pivot_longer(Control_1:Micit_FU_3, 
               names_to = 'Condition', 
               values_to = 'value') %>% 
  mutate(Replicate = str_sub(Condition, start = -1),
         Condition = str_sub(Condition, end = -3),
         Condition = factor(Condition, levels = c('Control','Micit',
                                                  'FU', 'Micit_FU'))) %>% 
  select(Elapsed, Condition, Replicate, value) %>% 
  mutate(Cell = 'WT', .before = Replicate)


p53_cells =  read_excel("data_for_R.xlsx", sheet = 'p53') %>% 
  pivot_longer(Control_1:Micit_FU_3, 
               names_to = 'Condition', 
               values_to = 'value') %>% 
  mutate(Replicate = str_sub(Condition, start = -1),
         Condition = str_sub(Condition, end = -3),
         Condition = factor(Condition, levels = c('Control','Micit',
                                                  'FU', 'Micit_FU'))) %>% 
  select(Elapsed, Condition, Replicate, value) %>% 
  mutate(Cell = 'p53', .before = Replicate)


cells = bind_rows(wt_cells, p53_cells) %>% 
  mutate(Cell = factor(Cell, levels = c('WT', 'p53')))


# plot growth curves ------------------------------------------------------


cells_sum = cells %>% 
  group_by(Elapsed, Condition, Cell) %>% 
  summarise(Mean = mean(value),
            SD = sd(value))


cells_sum %>% 
  ggplot(aes(x = Elapsed, y = Mean, color = Condition, fill = Condition)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth'
  ) +
  scale_fill_manual(values = c(
    '#3AD535', # green
    '#2065D6', # blue
    '#D69E2B', # orange
    '#D62055' # red
  )) +
  scale_color_manual(values = c(
    '#3AD535', # green
    '#2065D6', # blue
    '#D69E2B', # orange
    '#D62055' # red
  )) +
  facet_wrap(~Cell)

ggsave('plots/growth_curves_1000cells.pdf',
       height = 4.5, width = 9)


# PAPER VERSION -----------------------------------------------------------

cells_sum %>% 
  filter(!(str_detect(Condition, "FU"))) %>% 
  unite(sample, Cell, Condition, sep = ' ') %>% 
  ggplot(aes(x = Elapsed, y = Mean, 
             color = sample, 
             fill = sample)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.5,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth (Confluence %)'
  ) +
  scale_fill_manual(values = c(
    "WT Control" = '#3221E5', # green
    "p53 Control" = '#E52421', # blue
    "WT Micit" = '#2FC4E5', # orange
    "p53 Micit" = '#E5983D' # red
  )) +
  scale_color_manual(values = c(
    "WT Control" = '#3221E5', # green
    "p53 Control" = '#E52421', # blue
    "WT Micit" = '#2FC4E5', # orange
    "p53 Micit" = '#E5983D' # red
  )) +
  theme_cowplot(15, font_family = "Arial")

ggsave('plots/growth_curves_PAPER.pdf',
       height = 4.5, width = 5.5)




# calculate AUCs  ---------------------------------------------------------

library(MESS)


cells_auc = cells  %>% 
  group_by(Condition, Cell, Replicate) %>% 
  summarise(AUC = auc(Elapsed, value))


cells_auc %>% 
  ggplot(aes(x = Condition, y = AUC, fill = Condition)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~Cell) +
  scale_fill_manual(values = c(
    '#3AD535', # green
    '#2065D6', # blue
    '#D69E2B', # orange
    '#D62055' # red
  )) +
  labs(
    y = 'AUC (a.u.)'
  ) +
  ylim(3000, 10500) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  
ggsave('plots/AUC_boxplot.pdf',
       height = 5.5, width = 6)




# stats -------------------------------------------------------------------


library(openxlsx)

stats_by_cell = cells_auc %>% 
  group_by(Cell) %>% 
  t_test(AUC ~ Condition, detailed = T)

stats_by_condition = cells_auc %>% 
  group_by(Condition) %>% 
  t_test(AUC ~ Cell, detailed = T) %>% 
  add_significance('p')

  
interaction = cells_auc %>% 
  ungroup %>% 
  anova_test(AUC ~ Cell*Condition)





list_of_datasets = list(
  'stats by cell' = stats_by_cell,
  'stats by condition'=stats_by_condition,
  'two-way ANOVA' = interaction
)

write.xlsx(list_of_datasets, 'stats_HCT116_WT_p53.xlsx')


# pair by pair interactions
micit_control = cells_auc %>% 
  filter(Condition %in% c('Micit', 'Control')) %>% 
  ungroup %>% 
  anova_test(AUC ~ Cell*Condition)

fu_control = cells_auc %>% 
  filter(Condition %in% c('FU', 'Control')) %>% 
  ungroup %>% 
  anova_test(AUC ~ Cell*Condition)

micitFU_control = cells_auc %>% 
  filter(Condition %in% c('Micit_FU', 'Control')) %>% 
  ungroup %>% 
  anova_test(AUC ~ Cell*Condition)

micitFU_micit = cells_auc %>% 
  ungroup %>% 
  filter(Condition %in% c('Micit_FU', 'Micit')) %>% 
  anova_test(AUC ~ Cell*Condition)

micitFU_FU = cells_auc %>% 
  ungroup %>% 
  filter(Condition %in% c('Micit_FU', 'FU')) %>% 
  anova_test(AUC ~ Cell*Condition)

list_of_datasets = list(
  'Micit vs control' = micit_control,
  '5FU vs control'=fu_control,
  'Micit+FU vs control' = micitFU_control,
  'Micit+FU vs micit' = micitFU_micit, 
  'Micit+FU vs FU' = micitFU_FU
)

write.xlsx(list_of_datasets, 'interactions_HCT116_WT_p53.xlsx')




# growth rates ------------------------------------------------------------

library(growthrates)

# helper function
easyfit = function(df, h){
  fit = fit_easylinear(df$Elapsed, df$value, h = h)
  return(coef(fit))
}

test = cells %>% 
  filter(Condition == 'Control', Cell == 'WT', Replicate == 1)

easyfit(test, h = 10)
fit_names = names(easyfit(test, h = 10))


mumax = cells %>% 
  # remove first timepoints
  filter(Elapsed > 20) %>% 
  # filter(Cell == 'WT', Condition == 'Control') %>% 
  group_by(Condition, Cell, Replicate) %>%
  nest() %>%
  # modify the value of h to take a larger or smaller time window
  summarise(numax = map(data, easyfit, 30)) %>% 
  unnest_wider(numax) %>% 
  ungroup

mumax %>% 
  filter(!(str_detect(Condition, "FU"))) %>% 
  unite(sample, Cell, Condition, sep = ' ') %>% 
  mutate(sample = factor(sample, 
                         levels = c("WT Control",
                                    "p53 Control",
                                    "WT Micit",
                                    "p53 Micit"))) %>% 
  ggplot(aes(x = sample, y = mumax, fill = sample)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge(),
             show.legend = F) +
  scale_fill_manual(values = c(
    "WT Control" = '#3221E5', # green
    "p53 Control" = '#E52421', # blue
    "WT Micit" = '#2FC4E5', # orange
    "p53 Micit" = '#E5983D' # red
  )) +
  labs(
    x = NULL, 
    y = 'Max growth rate (Hour^-1)'
  ) +
  theme_cowplot(15, font_family = "Arial")

ggsave('plots/growth_rates_PAPER.pdf', height = 5, width = 5.5)

mumax %>% 
  filter(!(str_detect(Condition, "FU"))) %>% 
  unite(sample, Cell, Condition, sep = ' ') %>% 
  mutate(sample = factor(sample, 
                         levels = c("WT Control",
                                    "p53 Control",
                                    "WT Micit",
                                    "p53 Micit"))) %>% 
  write_csv("plots/growth_rates_PAPER.csv")


cells_sum %>% 
  filter(!(str_detect(Condition, "FU"))) %>% 
  unite(sample, Cell, Condition, sep = ' ') %>% 
  ggplot(aes(x = Elapsed, y = Mean, 
             color = sample, 
             fill = sample)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.5,
              color = NA) +
  geom_line() +
  labs(
    x = 'Time elapsed (in hours)',
    y = 'Growth (Confluence %)'
  ) +
  scale_fill_manual(values = c(
    "WT Control" = '#3221E5', # green
    "p53 Control" = '#E52421', # blue
    "WT Micit" = '#2FC4E5', # orange
    "p53 Micit" = '#E5983D' # red
  )) +
  scale_color_manual(values = c(
    "WT Control" = '#3221E5', # green
    "p53 Control" = '#E52421', # blue
    "WT Micit" = '#2FC4E5', # orange
    "p53 Micit" = '#E5983D' # red
  )) +
  theme_cowplot(15, font_family = "Arial")

ggsave('plots/growth_curves_PAPER.pdf',
       height = 4.5, width = 5.5)


cells_sum %>% 
  filter(!(str_detect(Condition, "FU"))) %>% 
  write_csv("plots/growth_curves_PAPER.csv")





numax %>% 
  group_by(Cell) %>% 
  rstatix::t_test(mumax ~ Condition, detailed = TRUE) %>% 
  adjust_pvalue(method = 'fdr') %>% 
  add_significance("p.adj")  %>% 
  write_csv('ttest_growth_rate.csv')


numax %>% 
  filter(!(str_detect(Condition, "FU"))) %>% 
  rstatix::anova_test(mumax ~ Condition * Cell, 
                      detailed = TRUE, 
                      effect.size = "pes") %>% 
  write_csv('interaction_growth_rate.csv')


many_splines = all_splines(value ~ Elapsed | Replicate + Cell + Condition, 
                           data = cells %>% 
                             filter(Elapsed > 30), optgrid = 30)



par(mfrow = c(8, 3))
par(mar = c(2.5, 4, 2, 1))
plot(many_splines)





# racemic mixture growth curves -------------------------------------------




library(readxl)
racemic_mixture = read_excel("racemic_mixture.xlsx")


mix.sum = racemic_mixture %>% 
  pivot_longer(Control_1:Batch_4, names_to = 'Condition', 
               values_to = 'Confluence') %>% 
  separate(Condition, into = c('Condition', 'Replicate'), sep = '_') %>% 
  mutate(Condition = case_when(Condition == 'Mixture' ~ 'Racemic mixture',
                               Condition == 'Batch' ~ '2R, 3S', 
                               TRUE ~ Condition),
         Condition = factor(Condition, 
                            levels = 
                              c('Control', 
                                'Racemic mixture',
                                '2R, 3S'))) %>% 
  group_by(Elapsed, Condition) %>% 
  summarise(mean_conf = mean(Confluence),
            sd_conf = sd(Confluence)) %>% 
  ungroup


mix.sum %>% 
  ggplot(aes(x = Elapsed, y = mean_conf, colour = Condition, fill = Condition)) +
  geom_ribbon(aes(ymin = mean_conf - sd_conf, ymax = mean_conf + sd_conf),
              alpha = 0.3, color = NA) +
  geom_line() +
  scale_fill_manual(values = c(
    '#196FE0',
    '#0AA652',
    '#F01853'
  )) +
  scale_color_manual(values = c(
    '#196FE0',
    '#0AA652',
    '#F01853'
  )) +
  labs(
    x = 'Time (h)',
    y = 'Confluency (%)'
  ) +
  theme_cowplot(15) +
  theme(
    legend.position = c(0.06, 0.8),
    legend.title = element_blank()
  )

ggsave('plots/racemic_mixture.pdf',
       height = 4.5, width = 6 
       )  
  
  
  
  

