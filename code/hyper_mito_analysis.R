library(tidyverse)
library(glue)
library(cowplot)
library(rstatix)
library(readxl)
library(openxlsx)
library(extrafont)
library(rio)
library(lubridate)
library(gdata)
library(MESS)
library(patchwork)
library(viridis)

loadfonts()

theme_set(theme_cowplot(14, font_family = "Arial"))


hyper4 = read_excel("supplementary_tables/Table S5.xlsx", 
           sheet = "H7mito_data")

well_bad = c("E4", "F2", "B6", "C4", "C9", "D2", "G3")

# filter out failed experiments
hyper4 = hyper4 %>%
  filter(!Well %in% well_bad) 

hyper4

treat_cols = c(
  "0mM" = "grey50",
  "10mM" = "#F5DF18",
  "20mM" = "#EB352A"
)

# summarise data

hyper4_sum = hyper4 %>% 
  group_by(Well, genotype, Treatment, Time) %>% 
  summarise(mean = mean(Ratio),
            sd = sd(Ratio))


hyper4_sum = hyper4_sum %>% 
  group_by(Treatment, genotype, Time) %>% 
  mutate(Replicate = 1:n()) 



hyper4_sum_condition = hyper4 %>% 
  group_by(genotype, Treatment, Time) %>% 
  summarise(mean = mean(Ratio),
            sd = sd(Ratio)) %>% 
  ungroup

hyper4_sum_condition %>% 
  # filter(genotype == "WT") %>% 
  ggplot(aes(x = Time, 
             y = mean, 
             color = Treatment,
             fill = Treatment)) +
  geom_point(data = hyper4, 
             aes(x = Time, 
                 y = Ratio, 
                 color = Treatment),
             alpha = 0.05,
             color = "black",
             size = 0.5) +
  geom_ribbon(aes(ymin = mean - sd, 
                  ymax = mean + sd), 
              alpha = 0.2) +
  geom_line() +
  scale_fill_manual(values = treat_cols) +
  facet_grid(genotype ~ Treatment) +
  scale_color_manual(values = treat_cols) +
  background_grid(major = "xy", minor = "none") +
  panel_border() +
  # breaks for x axis every 20
  scale_x_continuous(breaks = seq(0, 160, 20)) +
  ylim(-1, 6) +
  labs(x = "Time (min)", 
       y = "Ratio GFP/roGFP", 
       title = "Mean and SD of Ratio over Time") +
  theme(legend.position = "bottom")





## Track changes hyper 4 -----------------------------------------------------------

## treatment ---------------------------

# micit treatment

micit_auc = hyper4_sum %>% 
  filter(Time > 30 & Time < 91) %>% 
  group_by(Well, genotype, Treatment) %>% 
  summarise(AUC = auc(Time, mean))


micit_plot = micit_auc %>%
  ggplot(aes(x = Treatment, y = AUC, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~genotype) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_manual(values = treat_cols) +
  labs(x = "2-MiCit treatment", 
       y = "Cumulative ROS production") +
  theme(legend.position = "none")

micit_plot



# 2way anova stats

micit_anova = aov(micit_auc$AUC ~ micit_auc$Treatment * micit_auc$genotype) %>% 
  broom::tidy()

micit_ttest = micit_auc %>% 
  group_by(genotype) %>% 
  t_test(AUC ~ Treatment, p.adjust.method = "fdr", 
         ref.group = "0mM",
         detailed = TRUE)
