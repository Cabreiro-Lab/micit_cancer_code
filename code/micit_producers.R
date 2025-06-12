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
                     sheet = "cohort_production", skip = 1) %>% 
  mutate(Group = as.factor(Group),
         MicitFlux = as.numeric(MicitFlux))%>% 
  drop_na(MicitFlux)




cohorts %>% 
  ggplot(aes(x = Group, y = MicitFlux, fill = Group)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3),
             show.legend = F) +
  # facet_wrap(~cohort, scales = 'free_x', nrow = 1) +
  facet_grid(~cohort, scales = "free", space = "free") +
  scale_fill_manual(values = c(
    "Adenocarcinoma" = "#EB403D",
    "Adeno\r\ncarcinoma (III-IV)" = "#EB403D",
    "Adeno\r\ncarcinoma (I-II)" = "#EB403D",
    "Adeno\r\ncarcinoma (0)" = "#EB403D",
    "Adeno\r\ncarcinoma (young)" = "#EB403D",
    "Adeno\r\ncarcinoma (elderly)" = "#EB403D",
    "Adenoma" = "#D3D93B",
    "Healthy" = "#749DCC",
    "Control" = "#749DCC",
    "Healthy\r\n(Post surgery)" = "#749DCC",
    "Healthy (young)" = "#749DCC",
    "Healthy (elderly)" = "#749DCC"
  )) +
  labs(x = NULL,
       y = "Predicted 2MiCit production (AU)") +
  theme_cowplot(10, font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



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





