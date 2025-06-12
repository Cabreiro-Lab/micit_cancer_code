# script to generate the plots for the micit prediction part

library(tidyverse)
library(cowplot)
library(readxl)
library(glue)
library(broom)
library(ggrepel)
library(ggtern)
library(rstatix)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_cowplot(15))



# cell viability heatmap --------------------------------------------------

## read the data 
cell_viab = read_excel("supplementary_tables/Table S5.xlsx", 
                     sheet = "cell_viability")





# remove MDA-MB-453 because it didn't grow well
cell_metadata = read_xlsx("supplementary_tables/Table S5.xlsx", 
                          sheet = "cell_line_metadata") %>% 
  rename(cell = `Cell line`) %>% 
  mutate(cell = str_replace(cell, "CCD841 CoN", "CCD841_CoN")) # a point mutation in a name

cell_order = c("CCD841_CoN", "CCD-18Co", "HCT116", "SW1417", 
               "LoVo", "HT29", "SW837", "MCF7",
               "RKO", "T-47D", "SW48", "T84", "SK-CO-1", 
               "U2-OS", "Hs 578T", 
               "DLD-1",  "MDA-MB-231", "SW948", "LS411N")

cell_order %in% cell_metadata$cell

tissue_data = cell_metadata %>% 
  filter(cell %in% cell_order) %>% 
  select(cell, Tissue) %>% 
  column_to_rownames("cell")

tissue_data = tissue_data[cell_order,]



ha = HeatmapAnnotation(Tissue = tissue_data,
                       col = list(Tissue = c("Colon" = "#F0C73C", 
                                             "Breast" = "#EA43F0", 
                                             "Cecum" = "#3CF0CE",
                                             "Bone" = "#4B4BB5")))



# stats by biorep
auc_viab_stats = cell_viab %>% 
  mutate(Micit_mM = factor(Micit_mM)) %>% 
  group_by(cell) %>% 
  nest() %>% 
  mutate( model = map(data, ~aov(cell_proliferation~Micit_mM, data = .x)) ,
          tukey = map(.x = model, ~TukeyHSD(.x)),
          tukey = map(tukey, tidy))  %>% 
  select(-model,-data) %>% 
  unnest(cols = c(tukey)) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value)) 



auc_stats_matrix = auc_viab_stats %>% 
  filter(str_detect(contrast, "-0")) %>% 
  separate(contrast, into = c("Micit_mM", "control"), sep = "-") %>% 
  filter(cell %in% cell_order) %>% 
  select(cell, Micit_mM, p.stars) %>% 
  pivot_wider(names_from = Micit_mM, values_from = p.stars) %>% 
  mutate(`0` = "", .before = `1`) %>% 
  column_to_rownames("cell") %>% 
  as.matrix() 

auc_stats_matrix = t(auc_stats_matrix[cell_order,])

cell_matrix = cell_viab %>% 
  filter(cell %in% cell_order) %>% 
  group_by(cell, Micit_mM) %>% 
  summarise(mean_viab = mean(cell_proliferation)) %>% 
  ungroup %>% 
  pivot_wider(names_from = Micit_mM, values_from = mean_viab) %>% 
  column_to_rownames("cell") %>% 
  as.matrix()

cell_matrix = t(cell_matrix[cell_order,])

col_fun = colorRamp2(c(0.4, 1, 1.4), c("red", "white", "green"))

cell_matrix %>%
  Heatmap(cluster_columns = F,
          cluster_rows = F,
          col = col_fun,
          name = "Mean cell\nviability (%)", 
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%s", auc_stats_matrix[i, j]), 
                      x, y, 
                      gp = gpar(fontsize = 10))},
          top_annotation = ha,
          column_names_side = "top")




# cancer cell growth ------------------------------------------------------



cell_growth = read_excel("supplementary_tables/Table S5.xlsx", 
                       sheet = "cell_growth")


cell_growth %>% 
  mutate(Micit_mM = factor(Micit_mM, 
                           levels = c(0, 1,5,10)
                           )) %>% 
  filter(elapsed < 300) %>% 
  ggplot(aes(x = elapsed, y = mean_confluence, 
             color = Micit_mM, 
             fill = Micit_mM)) +
  geom_ribbon(aes(ymin = mean_confluence - SD, 
                  ymax = mean_confluence + SD), 
              color = NA, alpha = 0.2)+
  geom_line() +
  labs(y = 'Mean confluence (&plusmn; SD)',
       x = 'Time (h)') +
  scale_x_continuous(breaks = seq(0, 240, by = 60),
                     limits = c(0, 240)) +
  theme_cowplot(15,
                font_family = 'Arial') +
  panel_border() +
  theme(legend.position = "bottom", 
        strip.text = element_text(size = 13)) +
  scale_fill_viridis_d(option = "inferno", end = 0.9, direction = -1) + 
  scale_color_viridis_d(option = "inferno", end = 0.9, direction = -1) +
  facet_wrap(~cell, nrow = 4) +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = ggtext::element_markdown())


