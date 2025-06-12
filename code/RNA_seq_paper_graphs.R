## RNAseq enrichment + PCA


### libraries ####
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(openxlsx)
library(viridis)
library(glue)
library(cowplot)
library(ggtext)
library(extrafont)
library(readxl)

extrafont::loadfonts()
theme_set(theme_cowplot(14, font_family = "Arial"))

rna_enrich = read_excel("supplementary_tables/Table S6.xlsx", 
           sheet = "RNA_seq KEGG", skip = 1)


selected_kegg_enrich = rna_enrich %>% 
  group_by(description) %>% 
  count() %>% 
  filter(n >= 3) %>% 
  pull(description)


rna_enrich %>%
  filter((description %in% selected_kegg_enrich)) %>%
  mutate(direction=factor(direction, levels = c('DOWN', 'UP'))) %>%
  mutate(FDR = cut(fdr, 
                   labels = c('0.001', '0.01', '0.05', 'ns'),
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf)),
         .before = term) %>%
  mutate(description = str_wrap(description, width = 25)) %>%
  ggplot(aes(y = description, x = direction ,fill = FDR)) +
  geom_tile() +
  scale_fill_manual(values = c('#2432FF',
                               '#616BFF',
                               '#A3A9FF',
                               '#FFFFFF')) +
  labs(
    x = NULL,
    y = NULL,
    caption = "All categories present in 3 or more cell lines"
  ) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~cell, nrow = 1) +
  panel_border()





# volcano plots -----------------------------------------------------------


rna_stats = read_excel("supplementary_tables/Table S6.xlsx", 
           sheet = "RNA_seq_stats", skip = 0)


rna_stats


rna_stats %>% 
  drop_na(padj) %>% 
  # filter(Contrast == "HT29") %>%
  mutate(Gene = case_when(padj < 0.05 ~ "Significant",
                          TRUE ~ "Not significant"),
         labels = case_when(Gene == "Significant" ~ gene_name)) %>% 
  filter(log2FoldChange < 3) %>% # filter an outlier
  mutate(Contrast = factor(Contrast, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                      'HT29', 'SK-CO-1', 'SW1417', 
                                      'SW948'))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = Gene)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not significant" = "black"
    )
  ) +
  xlim(-4, 3) +
  ylim(0, 180) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 10,
  #                          size = 3) +
  facet_wrap(~Contrast, nrow = 2) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(14, font_family = "Arial") +
  panel_border() +
  # legend at bottom
  theme(legend.position = "top")



