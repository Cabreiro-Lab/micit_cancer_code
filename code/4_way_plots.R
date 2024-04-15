# script to generate the plots for the 4-way screen

library(tidyverse)
library(cowplot)
library(here)
library(rstatix)
library(readxl)
library(glue)
library(broom)
library(ggrepel)

theme_set(theme_cowplot(15))

# read data from Table S2 -------------------------------------------------

# data from BW
BW = read_excel("data/Table S2.xlsx", sheet = "4way_BW") %>% 
  mutate(worm_score_5 = as.numeric(worm_score_5),
         worm_score_1.5 = as.numeric(worm_score_1.5))

# data from pyrE
pyrE = read_excel("data/Table S2.xlsx", sheet = "4way_pyrE")

# data from the triple mutant
TM = read_excel("data/Table S2.xlsx", sheet = "4way_TM")

# model adjustments
adj = read_excel("data/Table S2.xlsx", sheet = "4way_adjustments")







# scatterplots ------------------------------------------------------------

worm_colors = c('red','orange', 'yellow','#71B83B') # worm phenotype

## BW case ----------------------------------------------------------------

# to plot different drug concentrations just change worm_score_5 for any
# of the other worm phenotypes at different drug concentrations

BW %>% 
  drop_na(worm_score_5) %>% 
  ggplot(aes(x = log2FC_control, 
             y = log2FC_5FU)) +
  geom_abline(intercept = 0,
              slope = 1, alpha = 1,
              color = 'grey', linetype = 'longdash') +
  geom_abline(aes(intercept = adj$intercept[1], slope = adj$slope[1]),
              alpha = 0.8, color = 'red') +
  geom_vline(aes(xintercept = 0), 
             alpha = 0.9, color = 'grey') +
  geom_hline(aes(yintercept = 0), 
             alpha = 0.9, color = 'grey') +
  geom_point(aes(fill = worm_score_5),
             shape = 21,
             size = 3) +
  scale_fill_gradientn(colours = worm_colors,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  labs(x = "E. coli growth (+metabolite/control) logFC\n No Drug",
       y = "E. coli growth (+metabolite/control) logFC\n 100 uM 5-FU")
  

## pyrE case ----------------------------------------------------------------


# to plot different drug concentrations just change worm_score_5 for any
# of the other worm phenotypes at different drug concentrations

pyrE %>% 
  drop_na(worm_score_5) %>% 
  ggplot(aes(x = log2FC_control, 
             y = log2FC_5FU)) +
  geom_abline(intercept = 0,
              slope = 1, alpha = 1,
              color = 'grey', linetype = 'longdash') +
  geom_abline(aes(intercept = adj$intercept[2], slope = adj$slope[2]),
              alpha = 0.8, color = 'red') +
  geom_vline(aes(xintercept = 0), 
             alpha = 0.9, color = 'grey') +
  geom_hline(aes(yintercept = 0), 
             alpha = 0.9, color = 'grey') +
  geom_point(aes(fill = worm_score_5),
             shape = 21,
             size = 3) +
  scale_fill_gradientn(colours = worm_colors,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  labs(x = "E. coli growth (+metabolite/control) logFC\n No Drug",
       y = "E. coli growth (+metabolite/control) logFC\n 100 uM 5-FU")

## TM case ----------------------------------------------------------------

# to plot different drug concentrations just change worm_score_250 for any
# of the other worm phenotypes at different drug concentrations

TM %>% 
  drop_na(worm_score_250) %>% 
  ggplot(aes(x = log2FC_control, 
             y = log2FC_5FU)) +
  geom_abline(intercept = 0,
              slope = 1, alpha = 1,
              color = 'grey', linetype = 'longdash') +
  geom_abline(aes(intercept = adj$intercept[3], slope = adj$slope[3]),
              alpha = 0.8, color = 'red') +
  geom_vline(aes(xintercept = 0), 
             alpha = 0.9, color = 'grey') +
  geom_hline(aes(yintercept = 0), 
             alpha = 0.9, color = 'grey') +
  geom_point(aes(fill = worm_score_250),
             shape = 21,
             size = 3) +
  scale_fill_gradientn(colours = worm_colors,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  labs(x = "E. coli growth (+metabolite/control) logFC\n No Drug",
       y = "E. coli growth (+metabolite/control) logFC\n 100 uM 5-FU")




# main scatterplots -------------------------------------------------------

# specify labels to plot 
labels_4way = c("Uracil_N","Cytidine_N","Uridine_N","Cytidine-3'-monophosphate",
                "Cytidine- 2',3'-cyclic monophosphate","Uridine-3'-monophosphate",
                "Uridine-2',3'-cyclic-monophosphate","Cytidine-2'-monophosphate",
                "Cytidine-5'-monophosphate","Uridine-2'-monophosphate",
                "Uridine-5'-monophosphate","Uridine","D-Galactose",
                "Glycerol","D-Sorbitol","D-Trehalose","Dulcitol",
                "Maltose","D-Mannose",'alpha-D-Glucose')

labels_4way_nogluc = c("Uracil_N","Cytidine_N","Uridine_N","Cytidine-3'-monophosphate",
                       "Cytidine- 2',3'-cyclic monophosphate","Uridine-3'-monophosphate",
                       "Uridine-2',3'-cyclic-monophosphate","Cytidine-2'-monophosphate",
                       "Cytidine-5'-monophosphate","Uridine-2'-monophosphate",
                       "Uridine-5'-monophosphate","Uridine","D-Galactose",
                       "Glycerol","D-Sorbitol","D-Trehalose","Dulcitol",
                       "Maltose","D-Mannose")

# specify nucleotide labels
labels_nuc = c("Uracil_N","Cytidine_N","Uridine_N","Cytidine-3'-monophosphate",
               "Cytidine- 2',3'-cyclic monophosphate","Uridine-3'-monophosphate",
               "Uridine-2',3'-cyclic-monophosphate","Cytidine-2'-monophosphate",
               "Cytidine-5'-monophosphate","Uridine-2'-monophosphate",
               "Uridine-5'-monophosphate","Uridine")
# specify sugar labels
labels_sug = c("D-Galactose","Glycerol","D-Sorbitol",
               "D-Trehalose","Dulcitol","Maltose","D-Mannose")

## BW case ----------------------------------------------------------------


BW %>% 
  drop_na(worm_score_5) %>% 
  mutate(score = abs(log2FC_5FU - (log2FC_control - adj$slope[1]))) %>% 
  arrange(score) %>% 
  mutate(new_labels = case_when(MetaboliteU %in% labels_4way_nogluc ~ MetaboliteU,
                                TRUE ~ ''),
         points_alpha = case_when(MetaboliteU %in% labels_4way ~ 1,
                                  TRUE ~ 0.15),
         label_size = case_when(worm_score_5 == 4 ~ 3,
                                TRUE ~ 2),
         label_color = case_when(MetaboliteU %in% labels_nuc ~ 'Nucleotide',
                                 MetaboliteU %in%  labels_sug ~ 'Sugars',
                                 TRUE ~ ''),
         new_labels = case_when(new_labels == "Cytidine-3'-monophosphate" ~ "3'-CMP",
                                new_labels == "Cytidine- 2',3'-cyclic monophosphate" ~ "2',3'-cCMP",
                                new_labels == "Uridine-3'-monophosphate" ~ "3'-UMP",
                                new_labels == "Uridine-2',3'-cyclic-monophosphate" ~ "2',3'-cUMP",
                                new_labels == "Cytidine-2'-monophosphate" ~ "2'-CMP",
                                new_labels == "Cytidine-5'-monophosphate" ~ "5'-CMP",
                                new_labels == "Uridine-2'-monophosphate" ~ "2'-UMP",
                                new_labels == "Uridine-5'-monophosphate" ~ "5'-UMP",
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                TRUE ~ new_labels)) %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_hline(aes(yintercept = 1), alpha = 0.4, color = 'grey') +
  geom_point(aes(fill = worm_score_5, alpha = points_alpha), pch=21, size = 5) +
  scale_fill_gradientn(colours = worm_colors,
                       breaks = c(1,2,3,4), limits = c(1,4), 
                       guide = "legend", name = 'Developmental\nscore') +
  scale_color_manual(values = c('black',
                                '#C70B00', # nucleotides 
                                '#310CB3' # sugars
  )) +
  geom_text_repel(aes(x = MetaboliteU, y = score, 
                      label= new_labels),
                  box.padding = 0.2,
                  size = 3.5,
                  segment.alpha = 0.4,
                  max.overlaps = 100) +
  labs(x = 'Metabolite',
       y = 'Normalised \nbacterial growth'
  ) +
  coord_cartesian(xlim = c(-10,415),
                  ylim = c(-0.5,4.5)) +
  scale_size(guide="none") +
  theme_classic() +
  theme(
    # axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_text(size=16, face="bold", color ='black'),
    axis.title.y = element_text(size=16, face="bold", color ='black'),
    plot.title = element_text(size=18),
    legend.text = element_text(size=18)
  ) +
  guides(color = FALSE,
         alpha = FALSE,
         fill = guide_legend(override.aes = list(size=8)))



## pyrE case ----------------------------------------------------------------


pyrE %>% 
  drop_na(worm_score_5) %>% 
  mutate(score = abs(log2FC_5FU - (log2FC_control - adj$slope[2]))) %>% 
  arrange(score) %>% 
  mutate(new_labels = case_when(MetaboliteU %in% labels_4way_nogluc ~ MetaboliteU,
                                TRUE ~ ''),
         points_alpha = case_when(MetaboliteU %in% labels_4way ~ 1,
                                  TRUE ~ 0.15),
         label_size = case_when(worm_score_5 == 4 ~ 3,
                                TRUE ~ 2),
         label_color = case_when(MetaboliteU %in% labels_nuc ~ 'Nucleotide',
                                 MetaboliteU %in%  labels_sug ~ 'Sugars',
                                 TRUE ~ ''),
         new_labels = case_when(new_labels == "Cytidine-3'-monophosphate" ~ "3'-CMP",
                                new_labels == "Cytidine- 2',3'-cyclic monophosphate" ~ "2',3'-cCMP",
                                new_labels == "Uridine-3'-monophosphate" ~ "3'-UMP",
                                new_labels == "Uridine-2',3'-cyclic-monophosphate" ~ "2',3'-cUMP",
                                new_labels == "Cytidine-2'-monophosphate" ~ "2'-CMP",
                                new_labels == "Cytidine-5'-monophosphate" ~ "5'-CMP",
                                new_labels == "Uridine-2'-monophosphate" ~ "2'-UMP",
                                new_labels == "Uridine-5'-monophosphate" ~ "5'-UMP",
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                TRUE ~ new_labels)) %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_hline(aes(yintercept = 1), alpha = 0.4, color = 'grey') +
  geom_point(aes(fill = worm_score_5, alpha = points_alpha), pch=21, size = 5) +
  scale_fill_gradientn(colours = worm_colors,
                       breaks = c(1,2,3,4), limits = c(1,4), 
                       guide = "legend", name = 'Developmental\nscore') +
  scale_color_manual(values = c('black',
                                '#C70B00', # nucleotides 
                                '#310CB3' # sugars
  )) +
  geom_text_repel(aes(x = MetaboliteU, y = score, 
                      label= new_labels),
                  box.padding = 0.2,
                  size = 3.5,
                  segment.alpha = 0.4,
                  max.overlaps = 100) +
  labs(x = 'Metabolite',
       y = 'Normalised \nbacterial growth'
  ) +
  coord_cartesian(xlim = c(-10,415),
                  ylim = c(-0.5,4.5)) +
  scale_size(guide="none") +
  theme_classic() +
  theme(
    # axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_text(size=16, face="bold", color ='black'),
    axis.title.y = element_text(size=16, face="bold", color ='black'),
    plot.title = element_text(size=18),
    legend.text = element_text(size=18)
  ) +
  guides(color = FALSE,
         alpha = FALSE,
         fill = guide_legend(override.aes = list(size=8)))




## pyrE case ----------------------------------------------------------------


TM %>% 
  drop_na(worm_score_250) %>% 
  mutate(score = abs(log2FC_5FU - (log2FC_control - adj$slope[3]))) %>% 
  arrange(score) %>% 
  mutate(new_labels = case_when(MetaboliteU %in% labels_4way_nogluc ~ MetaboliteU,
                                TRUE ~ ''),
         points_alpha = case_when(MetaboliteU %in% labels_4way ~ 1,
                                  TRUE ~ 0.15),
         label_size = case_when(worm_score_250 == 4 ~ 3,
                                TRUE ~ 2),
         label_color = case_when(MetaboliteU %in% labels_nuc ~ 'Nucleotide',
                                 MetaboliteU %in%  labels_sug ~ 'Sugars',
                                 TRUE ~ ''),
         new_labels = case_when(new_labels == "Cytidine-3'-monophosphate" ~ "3'-CMP",
                                new_labels == "Cytidine- 2',3'-cyclic monophosphate" ~ "2',3'-cCMP",
                                new_labels == "Uridine-3'-monophosphate" ~ "3'-UMP",
                                new_labels == "Uridine-2',3'-cyclic-monophosphate" ~ "2',3'-cUMP",
                                new_labels == "Cytidine-2'-monophosphate" ~ "2'-CMP",
                                new_labels == "Cytidine-5'-monophosphate" ~ "5'-CMP",
                                new_labels == "Uridine-2'-monophosphate" ~ "2'-UMP",
                                new_labels == "Uridine-5'-monophosphate" ~ "5'-UMP",
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                TRUE ~ new_labels)) %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_hline(aes(yintercept = 1), alpha = 0.4, color = 'grey') +
  geom_point(aes(fill = worm_score_250, alpha = points_alpha), pch=21, size = 5) +
  scale_fill_gradientn(colours = worm_colors,
                       breaks = c(1,2,3,4), limits = c(1,4), 
                       guide = "legend", name = 'Developmental\nscore') +
  scale_color_manual(values = c('black',
                                '#C70B00', # nucleotides 
                                '#310CB3' # sugars
  )) +
  geom_text_repel(aes(x = MetaboliteU, y = score, 
                      label= new_labels),
                  box.padding = 0.2,
                  size = 3.5,
                  segment.alpha = 0.4,
                  max.overlaps = 100) +
  labs(x = 'Metabolite',
       y = 'Normalised \nbacterial growth'
  ) +
  coord_cartesian(xlim = c(-10,415),
                  ylim = c(-0.5,4.5)) +
  scale_size(guide="none") +
  theme_classic() +
  theme(
    # axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_text(size=16, face="bold", color ='black'),
    axis.title.y = element_text(size=16, face="bold", color ='black'),
    plot.title = element_text(size=18),
    legend.text = element_text(size=18)
  ) +
  guides(color = FALSE,
         alpha = FALSE,
         fill = guide_legend(override.aes = list(size=8)))

