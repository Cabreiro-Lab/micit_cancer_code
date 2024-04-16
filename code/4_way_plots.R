# script to generate the plots for the 4-way screen

library(tidyverse)
library(cowplot)
library(here)
library(rstatix)
library(readxl)
library(glue)
library(broom)
library(ggrepel)
library(ggtern)

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




# ternary plot ------------------------------------------------------------

## calculate scores and join datasets 

BW_scores = BW %>% 
  drop_na(worm_score_5) %>% 
  mutate(score = abs(log2FC_5FU - (log2FC_control - adj$slope[1])) * worm_score_5) %>% 
  select(MetaboliteU, BW = score, wBW = worm_score_5)

pyrE_scores = pyrE %>% 
  drop_na(worm_score_5) %>% 
  mutate(score = abs(log2FC_5FU - (log2FC_control - adj$slope[2])) * worm_score_5) %>% 
  select(MetaboliteU, pyrE = score, wpyrE = worm_score_5)

TM_scores = TM %>% 
  drop_na(worm_score_250) %>% 
  mutate(score = abs(log2FC_5FU - (log2FC_control - adj$slope[3])) * worm_score_250) %>% 
  select(MetaboliteU, TM = score, wTM = worm_score_250)

scores_all = BW_scores %>% 
  left_join(pyrE_scores) %>% 
  left_join(TM_scores)


## find where metabolites rescue worm phenotype and annotate
scores_all = scores_all %>% 
  mutate(Phenotype = case_when(wTM == 4 ~ 'TM',
                               wBW == 4 & wpyrE == 4 ~ 'BW and pyrE',
                               wpyrE == 4 & wBW != 4 ~ 'pyrE',
                               wpyrE != 4 & wBW == 4 ~ 'BW'))

# lists of metabolites belonging to sugars and nucleotides
sugars = c("Dulcitol", "Maltose", "Glycerol",     
           "D-Galactose", "D-Trehalose", "D-Sorbitol")

nucleotides  = c("Uridine-2',3'-cyclic-monophosphate", 
                 "Cytidine-5'-monophosphate", "Uridine-2'-monophosphate", 
                 "Uridine-5'-monophosphate", "Cytidine-2'-monophosphate", 
                 "Uracil_N", "Uridine-3'-monophosphate", 
                 "Cytidine-2',3'-cyclic monophosphate", "Cytidine_N", "Uridine", 
                 "Cytidine-3'-monophosphate", "Thymidine", "Uridine_N")


# plot with ggtern
scores_all %>% 
  mutate(EcoCyc_Classes = case_when(MetaboliteU %in% sugars ~ "Sugars",
                                    MetaboliteU %in% nucleotides ~ "Nucleotides")) %>%
  filter(EcoCyc_Classes %in% c('Sugars','Nucleotides')) %>%
  ggtern(aes(pyrE,BW,TM,
             color = EcoCyc_Classes)) +
  stat_density_tern(geom = 'polygon',
                    size = 0.5,
                    bdl = 0.058,
                    n = 100,
                    aes(alpha=..level..,
                        fill = EcoCyc_Classes),
                    weight = 1,
                    base = "ilr") + 
  geom_point(data = scores_all, 
             pch = 21, 
             color = 'grey70',
             aes(group = Phenotype, 
                 fill = Phenotype)) +
  scale_fill_manual(
    values = c('#DB1D51', # red
               '#FA4F22', # orange
               '#FA754B',  # orange CLASS
               '#0CEB96', # light blue
               '#D6D613', # yellow
               '#3D9EF0', # blue CLASS
               '#5913F2' # dark blue
               
    )
  ) +
  scale_color_manual(
    values = c(
      '#FA754B', # orange
      '#3D9EF0'  # blue
    )
  ) + 
  theme_rgbw() +
  scale_size(range = c(1, 5)) +
  guides(alpha = "none", size = 'none') +
  theme_nomask() +
  labs(color = 'Rescued\nphenotype') + 
  theme(legend.key = element_rect(fill = NA, color = NA))





# Enrichment analysis -----------------------------------------------------


# Nutrient classes --------------------------------------------------------


# prepare info about metabolite classes and pathways
EC_classes = read.csv("KEGGEcoCycClass.csv") %>% 
  separate_rows(EcoCyc_Classes, sep = ';')

EC_path = read.csv("KEGGEcoCycClass.csv") %>%
  separate_rows(Pathways, sep = ';')

## read info about stats from the 4-way screen

screen_stats = read_excel("data/Table S2.xlsx", sheet = "4way_stats") 


# total metab
met = screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "BW_5FU") %>%
  pull(EcoCycID)

# significative metabolites 
sig.met.BW.up =   screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "BW_5FU", 
         FDR < 0.05, logFC > 0) %>% pull(EcoCycID)
sig.met.BW.down = screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "BW_5FU", 
         FDR < 0.05, logFC < 0) %>% pull(EcoCycID)

sig.met.pyrE.up =   screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "pyrE_5FU", 
         FDR < 0.05, logFC > 0) %>% pull(EcoCycID)
sig.met.pyrE.down = screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "pyrE_5FU", 
         FDR < 0.05, logFC < 0) %>% pull(EcoCycID)

sig.met.TM.up =   screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "TM_5FU", 
         FDR < 0.05, logFC > 0) %>% pull(EcoCycID)
sig.met.TM.down = screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "TM_5FU", 
         FDR < 0.05, logFC < 0) %>% pull(EcoCycID)




## Bacterial side ---------------------------------------

N = length(met %in% EC_classes$EcoCycID) # total marked elements

classes = EC_classes %>% 
  filter(Plate != 'extra') %>%
  distinct(EcoCyc_Classes) %>% 
  pull(EcoCyc_Classes)


## helper function 
hyper_classes = function(gene_list){
  enrich = c()
  for (class in classes){
    class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% pull(EcoCycID) 
    m = length(class.met)
    n = N - m
    k = length(gene_list)
    x = length(class.met[class.met %in% gene_list])
    fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
    enrich = c(enrich, fit)
  }
  enrich = p.adjust(enrich, method = 'fdr')
  names(enrich) = classes
  enrich = enrich[enrich < 0.05]
  return(enrich)
}

## perform enrichment ----------
# BW
enrich.BW.up = hyper_classes(sig.met.BW.up)
enrich.BW.down = hyper_classes(sig.met.BW.down)

# pyrE
enrich.pyrE.up = hyper_classes(sig.met.pyrE.up)
enrich.pyrE.down = hyper_classes(sig.met.pyrE.down)

# TM
enrich.TM.up = hyper_classes(sig.met.TM.up)
enrich.TM.down = hyper_classes(sig.met.TM.down)


enrich.BW.up      = data.frame(enrich.BW.up     ) %>% 
  mutate(Class = rownames(.), Direction =   'up', Strain =   'BW');        
enrich.BW.down    = data.frame(enrich.BW.down   ) %>% 
  mutate(Class = rownames(.), Direction = 'down', Strain =   'BW');       
enrich.pyrE.up    = data.frame(enrich.pyrE.up   ) %>% 
  mutate(Class = rownames(.), Direction =   'up', Strain = 'pyrE');       
enrich.pyrE.down  = data.frame(enrich.pyrE.down ) %>% 
  mutate(Class = rownames(.), Direction = 'down', Strain = 'pyrE');       
enrich.TM.up      = data.frame(enrich.TM.up     ) %>% 
  mutate(Class = rownames(.), Direction =   'up', Strain =   'TM');       
enrich.TM.down    = data.frame(enrich.TM.down   ) %>% 
  mutate(Class = rownames(.), Direction = 'down', Strain =   'TM');       


names(enrich.BW.up)[1] = 'p.value'
names(enrich.BW.down)[1] = 'p.value'
names(enrich.pyrE.up)[1] = 'p.value'
names(enrich.pyrE.down)[1] = 'p.value'
names(enrich.TM.up)[1] = 'p.value'
names(enrich.TM.down)[1] = 'p.value'


bac.enrich = rbind(enrich.BW.up, enrich.BW.down, enrich.pyrE.up, enrich.pyrE.down)

# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(bac.enrich$Class))
direct = c('up', 'down')

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes,
                 Direction = c('up', 'down'))

df = left_join(df, bac.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))

# enrichment procedure
enrbrks = c(0, -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.','<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 6)

p.theme = theme(axis.ticks = element_blank(), panel.border = element_blank(), 
                panel.background = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), 
                axis.line = element_line(colour = NA), 
                axis.line.x = element_line(colour = NA), 
                axis.line.y = element_line(colour = NA), 
                strip.text = element_text(colour = "black", 
                                          face = "bold", 
                                          size = 10),
                axis.text.x = element_text(face = "bold", 
                                           colour = "black", size = 10, 
                                           angle = 90, hjust = 1))


# plot enrichment p-values
df %>% 
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, 
                         labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T)),
         Direction = factor(Direction, 
                            levels = sort(direct, decreasing = T))) %>%
  ggplot(aes(x = Direction, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  facet_wrap(~Strain) +
  p.theme


## Worm side ---------------------------------------

N = length(met %in% EC_classes$EcoCycID) # total marked elements

classes = EC_classes %>% 
  filter(Plate != 'extra') %>%
  distinct(EcoCyc_Classes) %>% 
  pull(EcoCyc_Classes)




sig.w.BW = screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "BW_5FU", 
         Median == 4) %>% pull(EcoCycID)
sig.w.pyrE = screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "pyrE_5FU", 
         Median == 4) %>% pull(EcoCycID)
sig.w.TM = screen_stats %>% 
  filter(Metabolite != 'Negative Control', Contrast == "TM_5FU", 
         Median == 4) %>% pull(EcoCycID)


## helper function 
hyper_classes = function(gene_list){
  classes = unique(EC_classes$EcoCyc_Classes)
  enrich = c()
  for (class in classes){
    class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% pull(EcoCycID) 
    m = length(class.met)
    n = N - m
    k = length(gene_list)
    x = length(class.met[class.met %in% gene_list])
    fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
    enrich = c(enrich, fit)
  }
  enrich = p.adjust(enrich, method = 'fdr')
  names(enrich) = classes
  enrich = enrich[enrich < 0.05]
  return(enrich)
}



enrich.BW = hyper_classes(sig.w.BW)
enrich.pyrE = hyper_classes(sig.w.pyrE)
enrich.TM = hyper_classes(sig.w.TM)



enrich.BW      = data.frame(enrich.BW     ) %>% 
  mutate(Class = rownames(.), Strain =   'BW');        
enrich.pyrE    = data.frame(enrich.pyrE   ) %>% 
  mutate(Class = rownames(.), Strain = 'pyrE');       
enrich.TM      = data.frame(enrich.TM     ) %>% 
  mutate(Class = rownames(.), Strain =   'TM');       


names(enrich.BW)[1] = 'p.value'
names(enrich.pyrE)[1] = 'p.value'
names(enrich.TM)[1] = 'p.value'

worm.enrich = rbind(enrich.BW, enrich.pyrE, enrich.TM)

# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(worm.enrich$Class))

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes)

df = left_join(df, worm.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))


# plot enrichment p-values
df %>% 
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, 
                         labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T))) %>%
  ggplot(aes(x = Strain, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  p.theme




# KEGG Pathways --------------------------------------------------------


## Bacterial side ------------------------------------------------------

# helper function
classes = unique(EC_path$Pathways)
classes = classes[-16] # remove an NA
classes = classes[-89]

hyper_path = function(gene_list){
  enrich = c()
  for (class in classes){
    class.met = EC_path %>% filter(Pathways == class) %>% pull(EcoCycID) 
    m = length(class.met)
    n = N - m
    k = length(gene_list)
    x = length(class.met[class.met %in% gene_list])
    fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
    enrich = c(enrich, fit)
  }
  enrich = p.adjust(enrich, method = 'fdr')
  names(enrich) = classes
  enrich = enrich[enrich < 0.05]
  return(enrich)
}


## perform enrichment ----------
# BW
enrich.BW.up = hyper_path(sig.met.BW.up)
enrich.BW.down = hyper_path(sig.met.BW.down)

# pyrE
enrich.pyrE.up = hyper_path(sig.met.pyrE.up)
enrich.pyrE.down = hyper_path(sig.met.pyrE.down)

# TM
enrich.TM.up = hyper_path(sig.met.TM.up)
enrich.TM.down = hyper_path(sig.met.TM.down)

# let's workout the results from each enrichment

enrich.BW.up      = data.frame(enrich.BW.up     ) %>% 
  mutate(Class = rownames(.), Direction =   'up', Strain =   'BW');        
enrich.BW.down    = data.frame(enrich.BW.down   ) %>% 
  mutate(Class = rownames(.), Direction = 'down', Strain =   'BW');       
enrich.pyrE.up    = data.frame(enrich.pyrE.up   ) %>% 
  mutate(Class = rownames(.), Direction =   'up', Strain = 'pyrE');       
enrich.pyrE.down  = data.frame(enrich.pyrE.down ) %>% 
  mutate(Class = rownames(.), Direction = 'down', Strain = 'pyrE');       
enrich.TM.up      = data.frame(enrich.TM.up     ) %>% 
  mutate(Class = rownames(.), Direction =   'up', Strain =   'TM');       
enrich.TM.down    = data.frame(enrich.TM.down   ) %>% 
  mutate(Class = rownames(.), Direction = 'down', Strain =   'TM');       


names(enrich.BW.up)[1] = 'p.value'
names(enrich.BW.down)[1] = 'p.value'
names(enrich.pyrE.up)[1] = 'p.value'
names(enrich.pyrE.down)[1] = 'p.value'
names(enrich.TM.up)[1] = 'p.value'
names(enrich.TM.down)[1] = 'p.value'


bac.enrich = rbind(enrich.BW.up, enrich.BW.down, 
                   enrich.pyrE.up, enrich.pyrE.down)

# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(bac.enrich$Class))
direct = c('up', 'down')

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes,
                 Direction = c('up', 'down'))

df = left_join(df, bac.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))


# plot enrichment p-values
df %>% 
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T)),
         Direction = factor(Direction, levels = sort(direct, decreasing = T))) %>%
  ggplot(aes(x = Direction, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  facet_wrap(~Strain) +
  p.theme






## Worm side ---------------------------------------

# N = length(met %in% EC_classes$EcoCycID) # total marked elements

classes = unique(EC_path$Pathways)

enrich.BW = hyper_path(sig.w.BW)
enrich.pyrE = hyper_path(sig.w.pyrE)
enrich.TM = hyper_path(sig.w.TM)



enrich.BW      = data.frame(enrich.BW     ) %>% 
  mutate(Class = rownames(.), Strain =   'BW');        
enrich.pyrE    = data.frame(enrich.pyrE   ) %>% 
  mutate(Class = rownames(.), Strain = 'pyrE');       
enrich.TM      = data.frame(enrich.TM     ) %>% 
  mutate(Class = rownames(.), Strain =   'TM');       


names(enrich.BW)[1] = 'p.value'
names(enrich.pyrE)[1] = 'p.value'
names(enrich.TM)[1] = 'p.value'

worm.enrich = rbind(enrich.BW, enrich.pyrE, enrich.TM)

# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(worm.enrich$Class))

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes)

df = left_join(df, worm.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))


# plot enrichment p-values
df %>% 
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, 
                         labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T))) %>%
  ggplot(aes(x = Strain, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  p.theme












