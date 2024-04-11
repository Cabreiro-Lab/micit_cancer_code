# script to generate the plots for the 4-way screen

library(tidyverse)
library(cowplot)
library(here)
library(rstatix)
library(readxl)
library(glue)
library(broom)

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








