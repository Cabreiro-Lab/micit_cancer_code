# script to generate the plots for the drug-gene-nutrient screen

library(tidyverse)
library(cowplot)
library(readxl)
library(glue)
library(broom)
library(ggrepel)
library(ggtern)

theme_set(theme_cowplot(15))


# read data from the tables

gene_drug = read_excel("data/Table S1.xlsx", sheet = "gene_nutrient_screen",
           skip = 1)


ecocyc = read_csv("ecocyc_paths.csv")

# genes that have an absolute normalised effect of > 0.5
genes2use = c("acnB","gpt","pfkA","pgi","ppk","rpiB","fbaB","guaA","rpe","talA",
              "upp","add","amn","appC","aroB","codA","codB","cpdB","deoB","flgB","flgD",
              "flgF","fliH","fliO","fliQ","fliS","fliT","galE","galU","gcvH",
              "gcvP","gcvT","gpmA","gpp","guaB","holE","ndk","nrdD","nrdE",
              "nrdF","nudE","nudl","nuoB","nuoC","nuoG","nuoH","nuoI","nuoJ",
              "otsA","otsB","pdxA","pdxB","pdxJ","ppx","purA","purC",
              "purF","purL","purM","talB","thrC","tktB","tpiA",
              "udk","udp","ugd","ushA","yqeA","aceF","allC","apaH",
              "aroD","aroK","aspC","atpA","atpH","cyaA","cydB","cyoE",
              "cysC","cysD","cysN","dnaQ","fbp",
              "flhC","fliJ","fliP","folM","fumD","gltA","gnd","gsk","lpd",
              "nuoA","nuoF","pdxH","ppnP","prpC","purD",
              "purE","purH","purK","sdhB","serC","sucB",
              "sucC","sucD","ybhA","allB","aroA","aroC","aroG","carA","carB",
              "flgC","flgG","flgH","flgI","flgN","flhD","fliI",
              "fliR","frdC","holC","mazG","mdh",
              "motA","motB","nudF","nuoK","nuoL","pdxK","pheA",
              "purN","purT","purU","pyrB","pyrC","pyrD","pyrE",
              "pyrF","pyrI","sucA","trpA","trpC","xdhA","yahI",
              "yfbR","yjjG","zwf","nuoN")

paths2represent = c("superpathway of purine nucleotides de novo biosynthesis II",
                    "superpathway of histidine, purine, and pyrimidine biosynthesis",
                    "superpathway of guanosine nucleotides de novo biosynthesis II",
                    "NADH to trimethylamine N-oxide electron transfer",
                    "NADH to cytochrome bd oxidase electron transfer I",
                    "superpathway of adenosine nucleotides de novo biosynthesis II",
                    "flagellar assembly",
                    "salvage pathways of pyrimidine ribonucleotides",
                    "pentose phosphate pathway (non-oxidative branch)",
                    "TCA cycle I (prokaryotic)",
                    "superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
                    "UMP biosynthesis I",
                    "superpathway of pyrimidine ribonucleotides de novo biosynthesis",
                    "superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis",
                    "superpathway of aromatic amino acid biosynthesis",
                    "inosine-5'-phosphate biosynthesis I",
                    "chorismate biosynthesis from 3-dehydroquinate",
                    "salvage pathways of pyrimidine ribonucleotides"
)



drug_resistance_paths = gene_drug %>% 
  filter(Genes %in% genes2use) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>% 
  left_join(ecocyc, relationship = "many-to-many") %>% 
  filter(pathway %in% paths2represent) %>% 
  filter(Supplement_mM == 10) %>% 
  group_by(pathway) %>% 
  summarise(drug_resistance = mean(BW_norm)) %>% 
  distinct(pathway, .keep_all = T)


gene_drug %>% 
  filter(Genes %in% genes2use) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>% 
  left_join(ecocyc ,
            relationship = "many-to-many") %>% 
  filter(pathway %in% paths2represent) %>% 
  filter(Supplement_mM == 0) %>% 
  left_join(drug_resistance_paths) %>% 
  mutate(pathway =  factor(pathway,
                           levels = c(
                             # energy
                             "TCA cycle I (prokaryotic)",
                             "superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
                             "pentose phosphate pathway (non-oxidative branch)",
                             "pentose phosphate pathway",
                             "NADH to trimethylamine N-oxide electron transfer",
                             "NADH to cytochrome bd oxidase electron transfer I",
                             # nucleotides
                             "superpathway of pyrimidine ribonucleotides de novo biosynthesis",
                             "superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis",
                             "UMP biosynthesis I",
                             "superpathway of purine nucleotides de novo biosynthesis II",
                             "inosine-5'-phosphate biosynthesis I",
                             "superpathway of guanosine nucleotides de novo biosynthesis II",
                             "superpathway of adenosine nucleotides de novo biosynthesis II",
                             "salvage pathways of pyrimidine ribonucleotides",
                             "superpathway of histidine, purine, and pyrimidine biosynthesis",
                             "superpathway of aromatic amino acid biosynthesis",
                             "chorismate biosynthesis from 3-dehydroquinate",
                             "flagellar assembly"
                           ))) %>% 
  ggplot(aes(y = pathway, x = BW_norm, fill = drug_resistance)) +
  geom_violin(show.legend = T, 
              scale = "width",
              trim = T,
              alpha = 0.9) +
  geom_jitter(height = 0.2, width = 0.05, show.legend = F,
              shape = 21, size = 4) +
  scale_fill_gradient2(mid = "grey90", high = '#31A95E', low = '#FA2014') +
  labs(
    x = "C. elegans phenotype normalised scores",
    y = NULL,
    fill = "5-FU + Glucose\nresistance"
  ) +
  scale_y_discrete(labels = scales::label_wrap(40), limits=rev) +
  theme_cowplot(17, font_family = "Arial") 

