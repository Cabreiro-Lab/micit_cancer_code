# libraries -------------------

library(tidyverse)
library(readxl)
library(here)
library(openxlsx)
library(cowplot)
library(broom)
library(plotly)
library(cluster)
library(factoextra)
library(Rtsne)
library(extrafont)

theme_set(theme_cowplot(14))



# load data ---------------------------------------------------------------

rep1 = read_excel("rep1/Biolog 241120.xlsx", sheet = "Data") %>% 
  mutate(Micit = as.factor(Micit),
         Plate = as.factor(Plate))


rep2 = read_excel("rep2/Biolog 170721.xlsx", sheet = "Data") %>% 
  mutate(Micit = as.factor(Micit),
         Plate = as.factor(Plate))

rep3 = read_excel("rep3/Biolog_220721.xlsx", sheet = "Data") %>% 
  mutate(Micit = as.factor(Micit),
         Plate = as.factor(Plate))


drugs = read_excel("biolog_metabolites_cancer.xlsx") %>% 
  mutate(Plate = as.factor(Plate))


# load cancer drugs db from: https://www.anticancerfund.org/en/cancerdrugs-db
cancerdrugsdb = read_delim("cancerdrugsdb.txt", 
           delim = "\t", escape_double = FALSE, 
           trim_ws = TRUE)

# read file with drugbank compounds with FP 
drugbank_fp = read_csv("DrugBank_drugs_FP.csv")




# load metabolites dataset that I created
metabolites = read_excel("biolog_metabolites_cancer.xlsx")

# keep a unique list of metabolites for enrichment analysis
metaU = metabolites %>% 
  distinct(Drug,.keep_all = T) %>% 
  filter(Drug != 'Negative Control')



# metaU %>% 
#   separate_rows(Type, sep = ', ') %>% 
#   distinct(Type) %>% write_csv('mol_types.csv')


rep1 = left_join(rep1, drugs) %>% 
  mutate(Well = as.factor(Well),
         Drug = as.factor(Drug), 
         Replicate = 1, .before = Micit)

rep2 = left_join(rep2, drugs) %>% 
  mutate(Well = as.factor(Well),
         Drug = as.factor(Drug),
         Replicate = 2, .before = Micit)

rep3 = left_join(rep3, drugs) %>% 
  mutate(Well = as.factor(Well),
         Drug = as.factor(Drug),
         Replicate = 3, .before = Micit)


# rep3 fix
# change Micit 1 to 0
# rep3_micit_0 = rep3 %>% filter(Micit == 1) %>% 
#   mutate(Micit = 0, Micit = factor(Micit, levels = c(0,1,5,10)))
# # change Micit 0 to 1
# rep3_micit_1 = rep3 %>% filter(Micit == 0) %>% 
#   mutate(Micit = 1, Micit = factor(Micit, levels = c(0,1,5,10)))
# 
# rep3 = rep3 %>% 
#   filter(Micit %in% c(5,10)) %>% 
#   bind_rows(rep3_micit_0, rep3_micit_1)


# exploration -------------------------------------------------------------

data = bind_rows(rep1,rep2,rep3) %>% 
  mutate(Micit = factor(Micit, 
                        levels = c(0,1,5,10)))

data %>% 
  write_csv('raw_data_biolog.csv')

# how are the Neg Controls behaving?
data %>%
  filter(DrugU %in% c('Negative Control|1','Negative Control|2',
                     'Negative Control|3','Negative Control|4')) %>%
  ggplot(aes(y = Value, x = Micit, fill = Micit)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())

ggsave(here('exploration', 'controls_boxplots.pdf'), height = 10, width = 11)

# summarise controls
controls = data %>% 
  filter(DrugU %in% c('Negative Control|1','Negative Control|2',
                     'Negative Control|3','Negative Control|4')) %>% 
  group_by(Plate, Replicate, Micit) %>% 
  summarise(Mean_control = mean(Value, na.rm = TRUE),
            SD_control = sd(Value, na.rm = TRUE)) %>% 
  ungroup

controls

# control dispersion between neg controls
controls %>% 
  ggplot(aes(x = Micit, y = Mean_control, fill = Micit)) +
  geom_histogram(stat = 'identity') +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_control - SD_control, ymax = Mean_control + SD_control), width = 0.2) +
  facet_grid(vars(Plate), vars(Replicate))

ggsave(here('exploration', 'controls_barplots_replicates.pdf'), height = 10, width = 11)



data %>% 
  mutate(Micit = as.factor(Micit)) %>% 
  filter(DrugU %in% c('Negative Control|1','Negative Control|2',
                      'Negative Control|3','Negative Control|4')) %>% 
  group_by(Plate, Micit) %>% 
  summarise(Mean_control = mean(Value, na.rm = TRUE),
            SD_control = sd(Value, na.rm = TRUE)) %>% 
  ggplot(aes(x = Micit, y = Mean_control, fill = Micit)) +
  geom_histogram(stat = 'identity') +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_control - SD_control, ymax = Mean_control + SD_control), width = 0.2) +
  facet_wrap(vars(Plate))

ggsave(here('exploration', 'controls_barplots.pdf'), height = 10, width = 11)



#### calculate viability ####
### only use as a control the one of micit = 0, 
### as we want a global control for all the treatments. 
data_viability = data %>% 
  left_join(controls %>% filter(Micit == 0) %>% 
              select(Plate, Replicate, Mean_control)) %>% 
  mutate(Viability = (Value / Mean_control) * 100,
         Drug_conc = str_sub(DrugU, -1),
         Drug_conc = as.factor(Drug_conc)) %>% 
  select(Plate, Well, Replicate, Micit, Drug, Drug_conc, 
         Value, Viability, Mean_control) 


drug_list = unique(as.character(data_viability$Drug))


#### data summary ####

data_sum = data_viability %>%  
  group_by(Drug, Drug_conc, Micit) %>% 
  summarise(Mean = mean(Viability),
            SD = sd(Viability)) %>% 
  mutate(Micit = factor(Micit, 
                        levels = c(0,1,5,10)))



# plot mean of drug viability
drug = 'Fluorouracil'
data_sum %>% 
  filter(Drug == drug) %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#F57A20") +
  ggtitle(label = drug) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))





# Synnergy finder ---------------------------------------------------------


sfinder = data_viability %>% 
  rename(Response = Viability,
         Drug2 = Drug,
         Conc1 = Micit,
         Conc2 = Drug_conc) %>% 
  mutate(Drug1 = 'Micit',
         ConcUnit = 'A.U.',
         PairIndex = 0) %>% 
  select(-Well, -Value, -Mean_control, -Plate) %>% 
  arrange(Conc1, Drug2, Conc2) 


viab_ctr = sfinder %>% filter(Drug2 == 'Negative Control') %>%
  group_by(Conc1, Replicate) %>% 
  summarise(Mean = mean(Response))



drug_list = unique(as.character(sfinder$Drug2))

# sfinder_exp = expand_grid(Conc1 = c(0,1,5,10),
#                           Conc2 = c(0,1,2,3,4),
#                           Drug1 = 'Micit',
#                           Drug2 = drug_list,
#                           Replicate = c(1,2,3),
#                           ConcUnit = 'A.U.')

# put the controls in an expand grid with all drugs and concs
sfinder_ctrl_exp = expand_grid(Conc1 = c(0,1,5,10),
                          Conc2 = c(0),
                          Drug1 = 'Micit',
                          Drug2 = drug_list,
                          Replicate = c(1,2,3),
                          ConcUnit = 'A.U.')

sfinder_ctrl_exp = sfinder_ctrl_exp %>% 
  mutate(Conc1 = as.factor(Conc1),
         Conc2 = as.factor(Conc2)) %>% 
  left_join(viab_ctr %>% 
              rename(Response = Mean)) 

sfinder = sfinder %>% 
  bind_rows(sfinder_ctrl_exp) %>%
  filter(Drug2 != 'Negative Control')
  
  
# 
# 
# for (drug in drug_list) {
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 0,]$Response = viab_ctr[1,3]$Mean
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 1,]$Response = viab_ctr[2,3]$Mean
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 5,]$Response = viab_ctr[3,3]$Mean
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 10,]$Response = viab_ctr[4,3]$Mean
# 
# }


for (i in 1:length(drug_list)) {
  sfinder[sfinder$Drug2 == drug_list[i],]$PairIndex = i
}

sfinder = sfinder  %>% arrange(Drug2, Conc1, Conc2) %>% 
  select(-Replicate)

# Save statistical analysis results

list_of_datasets = list('PM-M1'= sfinder)

write.xlsx(list_of_datasets, here('exploration', 'SynergyFinder_grid_3reps_No_NC.xlsx'), 
           colNames = T, rowNames = F, overwrite = TRUE) 







# test how it looks in a heatmap


sfinder_sum = sfinder %>%  
  group_by(Drug2, Drug1, Conc2, Conc1) %>% 
  summarise(Mean = mean(Response),
            SD = sd(Response)) %>% 
  ungroup %>% 
  mutate(Conc2 = factor(Conc2, levels = c(0,1,2,3,4)))


drug_list = unique(as.character(sfinder_sum$Drug2))


# plot mean of drug viability
drug = 'Azaserine'
sfinder_sum %>% 
  filter(Drug2 == drug) %>% 
  mutate(Conc1 = factor(Conc1, levels = c(0,1,5,10)),
         Conc2 = factor(Conc2, levels = c(0,1,2,3,4))) %>% 
  ggplot(aes(x = Conc1, y = Conc2, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#F57A20") +
  ggtitle(label = drug) +
  labs(x = 'Micit (mM)',
       y = 'Query drug (A.U.)') +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))


  
  
# plot a heatmap per drug
# this represents the MEAN Viability
for (drug in drug_list){
    p = sfinder_sum %>% 
      mutate(Conc1 = factor(Conc1, levels = c(0,1,5,10)),
             Conc2 = factor(Conc2, levels = c(0,1,2,3,4))) %>% 
      filter(Drug2 == drug) %>% 
      ggplot(aes(x = Conc1, y = Conc2, fill = Mean)) +
      geom_tile() +
      scale_fill_gradient(name = "Viability",
                          low = "#FFFFFF",
                          high = "#F57A20") +
      ggtitle(label = drug) +
      labs(x = 'Micit (mM)',
           y = 'Query drug (A.U.)') +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
    
    ggsave(plot = p, here('exploration/heatmaps',paste0(drug,'_heatmap.pdf') ),device = 'pdf',
           width = 12, height = 11, units = 'cm')
  }
  
  
  sfinder_sum %>% 
    mutate(Conc1 = factor(Conc1, levels = c(0,1,5,10)),
           Conc2 = factor(Conc2, levels = c(0,1,2,3,4))) %>% 
    ggplot(aes(x = Conc1, y = Conc2, fill = Mean)) +
    geom_tile()  +
    scale_fill_gradient(name = "Viability",
                        low = "#FFFFFF",
                        high = "#F57A20",
                        limits = c(10,150)) +
    labs(x = 'Micit (mM)',
         y = 'Query drug (A.U.)') +
    facet_wrap(~Drug2, ncol = 9) +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
  
  
  ggsave(here('exploration','Complete_heatmap_MeanViability.pdf') ,device = 'pdf',
         width = 40, height = 26, units = 'cm')
  
  
  
  
# zip scores --------------------------------------------------------------


library(readxl)
# 
# result_ZIP = read_excel("exploration/zip_results_rep1/result_ZIP_2020-12-03.xlsx") %>% 
#   rename(Drug = Drug.combination) %>% 
#   mutate(Drug = str_sub(Drug, 1,-9))

result_ZIP = read_excel("exploration/zip_results/result_ZIP_2021-07-30.xlsx") %>% 
  rename(Drug = Drug.combination) %>% 
  mutate(Drug = str_sub(Drug, 1,-9))
  
# THE GOOD ONE
# result_ZIP = read_excel("exploration/zip_results/result_ZIP_2021-07-30.xlsx") %>% 
#   rename(Drug = Drug.combination) %>% 
#   mutate(Drug = str_sub(Drug, 1,-9))


result_ZIP %>% 
  mutate(S_left = Synergy.score - `95% CI`,
         S_right = Synergy.score + `95% CI`) %>% 
  ggplot(aes(y = fct_reorder(Drug, Synergy.score), x = Synergy.score)) + 
  geom_point()

data.zip = data_viability %>% left_join(result_ZIP)  %>% 
  mutate(Syn.direction = case_when(Synergy.score < 0 ~ 'Antagonic',
                                   Synergy.score >= 0 ~ 'Synergic'))

drug = 'Fluorouracil'
data.zip %>% filter(Drug == drug) %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#012345") +
  geom_text(aes(2,2, label=Synergy.score), size = 13) +
  ggtitle(label = drug) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))

for (drug in drug_list){
  p = data.zip %>% filter(Drug == drug) %>% 
    ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
    geom_tile() +
    # scale_fill_gradient(name = "Viability",
    #                     low = "#FFFFFF",
    #                     high = "#012345") +
    scale_fill_gradientn(name = 'Viability',
                         colours = c('#FFFFFF', '#012345'), limits = c(13,130)) +
    geom_text(aes(2,2, label=Synergy.score), size = 13) +
    ggtitle(label = drug) +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
  ggsave(plot = p, here('exploration/heatmaps_zip',filename = paste0(drug,'_heatmap.pdf') ),device = 'pdf',
         width = 12, height = 11, units = 'cm')
}


data.zip %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#012345") +
  facet_wrap(~Drug, ncol = 9) +
  geom_text(aes(2,2, label=Synergy.score, color = Syn.direction), size = 5) +
  scale_color_manual(name = 'Synergy\n direction', values = c('#00A325', '#F0011A')) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))


ggsave(here('exploration',filename = 'Complete_heatmap_ZIP.pdf') ,device = 'pdf',
       width = 40, height = 26, units = 'cm')



# nucleotide stats --------------------------------------------------------



result_ZIP = read_excel("exploration/zip_results/result_ZIP_2021-07-30.xlsx") %>% 
  rename(Drug = Drug.combination) %>% 
  mutate(Drug = str_sub(Drug, 1,-9))



nucl = c('Azathioprine', 'Fluorouracil',
         'Mercaptopurine', 'Thioguanine',
         "5-Fluoro-5'- DeoxyurIdine", 'Zidovudine', 
         'Azacytidine', 'Carmofur',
         'Zidovudine', 'Methotrexate',
         'Floxuridine', "Cytosine-Beta-DArabinofuranoside")

res_ZIP_cats = result_ZIP %>% 
  mutate(Category = case_when(Drug %in% nucl ~ 'Nucleotides',
                              TRUE ~ 'Other'))


res_ZIP_cats %>% filter(Category == 'Nucleotides')


#### Most syn stats ####
library(rstatix)

# homoscedasticity test
CG = subset(res_ZIP_cats, Category == "Nucleotides")$Most.synergistic.area.score
TG = subset(res_ZIP_cats, Category != "Nucleotides")$Most.synergistic.area.score
tidy(var.test(CG, TG))

# as they suffer from heterocedasticity, var.equal to FALSE
nucl_stats = t.test(CG, TG, var.equal = F)
library(broom)
pval = tidy(nucl_stats)$p.value

res_ZIP.sum = res_ZIP_cats %>%
  group_by(Category) %>%
  summarise(Mean = mean(Most.synergistic.area.score),
            SD = sd(Most.synergistic.area.score),
            SEM = SD/sqrt(n()))

## MAIN PLOT
res_ZIP_cats %>%
  ggplot(aes(x = Category, y = Most.synergistic.area.score, 
             fill = Category)) +
  geom_violin() +
  geom_jitter(size = 1.1, alpha = 0.5, height = 0, width = 0.1) +
  geom_point(data = res_ZIP.sum,
             size = 3.5,
             color = 'black',
             aes(x = Category, y = Mean)) +
  geom_errorbar(data = res_ZIP.sum,
                aes(x = Category, y = Mean,
                    ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.03) +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[1], xend = 1, 
                   yend = res_ZIP.sum$Mean[1]),
               linetype="dashed", colour = 'grey50') +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[2], xend = 2, 
                   yend = res_ZIP.sum$Mean[2]),
               linetype="dashed", colour = 'grey50') +
  # pval annotation
  geom_segment(aes(x = 1, xend = 2, y = 16.5, yend = 16.5)) +
  annotate("text", x = 1.5, y = 18,
           label = paste('P-value: ',round(pval, 3))) +
  # annotate("text", x = 1.5, y = 18, label = paste('P-value < ','0.0001')) +
  scale_fill_manual(values = c('#F5EC49',
                                '#3D9CE6'),
                    labels = c('Nucleotide \n antimetabolite',
                               'Other')) + 
  labs(x = 'Drug category',
       y = 'Most synergistic ZIP score') +
  scale_x_discrete(name = '',
                   labels = c("Nucleotides" = "Nucleotide \n antimetabolite",
                              "Other" = "Other")) +
  theme(axis.text.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.text.y = element_text(face = "bold", size = 13, color = 'black'),
        legend.position = 'none')

ggsave(here('exploration', 'violin_nucleotides_biolog_MostSynergyScore.pdf'),
       height = 8, width = 9)

ggsave(here('exploration', 'violin_nucleotides_biolog_MostSynergyScore_poster.pdf'), 
       height = 5, width = 4.5)

res_ZIP.sum %>%
  ggplot(aes(x = Category, y = Mean)) +
  geom_point(size = 5, aes(color = Category)) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1)




#### Syn scores stats ####

# homoscedasticity test
CG = subset(res_ZIP_cats, Category == "Nucleotides")$Synergy.score
TG = subset(res_ZIP_cats, Category != "Nucleotides")$Synergy.score
tidy(var.test(CG, TG))

# as they suffer from heterocedasticity, var.equal to FALSE
nucl_stats = t.test(CG, TG, var.equal = F)
library(broom)
pval = tidy(nucl_stats)$p.value

res_ZIP.sum = res_ZIP_cats %>%
  group_by(Category) %>%
  summarise(Mean = mean(Synergy.score),
            SD = sd(Synergy.score),
            SEM = SD/sqrt(n()))

## MAIN PLOT
res_ZIP_cats %>%
  ggplot(aes(x = Category, y = Synergy.score, fill = Category)) +
  geom_violin() +
  geom_jitter(size = 0.5,height = 0, width = 0.1) +
  geom_point(data = res_ZIP.sum,
             size = 3.5,
             color = 'black',
             aes(x = Category, y = Mean)) +
  geom_errorbar(data = res_ZIP.sum,
                aes(x = Category, y = Mean,
                    ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.03) +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[1], xend = 1, yend = res_ZIP.sum$Mean[1]),
               linetype="dashed", colour = 'grey50') +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[2], xend = 2, yend = res_ZIP.sum$Mean[2]),
               linetype="dashed", colour = 'grey50') +
  # pval annotation
  geom_segment(aes(x = 1, xend = 2, y = 16.5, yend = 16.5)) +
  annotate("text", x = 1.5, y = 18, label = paste('P-value: ',round(pval, 5))) +
  # annotate("text", x = 1.5, y = 18, label = paste('P-value < ','0.0001')) +
  scale_fill_manual(values = c('#F5EC49',
                               '#3D9CE6'),
                    labels = c('Nucleotide \n antimetabolite',
                               'Other')) + 
  labs(x = 'Drug category',
       y = 'ZIP score') +
  scale_x_discrete(name = '',
                   labels = c("Nucleotides" = "Nucleotide \n antimetabolite","Other" = "Other")) +
  theme(axis.text.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.text.y = element_text(face = "bold", size = 13, color = 'black'))

ggsave(here('exploration', 'violin_nucleotides_biolog.pdf'), height = 8, width = 9)




# Heatmap for Tanara ------------------------------------------------------



drug = 'Fluorouracil'
data %>% filter(Drug == drug) %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#FA7235") +
  # ggtitle(label = drug) +
  labs(x = 'Metabolite (mM)',
       y = '5-FU (A.U.)') +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.title.y = element_text(face = "bold", size = 13, color = 'black'))

ggsave(here('exploration', 'micit_heatmap_conference.pdf'), height = 8, width = 9)











# # # # # # # # # #
# # # # # # # # # #
# EXPLORATION ####
# # # # # # # # # #
# # # # # # # # # #



# Fingerprints ---------------------------------------------------

dir.create(here('exploration', 'Fingerprints_results'))

### libraries ####

library(broom)
library(plotly)
library(cluster)
library(factoextra)
library(Rtsne)

## load the Morgan fingerprints produced by the script: smiles2morgan.ipynb 

morgan = read_csv("drug_cancer_biolog_morganFP.csv")


### PCA ####

pca_fit = morgan %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F) # do PCA on scaled data


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5) +
  geom_label(aes(label = Drug), alpha = 0.6) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration/Fingerprints_results', 'PCA_morgan.pdf'), 
       height = 13, width = 14)


pca_fit %>%
  augment(morgan) %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug) %>% 
  add_markers()


pca_fit %>%
  tidy(matrix = "eigenvalues")

# barplot of percent explained by PC
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Fingerprints_results', 'percent_explained.pdf'), 
       height = 9, width = 10)

# cumulative expl by PC 
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = '#56B4E9', alpha = 0.6) +
  geom_text(aes(label = round(PC,2)), nudge_x = 0, nudge_y = 0.05) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Fingerprints_results', 'percent_explained_cumulative.pdf'), 
       height = 9, width = 10)



### t-SNE ####

tsne_fit = morgan %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  Rtsne(check_duplicates = FALSE, pca = FALSE, num_threads = 12,
        normalize = FALSE, max_iter = 2000, 
        perplexity = 20, theta = 0, dims = 3)

# generate data frame from tnse results
tsne.df = data.frame(tsne_fit$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df = tsne.df %>% 
  mutate(Drug = morgan$Drug)

tsne.df %>% 
  ggplot(aes(Dim1, Dim2)) + 
  geom_point(size = 1.5) +
  geom_label(aes(label = Drug), alpha = 0.6) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()




ggsave(here('exploration/Fingerprints_results', 't-SNE_fgroups.pdf'), 
       height = 13, width = 14)




tsne.df %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug) %>% 
  add_markers()



### distances ####

drug_dists = morgan %>% select(where(is.numeric)) %>%  data.frame %>% get_dist
fviz_dist(drug_dists, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

ggsave(here('exploration/Fingerprints_results', 'distances.pdf'), 
       height = 9, width = 10)

### K-means FP ####

morgan_mat = morgan %>% 
  select(where(is.numeric)) 



### silhouette plot ####

kclusts <- 
  tibble(k = 1:80) %>%
  mutate(
    kclust = map(k, ~kmeans(morgan_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, morgan_mat)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 90, nudge_x = 0.3)

ggsave(here('exploration/Fingerprints_results', 'silhouette_plot_kmeans.pdf'), 
       height = 9, width = 10)



### k-means plots ####

# perhaps 17 clusters...

kclust = morgan_mat %>% kmeans(centers = 17)


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster)) %>% 
  select(cluster, Drug) %>% 
  left_join(drugs) %>% 
  distinct(Drug, .keep_all=T) %>% 
  filter(Drug %in% nucl) %>% view




### K-means PCA ####

pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = '#56B4E9', alpha = 0.6) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)


pca_PC = pca_fit %>%
  augment(morgan) %>% 
  # select(Drug, .fittedPC1:.fittedPC89)
  select(Drug, .fittedPC1:.fittedPC30)


pca_pc_mat = pca_PC %>% 
  select(where(is.numeric)) 

kclust = pca_pc_mat %>% 
  kmeans(centers = 10)


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()



pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster)) %>% 
  select(cluster, Drug) %>% 
  left_join(drugs) %>% 
  distinct(Drug, .keep_all=T) %>% 
  filter(Drug %in% nucl) %>% view


kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(pca_pc_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, morgan_mat)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 90, nudge_x = 0.3)

ggsave(here('exploration/Fingerprints_results', 'silhouette_plot_kmeans_PCA.pdf'), 
       height = 9, width = 10)




### K-means t-SNE ####


tsne_mat = tsne.df %>% 
  select(where(is.numeric)) 

kclust = tsne_mat %>% 
  kmeans(centers = 10)


tsne.df %>%
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(Dim1,Dim2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 't-SNE 1',
    y = 't-SNE 2'
  ) +
  theme_half_open(12) + 
  background_grid()



kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(tsne_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 10000, nudge_x = 0.4)

ggsave(here('exploration/Fingerprints_results', 'silhouette_plot_kmeans_tSNE.pdf'), 
       height = 9, width = 10)



# Functional groups -------------------------------------------------------

dir.create(here('exploration', 'Functional_groups_results'))

## load the Morgan fingerprints produced by the script: pyMol2FuncGroup.py 
fgroups = read_csv("func_groups_biolog.csv")

# remove dactinomycin as it's very different from all the others
fgroups = fgroups %>% 
  filter(Drug != 'Dactinomycin')

### PCA ####

pca_fit = fgroups %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F) # do PCA on scaled data


pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5) +
  geom_label(aes(label = Drug)) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration/Functional_groups_results', 'PCA_fgroups.pdf'), 
       height = 13, width = 14)


pca_fit %>%
  augment(fgroups) %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug) %>% 
  add_markers()



pca_fit %>%
  tidy(matrix = "eigenvalues")

# barplot of percent explained by PC
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Functional_groups_results', 'percent_explained.pdf'), 
       height = 9, width = 10)

# cumulative expl by PC 
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = 'black', alpha = 0.7) +
  geom_text(aes(label = round(PC,2)), nudge_x = 0, nudge_y = 0.05) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Functional_groups_results', 'percent_explained_cumulative.pdf'), 
       height = 9, width = 10)




### t-SNE ####

tsne_fit = fgroups %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  Rtsne(check_duplicates = FALSE, pca = FALSE, num_threads = 12,
        normalize = FALSE, max_iter = 2000, 
        perplexity = 30, theta = 0.1, dims = 3)

# generate data frame from tnse results
tsne.df = data.frame(tsne_fit$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df = tsne.df %>% 
  mutate(Drug = fgroups$Drug)

tsne.df %>% 
  ggplot(aes(Dim1, Dim2)) + 
  geom_point(size = 1.5) +
  geom_label(aes(label = Drug), alpha = 0.6) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()


ggsave(here('exploration/Functional_groups_results', 't-SNE_fgroups.pdf'), 
       height = 13, width = 14)




tsne.df %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug) %>% 
  add_markers()






### distances ####

drug_dists = fgroups %>% select(where(is.numeric)) %>%  data.frame %>% get_dist
fviz_dist(drug_dists, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

ggsave(here('exploration/Functional_groups_results', 'distances.pdf'), 
       height = 9, width = 10)






### K-means  FP ####

morgan_mat = fgroups %>% 
  select(where(is.numeric)) 



### silhouette plot ####

kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(morgan_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 90, nudge_x = 0.3)



ggsave(here('exploration/Functional_groups_results', 'silhouette_plot_kmeans.pdf'), 
       height = 9, width = 10)






### k-means plots ####

# around 10 clusters is good!
kclust = morgan_mat %>% kmeans(centers = 9)

pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration/Functional_groups_results', 'PCA_fgroups_10_clusters.pdf'), 
       height = 9, width = 10)

# 3D plot 
pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()


tsne.df %>%
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color = ~cluster) %>% 
  add_markers()




### K-means t-SNE ####


tsne_mat = tsne.df %>% 
  select(where(is.numeric)) 

kclust = tsne_mat %>% 
  kmeans(centers = 8)


tsne.df %>%
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(Dim1,Dim2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 't-SNE 1',
    y = 't-SNE 2'
  ) +
  theme_half_open(12) + 
  background_grid()


tsne.df %>%
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color = ~cluster) %>% 
  add_markers()


kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(tsne_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 50000, nudge_x = 0.4)

ggsave(here('exploration/Functional_groups_results', 'silhouette_plot_kmeans_tSNE.pdf'), 
       height = 9, width = 10)



### SAVE DF CLUSTERS ####

# this object contains the clusters that could be useful to split the data

tsne_clusters = tsne.df %>%
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  as_tibble()




# BIG PCA -----------------------------------------------------------------



dir.create(here('exploration', 'DrugBank_biolog_drugs'))

## load the Morgan fingerprints produced by the script: pyMol2FuncGroup.py 

all_compounds = fgroups %>% 
  mutate(DB = 'Biolog', .before = 'Al_COO') %>% 
  bind_rows(drugbank_fp %>% 
              mutate(DB = 'DrugBank', .before = 'Al_COO'))

# remove dactinomycin as it's very different from all the others
all_compounds = all_compounds %>% 
  filter(Drug != 'Dactinomycin')

### PCA ####

pca_fit = all_compounds %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F) # do PCA on scaled data


pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5, alpha = 0.3) +
  # geom_label(aes(label = Drug)) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration/DrugBank_biolog_drugs', 'PCA_fgroups.pdf'), 
       height = 13, width = 14)


pca_fit %>%
  augment(all_compounds) %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug) %>% 
  add_markers()



pca_fit %>%
  tidy(matrix = "eigenvalues")

# barplot of percent explained by PC
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/DrugBank_biolog_drugs', 'percent_explained.pdf'), 
       height = 9, width = 10)

# cumulative expl by PC 
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = 'black', alpha = 0.7) +
  geom_text(aes(label = round(PC,2)), nudge_x = 0, nudge_y = 0.05) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/DrugBank_biolog_drugs', 'percent_explained_cumulative.pdf'), 
       height = 9, width = 10)




### t-SNE ####

tsne_fit = all_compounds %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  Rtsne(check_duplicates = FALSE, pca = FALSE, num_threads = 12,
        normalize = FALSE, max_iter = 2000, 
        perplexity = 20, theta = 0.1, dims = 3)

# generate data frame from tnse results
tsne.df = data.frame(tsne_fit$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df = tsne.df %>% 
  mutate(Drug = all_compounds$Drug)

tsne.df %>% 
  ggplot(aes(Dim1, Dim2)) + 
  geom_point(size = 1.5) +
  # geom_label(aes(label = Drug), alpha = 0.6) +
  labs(
    x = 't-SNE 1',
    y = 't-SNE 2'
  ) +
  theme_half_open(12) + 
  background_grid()


ggsave(here('exploration/DrugBank_biolog_drugs', 't-SNE_fgroups.pdf'), 
       height = 13, width = 14)



tsne.df %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug) %>% 
  add_markers()




### distances ####

# takes a lot of time

# drug_dists = all_compounds %>% select(where(is.numeric)) %>%  data.frame %>% get_dist
# fviz_dist(drug_dists, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
# 
# ggsave(here('exploration/DrugBank_biolog_drugs', 'distances.pdf'), 
#        height = 19, width = 21)


### K-means with morgan FP ####

morgan_mat = all_compounds %>% 
  select(where(is.numeric)) 



kclust = morgan_mat %>% kmeans(centers = 10)

pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()


# 3D plot 
pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()


tsne.df%>%
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()



# 3D plot with Biolog/DrugBank info
pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~DB) %>% 
  add_markers()

tsne.df%>%
  mutate(DB = all_compounds$DB) %>% 
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~DB) %>% 
  add_markers()





kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(morgan_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, morgan_mat)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 90, nudge_x = 0.3)

ggsave(here('exploration/Functional_groups_results', 'silhouette_plot_kmeans.pdf'), 
       height = 9, width = 10)




### get biolog compounds from general tsne ####

# explore the distribution
tsne.df%>%
  mutate(DB = all_compounds$DB) %>% 
  # mutate(cluster = kclust$cluster, .before = Drug,
  #        cluster = as.factor(cluster)) %>% 
  as_tibble() %>% 
  # filter(DB == 'Biolog') %>% 
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~DB) %>% 
  add_markers()


tsne.df.biolog = tsne.df %>%
  mutate(DB = all_compounds$DB) %>% 
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster)) %>% 
  as_tibble() %>% 
  filter(DB == 'Biolog') 


### k-clust of tSNE_biolog ####


tsne_mat_biolog = tsne.df.biolog %>% select(Dim1:Dim3) %>% data.frame

kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(tsne_mat_biolog, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 10000, nudge_x = 0.4)


n_clusters = 8
kclust = tsne_mat_biolog %>% 
  kmeans(centers = n_clusters)


# tsne.df.biolog %>%
#   mutate(cluster = kclust$cluster, .before = Drug,
#          cluster = as.factor(cluster))  %>% 
#   ggplot(aes(Dim1,Dim2, color = cluster)) + 
#   geom_point(size = 3.5) +
#   labs(
#     x = 't-SNE 1',
#     y = 't-SNE 2'
#   ) +
#   theme_half_open(12) + 
#   background_grid()

tsne.df.biolog = tsne.df.biolog %>%
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))

tsne.df.biolog %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()




# working with clusters ---------------------------------------------------

# a bit of exploration

# n_clusters = max(as.numeric(tsne.df.biolog$cluster))
# n_clusters = 9

tsne.df.biolog %>% 
  group_by(cluster) %>% 
  count() %>% 
  ggplot(aes(x = cluster, y = n)) +
  geom_col(aes(fill = cluster), color = 'black') +
  labs(x = 'Cluster',
       y = 'Elements in cluster',
       caption = 'Number of elements in each cluster as calculated in the t-SNE') +
  guides(fill = 'none')


ggsave(here('exploration/Functional_groups_results', 'elements_in_cluster.pdf'), 
       height = 9, width = 10)


result_ZIP %>% 
  left_join(tsne.df.biolog) %>% 
  drop_na(cluster) %>% 
  # left_join(drugs) %>% 
  ggplot(aes(Synergy.score)) +
  geom_density() +
  # geom_point(aes(y = cluster)) +
  facet_wrap(~cluster)



result_ZIP %>% 
  left_join(tsne.df.biolog) %>% 
  drop_na(cluster) %>% 
  ggplot(aes(y = Synergy.score, x = cluster, fill = cluster)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  guides(fill = 'none') +
  labs(x = 'Cluster',
       y = 'Synergy score', 
       caption = 'Boxplots from the synergy scores per cluster')

ggsave(here('exploration/Functional_groups_results', 'boxplot_cluster.pdf'), 
       height = 9, width = 10)





### stats ####

# let's do some magic 

results_cluster = result_ZIP %>% 
  left_join(tsne.df.biolog) %>% 
  drop_na(cluster)


results_cluster = results_cluster %>% 
  mutate(new_cat = 0,
         new_cat= factor(new_cat)) %>% # create a dummy category column
  bind_rows(results_cluster %>% 
              mutate(new_cat = cluster)) %>% 
  mutate(new_cat = as.factor(new_cat)) 



results_cluster %>% 
  mutate( new_cat = recode(new_cat, `0` = 'ALL')) %>%
  ggplot(aes(y = Synergy.score, x = new_cat, fill = new_cat)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  guides(fill = 'none') +
  labs(x = 'Cluster',
       y = 'Synergy score', 
       caption = 'Boxplots from the synergy scores per cluster')

ggsave(here('exploration/DrugBank_biolog_drugs', 'boxplot_cluster_ALL.pdf'), 
       height = 9, width = 10)



# helper function
adhoc_stats = function(cluster = 1, mode = "Synergy.score"){
  
  ### Takes the dataset I created, filters it, and calculates the stats
  
  form = paste(mode, '~ new_cat')
  form = formula(form)
  
  # print(glue::glue('Mode {mode} activated'))
  
  pair = c(0,cluster)
  
  res = results_cluster %>% 
    filter(new_cat %in% pair) %>% 
    # group_by(new_cat) %>% 
    nest(data = everything()) %>% 
    mutate(model = map(data, lm, formula = form),
           model_tidy = map(model, tidy)) %>% 
    select(model_tidy) %>% 
    unnest(cols = c(model_tidy)) %>% 
    filter(term != '(Intercept)') %>% 
    mutate(term = cluster)
  
  return(res)
    
}



# helper function for Welch test
adhoc_welch_stats = function(cluster = 1){
  
  ### Takes the dataset I created, filters it, and calculates the stats
  
  query = cluster
  basal = c(0)
  
  CG = results_cluster %>% 
    # mutate(new_cat = as.numeric(new_cat)) %>% 
    filter(new_cat == query) %>% pull(Synergy.score)
  TG = results_cluster %>% 
    # mutate(new_cat = as.numeric(new_cat)) %>% 
    filter(new_cat == basal) %>% pull(Synergy.score)
  
  # as they suffer from heterocedasticity, var.equal to FALSE
  welch_stats = tidy(t.test(CG, TG, var.equal = F))
  
  return(welch_stats)
  
}



adhoc_stats(cluster = n_clusters)


df_lm = tibble()
for (k in 1:n_clusters){
  
  res = adhoc_stats(cluster = k)
  
  df_lm = df_lm %>% bind_rows(res)

  }

df_lm

df_lm %>% 
  write_csv(here('exploration/DrugBank_biolog_drugs', 'clusters_stats.csv'))

# save the original file to see the clusters
tsne.df.biolog %>% 
  write_csv(here('exploration/DrugBank_biolog_drugs', 'biolog_with_clusters.csv'))


adhoc_welch_stats(cluster = 1)


df_welch = tibble()
for (k in 1:n_clusters){
  
  res = adhoc_welch_stats(cluster = k)
  
  df_welch = df_welch %>% bind_rows(res)
  
}

df_welch






# with metab classification -----------------------------------------------

library(Rtsne)
library(plotly)


# generate one-hot encoding of molecule type
onehot_mol = metaU %>% 
  separate_rows(Molecule, sep = ', ') %>% 
  mutate(val=1) %>% 
  select(Drug, Molecule, val) %>% 
  pivot_wider(names_from = Molecule, values_from = val, values_fill = 0)


# generate one-hot encoding of molecule type
onehot_type = metaU %>% 
  separate_rows(Type, sep = ', ') %>% 
  mutate(val=1) %>% 
  select(Drug, Type, val) %>% 
  pivot_wider(names_from = Type, values_from = val, values_fill = 0)


# generate one-hot encoding of process
onehot_proc = metaU %>% 
  separate_rows(Process, sep = ', ') %>% 
  mutate(val=1) %>% 
  select(Drug, Process, val) %>% 
  pivot_wider(names_from = Process, values_from = val, values_fill = 0)


# generate one-hot encoding of process
onehot_target = metaU %>% 
  separate_rows(Target, sep = ', ') %>% 
  mutate(val=1) %>% 
  select(Drug, Target, val) %>% 
  pivot_wider(names_from = Target, values_from = val, values_fill = 0)



### t-SNE ####

tsne_fit = onehot_type %>% 
  bind_cols(onehot_proc) %>%
  # bind_cols(onehot_mol) %>%
  bind_cols(onehot_target) %>%
  rename(Drug = `Drug...1`) %>%
  # filter(Drug != 'Dactinomycin') %>%
  # bind_cols(fgroups) %>%
  select(where(is.numeric)) %>% # retain only numeric columns
  Rtsne(check_duplicates = FALSE, pca = FALSE, num_threads = 8,
        normalize = FALSE, max_iter = 2000, 
        perplexity = 15, theta = 0.001, dims = 3)

# generate data frame from tnse results
tsne.df = data.frame(tsne_fit$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df = tsne.df %>% 
  mutate(Drug = onehot_type$Drug) %>% 
  mutate(Category = case_when(Drug %in% nucl ~ 'Nucleotide',
                              TRUE ~ 'Other'))


tsne.df %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color = ~Category) %>% 
  add_markers()



### k-means #####

# tsne_mat_biolog = scale(tsne_mat_biolog)


tsne_mat_biolog = tsne.df[,1:3] %>% as.data.frame

rownames(tsne_mat_biolog) = onehot_proc$Drug

res.km = eclust(tsne_mat_biolog, 
                "kmeans", 
                stand = F,
                nstart = 5,
                k.max = 12, 
                seed = 123)

fviz_gap_stat(res.km$gap_stat)

fviz_silhouette(res.km)
# 
# ggsave(here('exploration/Multivariate_final', 'silhouette_analysis.pdf'), 
#        height = 9, width = 10)

### exploration clusters ####

tsne.df.biolog  = tsne_mat_biolog %>% as_tibble(rownames = 'Drug')


tsne.df.biolog = tsne.df.biolog %>%
  mutate(cluster = res.km$cluster, .before = Drug,
         cluster = as.factor(cluster))

tsne.df.biolog %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()




### boxplots ####

# a bit of exploration

# n_clusters = max(as.numeric(tsne.df.biolog$cluster))
# n_clusters = res.km$nbclust

tsne.df.biolog %>% 
  group_by(cluster) %>% 
  count() %>% 
  ggplot(aes(x = cluster, y = n)) +
  geom_col(aes(fill = cluster), color = 'black') +
  labs(x = 'Cluster',
       y = 'Elements in cluster',
       caption = 'Number of elements in each cluster as calculated in the t-SNE') +
  guides(fill = 'none')


# ggsave(here('exploration/Multivariate_final', 'elements_in_cluster.pdf'), 
       # height = 9, width = 10)



### boxplots vs all measures together ####

results_cluster = result_ZIP %>% 
  left_join(tsne.df.biolog) %>% 
  drop_na(cluster)


results_cluster = results_cluster %>% 
  mutate(new_cat = 0,
         new_cat= factor(new_cat)) %>% # create a dummy category column
  bind_rows(results_cluster %>% 
              mutate(new_cat = cluster)) %>% 
  mutate(new_cat = as.factor(new_cat)) 



results_cluster %>% 
  mutate( new_cat = recode(new_cat, `0` = 'ALL')) %>%
  ggplot(aes(y = Most.synergistic.area.score, x = new_cat, fill = new_cat)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  guides(fill = 'none') +
  labs(x = 'Cluster',
       y = 'Most Synergistic area score', 
       caption = 'Boxplots from the synergy scores per cluster')

ggsave(here('exploration/Multivariate_final', 'boxplots_most_synergistic.pdf'), 
       height = 9, width = 10)




### STATS ####

n_clusters = res.km$nbclust


# stats for synergy score
df_lm = tibble()
for (k in 1:n_clusters){
  
  res = adhoc_stats(cluster = k)
  
  df_lm = df_lm %>% bind_rows(res)
  
}


df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value))

df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value)) %>% 
  write_csv(here('exploration/Multivariate_final', 'clusters_synergy_stats.csv'))


####

# stats for Most synergy region score
df_lm = tibble()
for (k in 1:n_clusters){
  
  res = adhoc_stats(cluster = k, mode = 'Most.synergistic.area.score')
  
  df_lm = df_lm %>% bind_rows(res)
  
}


df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value))

df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value)) %>% 
  write_csv(here('exploration/Multivariate_final', 'clusters_MostSynergy_stats.csv'))





# save the original file to see the clusters
# tsne.df.biolog %>% 
#   write_csv(here('exploration/Multivariate_final', 'biolog_with_clusters.csv'))



### mol enrichment cluster ####


## helper function (originally from )

meta_type = metaU %>% 
  separate_rows(Type, sep = ',') %>% 
  select(Drug, Type) %>% 
  drop_na(Drug)


cl1 = results_cluster %>% filter(cluster == 1) %>% pull(Drug)
cl1.enrich = enrich(cl1, db = meta_type) %>% 
  mutate(direction = 'cluster_1')

cl2 = results_cluster %>% filter(cluster == 2) %>% pull(Drug)
cl2.enrich = enrich(cl2, db = meta_type) %>% 
  mutate(direction = 'cluster_2')

cl3 = results_cluster %>% filter(cluster == 3) %>% pull(Drug)
cl3.enrich = enrich(cl3, db = meta_type) %>% 
  mutate(direction = 'cluster_3')

cl4 = results_cluster %>% filter(cluster == 4) %>% pull(Drug)
cl4.enrich = enrich(cl4, db = meta_type) %>% 
  mutate(direction = 'cluster_4')

cl5 = results_cluster %>% filter(cluster == 5) %>% pull(Drug)
cl5.enrich = enrich(cl5, db = meta_type) %>% 
  mutate(direction = 'cluster_5')

cl6 = results_cluster %>% filter(cluster == 6) %>% pull(Drug)
cl6.enrich = enrich(cl6, db = meta_type) %>% 
  mutate(direction = 'cluster_6')

cl7 = results_cluster %>% filter(cluster == 7) %>% pull(Drug)
cl7.enrich = enrich(cl7, db = meta_type) %>% 
  mutate(direction = 'cluster_7')

cl8 = results_cluster %>% filter(cluster == 8) %>% pull(Drug)
cl8.enrich = enrich(cl8, db = meta_type) %>% 
  mutate(direction = 'cluster_8')


cl1.enrich %>% 
  bind_rows(cl2.enrich, cl3.enrich, cl4.enrich, cl5.enrich,
            cl6.enrich, cl7.enrich, cl8.enrich) %>% 
  # filter(pval <= 0.05) %>% 
  view


# # # # # # # # # # # # # #
# # # # # # # # # # # # # #
# # # # # # # # # # # # # # 
# Final version ###########
# # # # # # # # # # # # # #
# # # # # # # # # # # # # #
# # # # # # # # # # # # # #

dir.create(here('exploration', 'Multivariate_final'))

library(broom)
library(plotly)
library(cluster)
library(factoextra)
library(Rtsne)



# the selected strategy is to use the functional groups from drugbank
# and do a t-SNE, select the optimal number of clusters and then get
# the stats

# this is a clean version of the exploration part


### get data ####

## load the Morgan fingerprints produced by the script: pyMol2FuncGroup.py 

all_compounds = fgroups %>% 
  mutate(DB = 'Biolog', .before = 'Al_COO') %>% 
  bind_rows(drugbank_fp %>% 
              mutate(DB = 'DrugBank', .before = 'Al_COO'))

# remove dactinomycin as it's very different from all the others
all_compounds = all_compounds %>% 
  filter(Drug != 'Dactinomycin')



### t-SNE ####

tsne_fit = all_compounds %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  Rtsne(check_duplicates = FALSE, pca = FALSE, num_threads = 8,
        normalize = FALSE, max_iter = 2000, 
        perplexity = 20, theta = 0.5, dims = 3)

# generate data frame from tnse results
tsne.df = data.frame(tsne_fit$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df = tsne.df %>% 
  mutate(Drug = all_compounds$Drug)


tsne.df %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug) %>% 
  add_markers()




### k-clust of tSNE_biolog ####

### get biolog compounds from the previous tSNE ####

# explore the distribution



t1 = list(family = 'Arial', size = 22)

tsne.df %>%
  mutate(DB = all_compounds$DB) %>% 
  as_tibble() %>% 
  mutate(opacity = case_when(DB == 'Biolog' ~ 1,
                             TRUE ~ 0.2)) %>% 
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~DB, 
          opacity=~opacity,
          colors = 'Set1') %>% 
  add_markers() %>% 
  layout(
    font = t1,
    scene = list(
      yaxis = list(gridwidth = 1.5,
                   title = ""),
      xaxis = list(gridwidth = 1.5,
                   title = ""),
      zaxis = list(gridwidth = 1.5,
                   title = "")
    )
  ) 










# keep only the compounds from biolog
tsne.df.biolog = tsne.df %>%
  mutate(DB = all_compounds$DB) %>% 
  as_tibble() %>% 
  filter(DB == 'Biolog') 


### silhouette plot ####
tsne_mat_biolog = tsne.df.biolog %>% select(Dim1:Dim3) %>% data.frame


# kclusts <- 
#   tibble(k = 1:40) %>%
#   mutate(
#     kclust = map(k, ~kmeans(tsne_mat_biolog, .x)),
#     tidied = map(kclust, tidy),
#     glanced = map(kclust, glance)
#   )
# 
# clusterings <- 
#   kclusts %>%
#   unnest(cols = c(glanced))

# average of several runs of silhouette plots
total_clusterings = tibble()
for (i in 1:30){
  kclusts <- 
    tibble(k = 1:40) %>%
    mutate(
      kclust = map(k, ~kmeans(tsne_mat_biolog, .x)),
      tidied = map(kclust, tidy),
      glanced = map(kclust, glance)
    )
  
  clusterings <- 
    kclusts %>%
    unnest(cols = c(glanced)) %>% 
    mutate(replicate = i)
  
  total_clusterings = total_clusterings %>%
    bind_rows(clusterings)

}

total_clusterings %>% 
  group_by(k) %>% 
  summarise(mean.withinss = mean(tot.withinss, na.rm = T),
            sd.withinss = sd(tot.withinss, na.rm = T)) %>% 
  ggplot(aes(k, mean.withinss)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean.withinss - sd.withinss, 
                    ymax = mean.withinss + sd.withinss,
                    width = 0.2)) +
  geom_text(aes(label = k), nudge_y = 10000, nudge_x = 0.4)

ggsave(here('exploration/Multivariate_final', 'silhouette_plot_30iters.pdf'), 
       height = 9, width = 10)


#### silhouette analysis ####

library(factoextra)

# keep only the compounds from biolog
tsne_mat_biolog = tsne.df.biolog %>% select(Dim1:Dim3) %>% data.frame

tsne_mat_biolog
rownames(tsne_mat_biolog) = tsne.df.biolog$Drug

# tsne_mat_biolog = scale(tsne_mat_biolog)

res.km = eclust(tsne_mat_biolog, 
                "kmeans", 
                stand = F,
                nstart = 5,
                k.max = 8, 
                seed = 1234567)

fviz_gap_stat(res.km$gap_stat)

fviz_silhouette(res.km)

ggsave(here('exploration/Multivariate_final', 'silhouette_analysis.pdf'), 
       height = 9, width = 10)

# reset the variable
tsne_mat_biolog = tsne.df.biolog %>% select(Dim1:Dim3) %>% data.frame


### k-means ####
# select the number of clusters

# TODO: apply the analysis from fviz instead of kmeans function

# n_clusters = 8
# kclust = tsne_mat_biolog %>% 
#   kmeans(centers = n_clusters)



tsne.df.biolog = tsne.df.biolog %>%
  mutate(cluster = res.km$cluster, .before = Drug,
         cluster = as.factor(cluster))

tsne.df.biolog %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()




### boxplots ####

# a bit of exploration

# n_clusters = max(as.numeric(tsne.df.biolog$cluster))
# n_clusters = res.km$nbclust

tsne.df.biolog %>% 
  group_by(cluster) %>% 
  count() %>% 
  ggplot(aes(x = cluster, y = n)) +
  geom_col(aes(fill = cluster), color = 'black') +
  labs(x = 'Cluster',
       y = 'Elements in cluster',
       caption = 'Number of elements in each cluster as calculated in the t-SNE') +
  guides(fill = 'none')


ggsave(here('exploration/Multivariate_final', 'elements_in_cluster.pdf'), 
       height = 9, width = 10)



### boxplots vs all measures together ####

results_cluster = result_ZIP %>% 
  left_join(tsne.df.biolog) %>% 
  drop_na(cluster)


results_cluster = results_cluster %>% 
  mutate(new_cat = 0,
         new_cat= factor(new_cat)) %>% # create a dummy category column
  bind_rows(results_cluster %>% 
              mutate(new_cat = cluster)) %>% 
  mutate(new_cat = as.factor(new_cat)) 



results_cluster %>% 
  mutate( new_cat = recode(new_cat, `0` = 'ALL')) %>%
  ggplot(aes(y = Most.synergistic.area.score, x = new_cat, fill = new_cat)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  guides(fill = 'none') +
  labs(x = 'Cluster',
       y = 'Most Synergistic area score', 
       caption = 'Boxplots from the synergy scores per cluster')

ggsave(here('exploration/Multivariate_final', 'boxplots_most_synergistic.pdf'), 
       height = 9, width = 10)



# ### plot with p vals inside
# 
# rstatix::t_test(results_cluster %>% 
#                   filter(new_cat != 0),
#                 Most.synergistic.area.score ~ cluster, 
#                 ref.group = "all",
#                 p.adjust.method = 'none')
# 
# 


### STATS ####

n_clusters = res.km$nbclust


# stats for synergy score
df_lm = tibble()
for (k in 1:n_clusters){
  
  res = adhoc_stats(cluster = k)
  df_lm = df_lm %>% bind_rows(res)
  
}


df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value))

df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value)) %>% 
  write_csv(here('exploration/Multivariate_final', 'clusters_synergy_stats.csv'))


####

# stats for Most synergy region score
df_lm = tibble()
for (k in 1:n_clusters){
  
  res = adhoc_stats(cluster = k, mode = 'Most.synergistic.area.score')
  
  df_lm = df_lm %>% bind_rows(res)
  
}


df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value))

df_lm %>% 
  mutate(p.stars = gtools::stars.pval(p.value)) %>% 
  write_csv(here('exploration/Multivariate_final', 'clusters_MostSynergy_stats.csv'))





# save the original file to see the clusters
tsne.df.biolog %>% 
  write_csv(here('exploration/Multivariate_final', 'biolog_with_clusters.csv'))



### mol enrichment cluster ####

# TODO: review the enrichment calculations!

## helper function (originally from )

enrich = function(gene, db){
  # initiate variables
  pval = c()
  m_total = c()
  x_total = c()
  k_total = c()
  gene_in_cat = c()
  db = as.data.frame(db)
  cats = unique(db[,2])
  
  for (cat in cats){
    subcat = db[db[,2] == cat,]
    N = (db %>% distinct(.[,1]) %>% count())$n
    m = dim(subcat)[1]
    n = N - m
    x = sum(gene %in% subcat[,1])
    k = sum(gene %in% db[,1]) # genes with at least 1 annotation!
    p = phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
    
    # save variables
    m_total = c(m_total, m)
    x_total = c(x_total, x)
    k_total = c(k_total, k)
    gene_in_cat = c(gene_in_cat)
    pval = c(pval, p)
  }
  
  # build the table
  table = tibble(categories = cats, N = N, elm_in_cat = m_total, gene_in_cat = x_total, k_tot = k, pval = pval) %>% 
    mutate(
      p.stars = gtools::stars.pval(pval),
      fdr = p.adjust(pval, method = 'fdr'),
      fdr.stars = gtools::stars.pval(fdr)
    ) %>% 
    arrange(pval)
  
  return(table)
  
}


enrich_metabolites = function(category = 'Target') {
      
      categories = c('Target', 'Process', 'Type', 'Molecule')  
  
      if (category %in% categories){
        
        cat(glue::glue("Calculating the enrichment for {category} category \n\n"))
        
        meta_type = metaU %>% 
          separate_rows(category, sep = ', ') %>% 
          select(Drug, category) %>% 
          drop_na(Drug) %>% 
          drop_na(category) 
      } else {
        print(glue::glue("You need to specify one of the following categories: {categories}"))
        break
      }
      
      # iterate over the clusters created
      enrichment_table = tibble()
      for (cls in 1:max(as.integer(results_cluster$cluster))) {
        print(glue::glue("Cluster {cls}"))
        # filter cluster
        cl = results_cluster %>% filter(cluster == cls) %>% pull(Drug)
        # print(cl)
        # calculate enrich
        cl.enrich = enrich(cl, db = meta_type) %>% 
          mutate(direction = glue::glue('cluster_{cls}'))
        # merge datasets
        enrichment_table = bind_rows(enrichment_table, cl.enrich)
      }
      
      return(enrichment_table)
}



# results_cluster_copy = results_cluster

results_cluster = read_csv('exploration/Multivariate_final/biolog_with_clusters.csv')

enrich_molecule = enrich_metabolites("Molecule")
enrich_target = enrich_metabolites("Target")
enrich_type = enrich_metabolites("Type")
enrich_process = enrich_metabolites("Process")



list_of_datasets = list(
  "molecule enrichment" = enrich_molecule,
  "target enrichment" = enrich_target,
  "type enrichment" = enrich_type,
  "process enrichment" = enrich_process
)


write.xlsx(list_of_datasets, "exploration/Multivariate_final/tSNE_embeddings_enrichment.xlsx")



# stats of each category vs others ----------------------------------------

clusters = read_csv('exploration/Multivariate_final/biolog_with_clusters.csv')

## Type ####

list_of_types = metaU %>% 
    separate_rows(Type, sep = ', ') %>% 
    distinct(Type) %>% 
    drop_na() %>% 
    pull(Type)

types_stats = tibble()
for (type in list_of_types) {
  
  temp_df = result_ZIP %>% 
    left_join(metaU, by = "Drug") %>% 
    mutate(group = case_when(str_detect(Type, type) ~ type,
                             TRUE ~ 'null'), 
           .before=Synergy.score) 
  
  if (length(temp_df$group[temp_df$group == type]) > 2){
    # print(glue::glue("Type {type} is valid"))
    temp_stats = temp_df %>%
      rstatix::t_test(Most.synergistic.area.score ~ group)
    
    types_stats = bind_rows(types_stats, temp_stats)
    
    temp_stats = temp_df %>%
      rstatix::t_test(Synergy.score ~ group)
    
    types_stats = bind_rows(types_stats, temp_stats)
    
    
  } 
  
}


## Molecule ####


list_of_types = metaU %>% 
  separate_rows(Molecule, sep = ', ') %>% 
  distinct(Molecule) %>% 
  drop_na() %>% 
  pull(Molecule)

mol_stats = tibble()
for (type in list_of_types) {
  
  temp_df = result_ZIP %>% 
    left_join(metaU, by = "Drug") %>% 
    mutate(group = case_when(str_detect(Molecule, type) ~ type,
                             TRUE ~ 'null'), 
           .before=Synergy.score) 
  
  if (length(temp_df$group[temp_df$group == type]) > 2){
    # print(glue::glue("Type {type} is valid"))
    temp_stats = temp_df %>%
      rstatix::t_test(Most.synergistic.area.score ~ group)
    
    mol_stats = bind_rows(mol_stats, temp_stats)
    
    temp_stats = temp_df %>%
      rstatix::t_test(Synergy.score ~ group)
    
    mol_stats = bind_rows(mol_stats, temp_stats)
    
  } 
  
}


## Process ####

list_of_types = metaU %>% 
  separate_rows(Process, sep = ', ') %>% 
  distinct(Process) %>% 
  drop_na() %>% 
  pull(Process)

proc_stats = tibble()
for (type in list_of_types) {
  
  temp_df = result_ZIP %>% 
    left_join(metaU, by = "Drug") %>% 
    mutate(group = case_when(str_detect(Process, type) ~ type,
                             TRUE ~ 'null'), 
           .before=Synergy.score) 
  print(length(temp_df$group[temp_df$group == type]))
        
  if (length(temp_df$group[temp_df$group == type]) > 2){
    # print(glue::glue("Type {type} is valid"))
    temp_stats = temp_df %>%
      rstatix::t_test(Most.synergistic.area.score ~ group)
    
    proc_stats = bind_rows(proc_stats, temp_stats)
    
    temp_stats = temp_df %>%
      rstatix::t_test(Synergy.score ~ group)
    
    proc_stats = bind_rows(proc_stats, temp_stats)
    
    
  } 
  
}

## Target ####

list_of_types = metaU %>% 
  separate_rows(Target, sep = ', ') %>% 
  distinct(Target) %>% 
  drop_na() %>% 
  pull(Target)

target_stats = tibble()
for (type in list_of_types) {
  
  temp_df = result_ZIP %>% 
    left_join(metaU, by = "Drug") %>% 
    mutate(group = case_when(str_detect(Target, type) ~ type,
                             TRUE ~ 'null'), 
           .before=Synergy.score) 
  
  if (length(temp_df$group[temp_df$group == type]) > 2){
    # print(glue::glue("Type {type} is valid"))
    temp_stats = temp_df %>%
      rstatix::t_test(Most.synergistic.area.score ~ group)
    
    target_stats = bind_rows(target_stats, temp_stats)
    
    temp_stats = temp_df %>%
      rstatix::t_test(Synergy.score ~ group)
    
    target_stats = bind_rows(target_stats, temp_stats)
    
    
  } 
  
}


result_ZIP %>% 
  left_join(metaU) %>% 
  mutate(group = case_when(str_detect(Process, 'mitosis') ~ 'mitosis',
                           TRUE ~ 'null'), 
         .before=Synergy.score) %>% 
  ggplot(aes(y= Most.synergistic.area.score, x = group)) +
  geom_boxplot() +
  geom_jitter(position = position_dodge(width = 1))
  

## save the datasets

types_stats = types_stats %>% 
  rstatix::add_significance("p") %>% 
  arrange(`.y.`)

mol_stats = mol_stats %>% 
  rstatix::add_significance("p") %>% 
  arrange(`.y.`)

types_stats = proc_stats %>% 
  rstatix::add_significance("p") %>% 
  arrange(`.y.`)

target_stats = target_stats %>% 
  rstatix::add_significance("p") %>% 
  arrange(`.y.`)

list_of_datasets = list(
  types_stats = types_stats,
  mol_stats = mol_stats,
  proc_stats = proc_stats,
  target_stats = target_stats
)

write.xlsx(list_of_datasets, 
           here('exploration','ZIP_categories_stats.xlsx'))




# description categories --------------------------------------------------

metaU %>% 
  separate_rows(Type, sep =  ', ') %>% 
  count(Type) %>%
  drop_na() %>% 
  ggplot(aes(y = fct_reorder(Type, n), x = n)) +
  geom_col(fill = 'grey50', color = 'white') +
  labs(
    x = 'Ocurrence',
    y = 'Molecule type'
  )

ggsave('exploration/molecule_type_metadata.pdf',
       height = 8, width = 5)


metaU %>% 
  separate_rows(Target, sep =  ', ') %>% 
  count(Target) %>%
  drop_na() %>% 
  ggplot(aes(y = fct_reorder(Target, n), x = n)) +
  geom_col(fill = 'grey50', color = 'white') +
  labs(
    x = 'Ocurrence',
    y = 'Target type'
  )

ggsave('exploration/molecule_target_metadata.pdf',
       height = 9, width = 5)

metaU %>% 
  separate_rows(Process, sep =  ', ') %>% 
  count(Process) %>%
  drop_na() %>% 
  ggplot(aes(y = fct_reorder(Process, n), x = n)) +
  geom_col(fill = 'grey50', color = 'white') +
  labs(
    x = 'Ocurrence',
    y = 'Process'
  )

ggsave('exploration/molecule_process_metadata.pdf',
       height = 9, width = 5)







# Synergy plot PAPER ------------------------------------------------------


#### Syn scores stats ####

# homoscedasticity test
CG = subset(res_ZIP_cats, Category == "Nucleotides")$Synergy.score
TG = subset(res_ZIP_cats, Category != "Nucleotides")$Synergy.score
tidy(var.test(CG, TG))

# as they suffer from heterocedasticity, var.equal to FALSE

nucl_stats = t.test(CG, TG, var.equal = F)

nucl_stats %>% tidy %>% write_csv("exploration/biolog_ZIP_stats_PAPER.csv")

library(broom)
pval = tidy(nucl_stats)$p.value

res_ZIP.sum = res_ZIP_cats %>%
  group_by(Category) %>%
  summarise(Mean = mean(Synergy.score),
            SD = sd(Synergy.score),
            SEM = SD/sqrt(n()))



clusters = read_csv('exploration/Multivariate_final/biolog_with_clusters.csv')

# clusters = results_cluster %>% distinct(Drug, .keep_all = T)


library(randomcoloR)

colors = randomColor(7)

colors = c("#ed2863", "#87ffe3", "#cd4cce", "#a57e11" ,"#fc33b2", "#eac885", "#1b6bea")

## MAIN PLOT
res_ZIP_cats %>%
  left_join(clusters) %>% 
  drop_na(cluster) %>% 
  mutate(cluster = factor(cluster),
         plot_labs = case_when(Category == "Nucleotides" ~ Drug)) %>% 
  ggplot(aes(x = Category, y = Synergy.score, fill = Category)) +
  geom_violin() +
  geom_jitter(size = 2.5,height = 0, width = 0.1,
              aes(color = cluster)) +
  scale_color_manual(values = colors,
                     name = "KNN clusters") +
  geom_point(data = res_ZIP.sum,
             size = 3.5,
             color = 'black',
             aes(x = Category, y = Mean)) +
  geom_errorbar(data = res_ZIP.sum,
                aes(x = Category, y = Mean,
                    ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.03) +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[1], xend = 1, yend = res_ZIP.sum$Mean[1]),
               linetype="dashed", colour = 'grey50') +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[2], xend = 2, yend = res_ZIP.sum$Mean[2]),
               linetype="dashed", colour = 'grey50') +
  # pval annotation
  geom_segment(aes(x = 1, xend = 2, y = 16.5, yend = 16.5)) +
  annotate("text", x = 1.5, y = 18, label = paste('P-value: ',round(pval, 5))) +
  # annotate("text", x = 1.5, y = 18, label = paste('P-value < ','0.0001')) +
  scale_fill_manual(values = c('#F5EC49',
                               '#3D9CE6'),
                    labels = c('Nucleotide\nantimetabolite',
                               'Other')) + 
  labs(x = 'Drug category',
       y = 'ZIP score') +
  ggrepel::geom_text_repel(aes(label = plot_labs),
                           box.padding = 0.5, 
                           nudge_x = -0.4,
                           max.overlaps = Inf) +
  scale_x_discrete(name = '',
                   labels = c("Nucleotides" = "Nucleotide \n antimetabolite","Other" = "Other")) +
  theme_cowplot(15, font_family = "Arial") +
  theme(axis.text.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.text.y = element_text(face = "bold", size = 13, color = 'black'))

ggsave(here('exploration', 'violin_nucleotides_biolog_PAPER.pdf'), height = 8, width = 9)




#### Most Syn scores stats ####

# homoscedasticity test
CG = subset(res_ZIP_cats, Category == "Nucleotides")$Most.synergistic.area.score
TG = subset(res_ZIP_cats, Category != "Nucleotides")$Most.synergistic.area.score
tidy(var.test(CG, TG))

# as they suffer from heterocedasticity, var.equal to FALSE

nucl_stats = t.test(CG, TG, var.equal = T)

wilcox.test(Most.synergistic.area.score ~ Category, data = res_ZIP_cats)
lm(Most.synergistic.area.score ~ Category, data = res_ZIP_cats) %>% summary


nucl_stats %>% tidy %>% write_csv("exploration/biolog_most_ZIP_stats_PAPER.csv")

library(broom)
pval = tidy(nucl_stats)$p.value

res_ZIP.sum = res_ZIP_cats %>%
  group_by(Category) %>%
  summarise(Mean = mean(Most.synergistic.area.score),
            SD = sd(Most.synergistic.area.score),
            SEM = SD/sqrt(n()))



clusters = read_csv('exploration/Multivariate_final/biolog_with_clusters.csv')

# clusters = results_cluster %>% distinct(Drug, .keep_all = T)


library(randomcoloR)

colors = randomColor(7)

colors = c("#ed2863", "#87ffe3", "#cd4cce", "#a57e11" ,"#fc33b2", "#eac885", "#1b6bea")

## MAIN PLOT
res_ZIP_cats %>%
  left_join(clusters) %>% 
  drop_na(cluster) %>% 
  mutate(cluster = factor(cluster),
         plot_labs = case_when(Category == "Nucleotides" ~ Drug)) %>% 
  ggplot(aes(x = Category, y = Most.synergistic.area.score, fill = Category)) +
  geom_violin() +
  geom_jitter(size = 2.5,height = 0, width = 0.1,
              aes(color = cluster)) +
  scale_color_manual(values = colors,
                     name = "KNN clusters") +
  geom_point(data = res_ZIP.sum,
             size = 3.5,
             color = 'black',
             aes(x = Category, y = Mean)) +
  geom_errorbar(data = res_ZIP.sum,
                aes(x = Category, y = Mean,
                    ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.03) +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[1], xend = 1, yend = res_ZIP.sum$Mean[1]),
               linetype="dashed", colour = 'grey50') +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[2], xend = 2, yend = res_ZIP.sum$Mean[2]),
               linetype="dashed", colour = 'grey50') +
  # pval annotation
  geom_segment(aes(x = 1, xend = 2, y = 16.5, yend = 16.5)) +
  annotate("text", x = 1.5, y = 18, label = paste('P-value: ',round(pval, 3))) +
  # annotate("text", x = 1.5, y = 18, label = paste('P-value < ','0.0001')) +
  scale_fill_manual(values = c('#D1D3D4',
                               '#E6E7E8'),
                    labels = c('Nucleotide\nantimetabolite',
                               'Other')) + 
  labs(x = 'Drug category',
       y = 'Most Synergistic ZIP area') +
  # ggrepel::geom_text_repel(aes(label = plot_labs),
  #                          box.padding = 0.5, 
  #                          nudge_x = -0.4,
  #                          max.overlaps = Inf) +
  scale_x_discrete(name = '',
                   labels = c("Nucleotides" = "Nucleotide \n antimetabolite","Other" = "Other")) +
  theme_cowplot(15, font_family = "Arial") +
  theme(axis.text.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.text.y = element_text(face = "bold", size = 13, color = 'black'))

ggsave(here('exploration', 'violin_nucleotides_biolog_PAPER_mostSyn.pdf'), height = 8, width = 9)

# stats between clusters ------

res_ZIP_cats %>%
  left_join(clusters) %>% 
  drop_na(cluster) %>% 
  mutate(cluster = factor(cluster)) %>% 
  # group_by(cluster) %>% 
  rstatix::pairwise_t_test(Synergy.score ~ cluster,
                           p.adjust.method = "fdr",
                           detailed = TRUE) %>% 
  arrange(p) %>% 
  write_csv('exploration/Multivariate_final/stats_synergy_score_PW.csv')


res_ZIP_cats %>%
  left_join(clusters) %>% 
  drop_na(cluster) %>% 
  mutate(cluster = factor(cluster)) %>% 
  # group_by(cluster) %>% 
  rstatix::pairwise_t_test(Most.synergistic.area.score ~ cluster, 
                           p.adjust.method = "fdr",
                           detailed = TRUE) %>% 
  arrange(p) %>% 
  write_csv('exploration/Multivariate_final/stats_Most_synergy_score_PW.csv')



## boxplots with cluster information -------------

res_ZIP_cluster = res_ZIP_cats %>%
  left_join(clusters) %>% 
  drop_na(cluster) %>% 
  mutate(cluster = factor(cluster))

res_ZIP_cluster %>% 
  bind_rows(res_ZIP_cluster %>% mutate(cluster = 'ALL')) %>% 
  mutate(cluster = factor(cluster, 
                          levels = c("ALL", 1,2,3,4,5,6,7))) %>% 
  ggplot(aes(y = Most.synergistic.area.score, 
             x = cluster, fill = cluster)) +
  geom_boxplot(show.legend = F) +
  scale_fill_manual(values = c("#aaaaaa" , colors),
                    name = NULL) +
  geom_point(position = position_jitterdodge(), 
             size = 2,
             show.legend = F) +
  labs(x = "Cluster",
       y = 'Most Synergistic area score')

ggsave(here('exploration/Multivariate_final', 
            'boxplots_most_synergistic_paper.pdf'), 
       height = 5, width = 8)




# t-SNE 2D biolog drugs ---------------------------------------------------




### t-SNE ####
# perplexity = 20, theta = 0.5, dims = 2

tsne_fit = all_compounds %>% 
  filter(DB == 'Biolog') %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  Rtsne(check_duplicates = FALSE, pca = T, num_threads = 8,
        normalize = FALSE, max_iter = 5000, 
        perplexity = 20, theta = 0.7, dims = 2)

# generate data frame from tnse results
tsne.df = data.frame(tsne_fit$Y)
colnames(tsne.df) = c('Dim1', 'Dim2')

tsne.df = tsne.df %>% 
  as_tibble() %>% 
  mutate(Drug = all_compounds %>% filter(DB == 'Biolog') %>% pull(Drug)) %>% 
  left_join(clusters %>% select(Drug, cluster)) %>% 
  mutate(cluster = as.factor(cluster))
  


tsne.df %>%
  ggplot(aes(Dim1, Dim2, color = cluster)) +
  geom_point(size = 4) +
  scale_color_manual(values = colors,
                    name = NULL)




ggsave("exploration/Multivariate_final/2D_tSNE_biolog.pdf", 
       height = 7, width = 9)
tsne.df %>% 
  write_csv("exploration/Multivariate_final/2D_tSNE_biolog.csv")




# 3D version 

t1 = list(family = 'Arial', size = 5)

clusters %>%
  mutate(cluster = as.factor(cluster)) %>% 
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~cluster, 
          colors = colors) %>% 
  add_markers() %>% 
  layout(
    font = t1,
    scene = list(
      yaxis = list(gridwidth = 1.5,
                   title = ""),
      xaxis = list(gridwidth = 1.5,
                   title = ""),
      zaxis = list(gridwidth = 1.5,
                   title = "")
    )
  ) 








