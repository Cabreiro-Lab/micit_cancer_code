
# libraries ---------------------------------------------------------------



library(readr)
library(multcomp)
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(openxlsx)
# library(PFun)
library(broom)
# library(multidplyr)
library(tau)
library(fmsb)
library(readxl)
library(cowplot)
library(rstatix)
library(extrafont)


# theme_set(theme_classic())

theme_set(theme_cowplot(15))

# load the data -----------------------------------------------------------

data = read_delim("stats_perseus_full.tsv", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)

names(data)

names(data) = gsub(' ', '_', names(data))

data




### TABLE TIDYING
# create groups depending on their DE profiles

data = data %>% 
  mutate(FU_C = case_when(`Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` < 0 ~ -1,
                          `Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` > 0 ~ 1,
                          TRUE ~ 0),
         micit1mM_C = case_when(`Student's_T-test_Significant_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_1mM_Micit_Control` < 0 ~ -1,
                                `Student's_T-test_Significant_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_1mM_Micit_Control` > 0 ~ 1,
                                 TRUE ~ 0),
         micit10mM_C = case_when(`Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` < 0 ~ -1,
                                 `Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` > 0 ~ 1,
                                  TRUE ~ 0),
         micit1mM_5FU_C = case_when(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` < 0 ~ -1,
                                    `Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` > 0 ~ 1,
                                     TRUE ~ 0),
         micit10mM_5FU_C = case_when(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` < 0 ~ -1,
                                     `Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` > 0 ~ 1,
                                     TRUE ~ 0),
         .after = "KEGG_name") 


# let's rename columns and remove the ones that are useless

names(data)

data_tidy = data %>% 
  select(Protein_IDs:Gene_names, Control:KEGG_name, -`Q-value`,
         mol_weight = `Mol._weight_[kDa]`,
         # 5FU_control
         logPval_FU_C = `-Log_Student's_T-test_p-value_5FU_Control`,
         Qval_FU_C = `Student's_T-test_q-value_5FU_Control`,
         difference_FU_C = `Student's_T-test_Difference_5FU_Control`,
         # 1 mM micit
         logPval_micit1mM_C = `-Log_Student's_T-test_p-value_1mM_Micit_Control`,
         Qval_micit1mM_C = `Student's_T-test_q-value_1mM_Micit_Control`,
         difference_micit1mM_C = `Student's_T-test_Difference_1mM_Micit_Control`,
         # 10 mM micit
         logPval_micit10mM_C = `-Log_Student's_T-test_p-value_10mM_Micit_Control`,
         Qval_micit10mM_C = `Student's_T-test_q-value_10mM_Micit_Control`,
         difference_micit10mM_C = `Student's_T-test_Difference_10mM_Micit_Control`,
         # 5FU + 1 mM micit
         logPval_micit1mM_5FU_C = `-Log_Student's_T-test_p-value_5FU_1mM_Micit_Control`,
         Qval_micit1mM_5FU_C = `Student's_T-test_Test_statistic_5FU_1mM_Micit_Control`,
         difference_micit1mM_5FU_C = `Student's_T-test_Difference_5FU_1mM_Micit_Control`,
         # 5FU + 10 mM micit
         logPval_micit10mM_5FU_C = `-Log_Student's_T-test_p-value_5FU_10mM_Micit_Control`,
         Qval_micit10mM_5FU_C = `Student's_T-test_q-value_5FU_10mM_Micit_Control`,
         difference_micit10mM_5FU_C = `Student's_T-test_Difference_5FU_10mM_Micit_Control`
  )


write_csv(data_tidy, here('summary','Stats_with_differences_full.csv'))




# convert data to long format

names(data)

data_long = data_tidy %>% 
  pivot_longer(Control:`5FU_10mM_Micit_4`, names_to = 'Sample', values_to = 'Intensity') 

# fix variable names
data_long = data_long %>% select(Gene_names, Sample, Intensity, everything()) %>% 
  mutate(Sample = case_when(Sample == 'Control_1' | Sample == 'Control_2' | Sample == 'Control_3' | Sample == 'Control_4' ~ 'Control',
                            Sample == '5FU_1' | Sample == '5FU_2' | Sample == '5FU_3' | Sample == '5FU_4' ~ '5FU',
                            Sample == '1mM_Micit_1' | Sample == '1mM_Micit_2' | Sample == '1mM_Micit_3' | Sample == '1mM_Micit_4' ~ '1mM_Micit',
                            Sample == '5FU_1mM_Micit_1' | Sample == '5FU_1mM_Micit_2' | Sample == '5FU_1mM_Micit_3' | Sample == '5FU_1mM_Micit_4' ~ '5FU_1mM_Micit',
                            Sample == '10mM_Micit_1' | Sample == '10mM_Micit_2' | Sample == '10mM_Micit_3' | Sample == '10mM_Micit_4' ~ '10mM_Micit',
                            Sample == '5FU_10mM_Micit_1' | Sample == '5FU_10mM_Micit_2' | Sample == '5FU_10mM_Micit_3' | Sample == '5FU_10mM_Micit_4' ~ '5FU_10mM_Micit',
                            TRUE ~ Sample),
         Sample = factor(Sample)) %>% 
  arrange(desc(Sample))

data_long

data_long %>% 
  select(Gene_names:Protein_IDs, KEGG_name) %>% 
  group_by(Gene_names) %>% 
  mutate(Replicate = 1:n(), .before = Protein_IDs) %>% 
  write_csv(here('summary','data_long.csv'))




data_long %>% 
  select(Gene_names, KEGG_name, GOBP_name, GOMF_name, GOCC_name) %>% 
  distinct(Gene_names,.keep_all = T) %>% 
  separate_rows(KEGG_name, sep = ';') %>% 
  write_csv(here('summary', 'protein_categories.csv'))


# UpSet general -----------------------------------------------------------

library(UpSetR)


### global differences ####

data

FU_genes = data %>% 
  filter(`Student's_T-test_Significant_5FU_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_genes = data %>% 
  filter(`Student's_T-test_Significant_10mM_Micit_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_1mM_5FU_genes = data %>% 
  filter(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_5FU_genes = data %>% 
  filter(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character




list_full = list(
  '5FU' = FU_genes,
  'Micit 10mM' = micit_10mM_genes,
  '5FU+Micit 1mM' = micit_1mM_5FU_genes,
  '5FU+Micit 10mM' = micit_10mM_5FU_genes
)


# generate a combination matrix
m1 = make_comb_mat(list_full)

m2 = make_comb_mat(list_full, mode = "intersect")

# size of the groups
comb_size(m1)
comb_size(m2)


UpSet(m1)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_full_distinct.pdf'),
             width = 8, height = 6, useDingbats = FALSE)


UpSet(m2)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_full_intersect.pdf'),
             width = 8, height = 6, useDingbats = FALSE)







# upset up/down -----------------------------------------------------------


# 5FU
FU_genes_down = data %>% 
  filter(`Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

FU_genes_up = data %>% 
  filter(`Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character


# micit 10 mM
micit_10mM_genes_down = data %>% 
  filter(`Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_genes_up = data %>% 
  filter(`Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character


# 5FU + 1mM micit
micit_1mM_5FU_genes_down = data %>% 
  filter(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_1mM_5FU_genes_up = data %>% 
  filter(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character


# 5FU + 10mM micit
micit_10mM_5FU_genes_down = data %>% 
  filter(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_5FU_genes_up = data %>% 
  filter(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character




list_directions = list(
  '5FU down' = FU_genes_down,
  '5FU up' = FU_genes_up,
  'Micit 10 mM down' = micit_10mM_genes_down,
  'Micit 10 mM up' = micit_10mM_genes_up,
  '5FU + Micit 1 mM down' = micit_1mM_5FU_genes_down,
  '5FU + Micit 1 mM up' = micit_1mM_5FU_genes_up,
  '5FU + Micit 10 mM down' = micit_10mM_5FU_genes_down,
  '5FU + Micit 10 mM up' = micit_10mM_5FU_genes_up
)


# generate a combination matrix
m1 = make_comb_mat(list_directions)

m2 = make_comb_mat(list_directions, mode = "intersect")




UpSet(m1)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_directions_distinct.pdf'),
             width = 8, height = 6, useDingbats = FALSE)


UpSet(m2)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_directions_intersect.pdf'),
             width = 8, height = 6, useDingbats = FALSE)
  




## extracting gene sets ####-----------------------------------------------

# gene order:
# 5FU || 10mM Micit || 5FU+1mM Micit || 5FU+10mM Micit
# down/up

comb_name(m2, readable = T)

# 5FU ∩ 10 mM Micit
FU_Micit1_down = extract_comb(m2, "10100000")
FU_Micit1_up = extract_comb(m2, "01010000")

# 5FU ∩ 5FU+1mM Micit
FU_Micit1FU_down = extract_comb(m2, "10001000")
FU_Micit1FU_up = extract_comb(m2, "01000100")

# 5FU ∩ 5FU+10mM Micit
FU_Micit10FU_down = extract_comb(m2, "10000010")
FU_Micit10FU_up = extract_comb(m2, "01000001")

# 10mM micit ∩ 5FU+10mM Micit
Micit10_Micit10FU_down = extract_comb(m2, "00100010")
Micit10_Micit10FU_up = extract_comb(m2, "00010001")

## Specific combinations
# 5FU+micit 10mM ∩ 5FU+micit 1mM
Micit10FU_Micit1FU_down = extract_comb(m2, "00001010")
Micit10FU_Micit1FU_up = extract_comb(m2, "00000101")

# 5FU+10mM micit ∩ 5FU+1mM micit U 10mM Micit
trio_down = extract_comb(m2, "00101010")
trio_up = extract_comb(m2, "00010101")


# 5FU + 1mM ! 5FU


setdiff(micit_1mM_5FU_genes_up, FU_genes_up)
length(setdiff(micit_1mM_5FU_genes_up, FU_genes_up))
length(micit_1mM_5FU_genes_up)
length(FU_genes_up)

setdiff(micit_1mM_5FU_genes_down, FU_genes_down)
length(setdiff(micit_1mM_5FU_genes_down, FU_genes_down))
length(micit_1mM_5FU_genes_down)
length(FU_genes_down)



FU_Micit1_diff_down = setdiff(micit_1mM_5FU_genes_down, FU_genes_down)
FU_Micit1_diff_up = setdiff(micit_1mM_5FU_genes_up, FU_genes_up)




# make a file of genes in each contrast
explanations = c('5FU ∩ 10 mM Micit', '5FU ∩ 5FU+1mM Micit',
                 ' 5FU ∩ 5FU+10mM Micit', '10mM micit ∩ 5FU+10mM Micit',
                 '5FU+micit 10mM ∩ 5FU+micit 1mM', 
                 '5FU+10mM micit ∩ 5FU+1mM micit ∩ 10mM Micit',
                 '5FU + 1mM micit ! 5FU')

table_names = c('FU_Micit1', 'FU_Micit1FU', 'FU_Micit10FU', 'Micit10_Micit10FU',
                'Micit10FU_Micit1FU', 'trio', 'FU_Micit1_diff')

exp_df = tibble(explanations, table_names)

library(openxlsx)

list_of_tables = list(
  'Metadata' = exp_df,
  FU_Micit1_down = unlist(str_split(FU_Micit1_down, pattern = ';')),
  FU_Micit1_up = unlist(str_split(FU_Micit1_up, pattern = ';')),
  FU_Micit1FU_down = unlist(str_split(FU_Micit1FU_down, pattern = ';')),
  FU_Micit1FU_up = unlist(str_split(FU_Micit1FU_up, pattern = ';')),
  FU_Micit10FU_down = unlist(str_split(FU_Micit10FU_down, pattern = ';')),
  FU_Micit10FU_up = unlist(str_split(FU_Micit10FU_up, pattern = ';')),
  Micit10_Micit10FU_down = unlist(str_split(Micit10_Micit10FU_down, pattern = ';')),
  Micit10_Micit10FU_up = unlist(str_split(Micit10_Micit10FU_up, pattern = ';')),
  Micit10FU_Micit1FU_down = unlist(str_split(Micit10FU_Micit1FU_down, pattern = ';')),
  Micit10FU_Micit1FU_up =unlist(str_split(Micit10FU_Micit1FU_up, pattern = ';')),
  trio_down = unlist(str_split(trio_down, pattern = ';')),
  trio_up = unlist(str_split(trio_up, pattern = ';')),
  FU_Micit1_diff_down = unlist(str_split(FU_Micit1_diff_down, pattern = ';')),
  FU_Micit1_diff_up = unlist(str_split(FU_Micit1_diff_up, pattern = ';'))
)

write.xlsx(list_of_tables, here('summary','Protein_sets.xlsx'))


for (i in 2:length(list_of_tables)) {
  file_name = names(list_of_tables)[i]
  write.table(list_of_tables[[i]],here('summary/protein_groups_intersection', 
                                     paste0(file_name,'.txt')),
              quote=F,col.names = F,row.names = F)
}



# saving the original protein differences
list_of_tables = list(
  FU_DOWN = unlist(str_split(FU_genes_down, pattern = ';')),
  FU_UP = unlist(str_split(FU_genes_up, pattern = ';')),
  micit10_DOWN = unlist(str_split(micit_10mM_genes_down, pattern = ';')),
  micit10_UP = unlist(str_split(micit_10mM_genes_up, pattern = ';')),
  micit1_5FU_DOWN = unlist(str_split(micit_1mM_5FU_genes_down, pattern = ';')),
  micit1_5FU_UP = unlist(str_split(micit_1mM_5FU_genes_up, pattern = ';')),
  micit10_5FU_DOWN = unlist(str_split(micit_10mM_5FU_genes_down, pattern = ';')),
  micit10_5FU_UP = unlist(str_split(micit_10mM_5FU_genes_up, pattern = ';'))
)

write.xlsx(list_of_tables, here('summary','Proteins_directions_MainGroups.xlsx'))

# loop to save everything in a folder
for (i in 1:length(list_of_tables)) {
  file_name = names(list_of_tables)[i]
  write.table(list_of_tables[[i]],here('summary/protein_groups', 
                                       paste0(file_name,'.txt')),
              quote=F,col.names = F,row.names = F)
}



# plot examples of genes
gene = c('RPL5','EIF5A;EIF5AL1')
data_long %>% filter(Gene_names %in%  gene) %>% 
  # filter(Sample != 'Control') %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  labs(title = gene) +
  facet_wrap(~Gene_names)





# radarplot ---------------------------------------------------------------


### LOAD DATA
# 1mM micit + 5FU
micit_1mM_5FU_enrich_down = read_excel("summary/protein_groups/output/micit_1mM_5FU_genes_down_output.xlsx", 
                                              sheet = "Process") %>% 
  select(-`...1`)

micit_1mM_5FU_enrich_up = read_excel("summary/protein_groups/output/micit_1mM_5FU_genes_up_output.xlsx", 
                                       sheet = "Process") %>% 
  select(-`...1`)

# 10mM micit + 5FU
micit_10mM_5FU_enrich_down = read_excel("summary/protein_groups/output/micit_10mM_5FU_genes_down_output.xlsx", 
                                               sheet = "Process")%>% 
  select(-`...1`)

micit_10mM_5FU_enrich_up = read_excel("summary/protein_groups/output/micit_10mM_5FU_genes_up_output.xlsx", 
                                        sheet = "Process")%>% 
  select(-`...1`)

# 10mM micit
micit_10mM_enrich_down = read_excel("summary/protein_groups/output/micit_10mM_genes_down_output.xlsx", 
                                        sheet = "Process")%>% 
  select(-`...1`)

micit_10mM_enrich_up = read_excel("summary/protein_groups/output/micit_10mM_genes_up_output.xlsx", 
                                      sheet = "Process")%>% 
  select(-`...1`)

# 5FU
FU_enrich_down = read_excel("summary/protein_groups/output/FU_genes_down_output.xlsx", 
                                   sheet = "Process") %>% 
  select(-`...1`)

FU_enrich_up = read_excel("summary/protein_groups/output/FU_genes_up_output.xlsx", 
                             sheet = "Process") %>% 
  select(-`...1`)


### FUNCTIONS 
terms_in_words = function(terms, elms=10){
  word_count = textcnt(terms,  split = ' ',method = "string", 
                       n = 1L)
  
  df = data.frame(matrix(word_count)) %>%
    mutate(words = names(word_count)) 
  
  names(df) = c('count', 'words')
  
  df = tibble(df)
  
  df = df %>% arrange(desc(count)) %>% 
    filter(!words %in% c('to', 'of', 'activity','process', 'metabolic', 
                         'regulation', 'cellular', 'cell', 'compound',
                         'in','to','from','g1/s','in.')) %>% 
    head(elms) %>% 
    mutate(count = count/sum(count))
  
  return(df)
}


down = terms_in_words(micit_1mM_5FU_enrich_down$description) %>% mutate(direction = 'down')
up = terms_in_words(micit_1mM_5FU_enrich_up$description) %>% mutate(direction = 'up')


radar_plot = function(down,up){
  test_radar = down %>% bind_rows(up)
  test_radar = test_radar %>% 
    pivot_wider(names_from = words, values_from = count) %>% 
    replace(is.na(.), 0) %>% 
    data.frame
  
  rownames(test_radar) = test_radar[,1]
  test_radar[,1] = NULL
  
  test_radar['Max',] = max(test_radar)
  test_radar['Min',] = 0
  
  test_radar = test_radar[c(3,4,1,2),]
  
  radarchart(test_radar)
  op <- par(mar = c(1, 2, 2, 1))
  create_beautiful_radarchart(test_radar, caxislabels = c(0, 0.25, 0.5, 0.75, 1),
                              color = c("#00AFBB", "#E7B800"))
  legend(
    x = "bottom", legend = rownames(test_radar[-c(1,2),]), horiz = TRUE,
    bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800"),
    text.col = "black", cex = 1, pt.cex = 1.5
  )
  par(op)
  
}


create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


### PLOT RESULTS
# 1mM micit + 5FU
down = terms_in_words(micit_1mM_5FU_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(micit_1mM_5FU_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'micit1mM_5FU_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)

# 10mM micit + 5FU
down = terms_in_words(micit_10mM_5FU_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(micit_10mM_5FU_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'micit10mM_5FU_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)

# 10mM micit
down = terms_in_words(micit_10mM_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(micit_10mM_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'micit10mM_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)


# 5FU
down = terms_in_words(FU_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(FU_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', '5FU_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)





# t-test with interactions ------------------------------------------------


# first do a test

test = data_long %>% 
  filter(Gene_names == 'RPL5') %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit')))

model1 = lm(Intensity ~ 0 +Sample, data = test)
summary(model1)

test %>% 
  group_by(Sample) %>% 
  summarise(Mean = mean(Intensity))

means = test %>% 
  group_by(Sample) %>% 
  summarise(Mean = mean(Intensity))

K = matrix(c(-1,   # Control
             1,  # 5FU
             0,  # 1mM_Micit
             0,  # 10mM_Micit
             0,  # 5FU_1mM_Micit
             0),  # 5FU_10mM_Micit
           1)


t <- glht(model1, linfct = mcp(Sample = K), test = adjusted('none'))
tidy(summary(t, test = adjusted('none')))


# plot examples of genes
gene = c('RRM2')
data_long %>% filter(Gene_names %in%  gene) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  # filter(Sample != 'Control') %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  labs(title = gene) +
  facet_wrap(~Gene_names)




### start building our test
# first let's remove things that have n <= 2

data_long %>% 
  group_by(Gene_names) %>% 
  drop_na(Intensity, Gene_names) %>% 
  count(Sample) %>% 
  arrange(desc(n)) %>% 
  ggplot(aes(n)) +
  geom_histogram()



## FILTER LIST
# get the list of genes with only 1 datapoint or less in Control 

removals = data_long %>% 
  group_by(Gene_names) %>% 
  drop_na(Intensity, Gene_names) %>% 
  count(Sample) %>% 
  filter(n < 2) %>% 
  distinct(Gene_names) %>% t %>%  as.character

# removing NA groups
removals_na = data_long %>% group_by(Gene_names, Sample) %>% 
  summarise(na_count = sum(is.na(Intensity))) %>% 
  arrange(desc(na_count)) %>% 
  filter(na_count > 3) %>% distinct(Gene_names) %>% t %>% as.character

# final list of removals
removals = unique(c(removals,removals_na))

## read contrasts
library(readxl)
Contrasts = read_excel("Contrasts.xlsx")

contr_mat = as.matrix(Contrasts[,5:10])

contr_factors = Contrasts[,1:4]


row.names(contr_mat) = Contrasts$contrast


statsR_raw = data_long %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  filter(!(Gene_names %in% removals)) %>% 
  drop_na(Intensity) %>% 
  group_by(Gene_names) %>% 
  nest() %>% 
  mutate(model = map(data, lm, formula = 'Intensity ~ 0 + Sample'),
         inter = map(model, glht, linfct = contr_mat, test = adjusted('none')),
         inter_sum = map(inter, summary, test = adjusted('none')),       ## this way I can extract raw p-values
         inter_tidy = map(inter_sum, tidy))



# summary(statsR_raw$inter[[1]])
# statsR_raw$inter_sum[[1]]



statsR = statsR_raw %>% 
  select(Gene_names, inter_tidy) %>% 
  unnest(cols = c(inter_tidy)) %>% 
  mutate(p.stars = gtools::stars.pval(p.value)) %>% 
  group_by(contrast) %>% 
  mutate(FDR = p.adjust(p.value, method = 'fdr')) %>% 
  ungroup %>% 
  mutate(FDR_stars = gtools::stars.pval(FDR))


statsR = statsR %>% left_join(contr_factors) %>% 
  select(Gene_names, Contrast_type:Target,contrast,everything())

write.xlsx(statsR %>% arrange(Gene_names), here('summary','Stats_R_interactions.xlsx'))





### protein groups from R stats ####

# this xlsx will be fed to the python script to extract categories and so on
# save protein groups
FU_UP = statsR %>% filter(contrast == '5FU-Ctrl', FDR < 0.05, estimate > 0) %>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
FU_DOWN = statsR %>% filter(contrast == '5FU-Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')

micit10_UP = statsR %>% filter(contrast == '10mMmicit - Ctrl', FDR < 0.05, estimate > 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
micit10_DOWN = statsR %>% filter(contrast == '10mMmicit - Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names)%>% separate_rows(genes, sep = ';')

micit1_5FU_UP = statsR %>% filter(contrast == '5FU_1mM_Micit - Ctrl', FDR < 0.05, estimate > 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
micit1_5FU_DOWN = statsR %>% filter(contrast == '5FU_1mM_Micit - Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names)%>% separate_rows(genes, sep = ';')

micit10_5FU_UP = statsR %>% filter(contrast == '5FU_10mM_Micit - Ctrl', FDR < 0.05, estimate > 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
micit10_5FU_DOWN = statsR %>% filter(contrast == '5FU_10mM_Micit - Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')

list_of_tables = list(
  FU_UP = FU_UP,
  FU_DOWN = FU_DOWN,
  micit10_UP = micit10_UP,
  micit10_DOWN = micit10_DOWN,
  micit1_5FU_UP = micit1_5FU_UP,
  micit1_5FU_DOWN = micit1_5FU_DOWN,
  micit10_5FU_UP = micit10_5FU_UP,
  micit10_5FU_DOWN = micit10_5FU_DOWN
)

write.xlsx(list_of_tables, here('summary','Protein_sets.xlsx'))



# exploratory plots -------------------------------------------------------



# plot examples of genes
gene = c('TYMS')
data_long %>% filter(Gene_names %in%  gene) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  # filter(Sample != 'Control') %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  # labs(title = gene) +
  facet_wrap(~Gene_names) +
  theme(axis.text.x = element_text(hjust = 1, angle=45))




gene = c('PTPRF','PDCD4','TYMS','RRM2','KIAA0101','CDKN1A')
for (gen in gene){
  data_long %>% filter(Gene_names %in%  gen) %>% 
    mutate(Sample = factor(Sample, levels = c('Control', 
                                              '5FU',
                                              '1mM_Micit',
                                              '10mM_Micit',
                                              '5FU_1mM_Micit', 
                                              '5FU_10mM_Micit'))) %>% 
    # filter(Sample != 'Control') %>% 
    ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    # labs(title = gene) +
    facet_wrap(~Gene_names) +
    theme(axis.text.x = element_text(hjust = 1, angle=45))
  ggsave(here('summary',paste0(gen,'_boxplot.pdf')), height = 6, width = 8)
}






# unique protein plots ----------------------------------------------------


stats_control = statsR %>% 
  filter(Contrast_type %in% c('Treatment', 
                              'Interaction'), 
         Reference %in% c('Control', 
                          '1mM_Micit,Control')) %>% 
  mutate(FDR_stars = case_when(Contrast_type == 'Interaction' ~ '',
                               TRUE ~ FDR_stars),
         estimate = case_when(Contrast_type == 'Interaction' ~ 0,
                            TRUE ~ estimate),
         Target = case_when(Contrast_type == 'Interaction' ~ 'Control',
                            TRUE ~ Target))

gene_plotter = function(gene) {
  
  data = data_long %>% 
    filter(Gene_names %in%  gene) %>% 
    left_join(stats_control %>% 
                filter(Gene_names %in% gene) %>% 
                rename(Sample = Target)) %>% 
    mutate(Sample = case_when(
      Sample == '1mM_Micit' ~ 'Micit\n(1mM)',
      Sample == '10mM_Micit' ~ 'Micit\n(10mM)',
      Sample == '5FU_1mM_Micit' ~ '5FU \n+\n Micit\n(1mM)',
      Sample == '5FU_10mM_Micit' ~ '5FU \n+\n Micit\n(10mM)',
      TRUE ~ Sample)) %>%
    mutate(Sample = factor(Sample, levels = c('Control',
                                              '5FU',
                                              'Micit\n(1mM)',
                                              'Micit\n(10mM)',
                                              '5FU \n+\n Micit\n(1mM)', 
                                              '5FU \n+\n Micit\n(10mM)'))) 
  
  minval = data %>% 
    drop_na(Intensity) %>% 
    mutate(minval = min(Intensity)) %>% 
    distinct(minval) %>% 
    pull(minval)
  
  maxval = data %>% 
    drop_na(Intensity) %>% 
    mutate(maxval = min(Intensity)) %>% 
    distinct(maxval) %>% 
    pull(maxval)
  
  
  data %>%
    ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
    geom_boxplot() +
    geom_text(aes(label = FDR_stars), y = maxval * 1.11) +
    geom_text(aes(label = paste0('log2FC:\n ',round(estimate,2)), 
                  y = maxval * 1.095)) +
    ylim(minval * 0.98, maxval * 1.12) +
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~Gene_names) 
  
  ggsave(here('summary/protein_individual_plots', 
              glue::glue('{gene}_plot.pdf')),
         height = 8, width = 10)
  
}

gene_plotter('TYMS')



gene_list = stats_control %>% 
  arrange(estimate) %>% 
  # head(50) %>% 
  distinct(Gene_names) %>%  pull(Gene_names)
  

for (gene in gene_list){
  gene_plotter(gene)
}

#loading packages

list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "tidyverse",
  "kableExtra"
)

for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

n.cores = parallel::detectCores() - 2

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

clusterEvalQ(cl,  library(tidyverse))

#check if it is registered (optional)
foreach::getDoParRegistered()

foreach(gene %in% gene_list) %dopar% {
  gene_plotter(gene)
}

parallel::stopCluster(cl = my.cluster)


# enrichment exploration --------------------------------------------------

library(readxl)
FU_enrich = read_excel("summary/Enrichment_RStats/FU/FU_output.xlsx", 
                        sheet = "Process") %>% 
  select(-`...1`)


FU_enrich %>%  
  mutate(strength = -log10(number_of_genes/number_of_genes_in_background)) %>% 
  filter(fdr <= 0.01)


FU_enrich = read_excel("summary/Enrichment_RStats/FU/FU_output.xlsx", 
                       sheet = "KEGG") %>% 
  select(-`...1`) %>% 
  mutate(Comparison = "5-FU vs Control")

micit10_enrich = read_excel("summary/Enrichment_RStats/micit10/micit10_output.xlsx", 
                       sheet = "KEGG") %>% 
  select(-`...1`) %>% 
  mutate(Comparison = "2-MiCit vs Control")

micit10_5FU_enrich = read_excel("summary/Enrichment_RStats/micit10_5FU/micit10_5FU_output.xlsx", 
                            sheet = "KEGG") %>% 
  select(-`...1`) %>% 
  mutate(Comparison = "2-MiCit + 5-FU vs Control")





kegg_string_enrich = FU_enrich %>% 
  bind_rows(micit10_enrich, micit10_5FU_enrich) %>% 
  select(Comparison, term, number_of_genes, number_of_genes_in_background, 
         p_value, fdr, description, direction)


kegg_string_enrich %>% write_csv("summary/FIGURES_paper/proteomics_KEGG_enrichment.csv")


# Heatmaps ----------------------------------------------------------------

# this specifies the groups to filter in the Rstats
contr_groups = c('5FU - 1mM_Micit', '5FU-Ctrl','10mMmicit - Ctrl','1mMmicit - Ctrl',
                 '5FU_1mM_Micit - Ctrl', '5FU_10mM_Micit - Ctrl')

### top 100 ####


large_effect_prots = statsR %>% 
  mutate(estimate = abs(estimate)) %>% 
  arrange(desc(estimate)) %>% 
  filter(FDR < 0.05) %>% 
  distinct(Gene_names, .keep_all = TRUE) %>% 
  head(100) %>% 
  separate_rows(Gene_names, sep = ';') %>% 
  select(Gene_names) %>% t %>% as.character()

write.table(large_effect_prots, 'large_effect_prots.txt',
            quote = F, row.names = F, col.names = F)

# KEGG pathways
data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% large_effect_prots) %>% 
  distinct(Gene_names,.keep_all = T) %>% 
  select(Gene_names, KEGG_name)  %>% 
  # separate(KEGG_name,sep=';',into=)
  separate_rows(KEGG_name, sep = ';') %>% 
  filter(!(KEGG_name %in% c('Phototransduction - fly','Cell cycle - yeast',
                          'Meiosis - yeast','Amoebiasis','Toxoplasmosis',
                          'Vascular smooth muscle contraction',
                          'Proximal tubule bicarbonate reclamation','Malaria'))) %>% 
  count(Gene_names) %>% arrange(desc(n))


heat_zscore_wide = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% large_effect_prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_top = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% large_effect_prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat_zscore_wide[,2:7])
row.names(heat_mat) = heat_zscore_wide$Gene_names

pval_mat = as.matrix(pval_top[,2:7])
row.names(pval_mat) = pval_top$Gene_names


Heatmap(heat_mat,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 5),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 4))


dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_top100.pdf'),
             width = 5, height = 9, useDingbats = FALSE)



### global heatmap ####


heat_zscore_wide_all = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  # filter(Gene_names %in% large_effect_prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)




heat_mat_all = as.matrix(heat_zscore_wide_all[,2:7])

row.names(heat_mat_all) = heat_zscore_wide_all$Gene_names


Heatmap(heat_mat_all,
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        na_col = "black",
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        show_row_names = FALSE)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_total.pdf'),
             width = 5, height = 6, useDingbats = FALSE)




### save dataset for comparison ####


comparison_df = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  pivot_wider(names_from = Sample, values_from = Mean)


colnames(comparison_df) = c("Gene_names","Micit_10",
                                   "Micit_1","FU","FU_10mM_Micit",
                                   "FU_1mM_Micit","Control")


write_csv(comparison_df,here('summary', 'means_samples.csv'))
write_csv(comparison_df,here('dataset_comparison', 'means_samples.csv'))




### TCA genes ####
tca = data %>% 
  filter(str_detect(KEGG_name,'TCA') ) %>% 
  select(Gene_names) %>% 
  t %>% as.character()


# calculate means and zscore and put it wider
heat_tca = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% tca) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_top = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% tca, 
         contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)


heat_mat = as.matrix(heat_tca[,2:7])
row.names(heat_mat) = heat_tca$Gene_names

pval_mat = as.matrix(pval_top[,2:7])
row.names(pval_mat) = pval_top$Gene_names


Heatmap(heat_mat,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 8))

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_TCA.pdf'),
             width = 5, height = 7, useDingbats = FALSE)



### mTOR genes ####

mtor = data %>% 
  filter(str_detect(KEGG_name,'mTOR') | str_detect(GOBP_name, 'mTOR') |  
           str_detect(GOMF_name, 'mTOR')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()



# calculate means and zscore and put it wider
heat_mtor = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% mtor) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% mtor, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)
  
# generate matrices
heat_mat = as.matrix(heat_mtor[,2:7])
row.names(heat_mat) = heat_mtor$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names


Heatmap(heat_mat, 
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_mTOR.pdf'),
             width = 5, height = 6, useDingbats = FALSE)







### p53 genes ####

colnames(heat_pres) = c('10mM', "1mM", '5-FU',
                        '5FU +\n10mM', '5FU +\n1mM', 
                        'Control')


heat_pres = heat_pres[,c(6,3,2,1,5,4)]

ha = HeatmapAnnotation(
  Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
              '5FU + Metabolite', '5FU + Metabolite'),
  col = list(Samples = c("Control" = "grey60", 
                         '5-FU' = 'red',
                         "Metabolite" = "orange", 
                         "5FU + Metabolite" = "blue")
  ),
  border = TRUE)


p53 = data %>% 
  filter(str_detect(KEGG_name,'p53') | str_detect(GOBP_name, 'p53') |  
           str_detect(GOMF_name, 'p53')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()


p53 =  data %>% 
  filter(str_detect(KEGG_name, 'p53 signaling pathway')) %>% 
  select(Gene_names) %>% pull


# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

heat_mat = heat_mat[,c('Control','5FU','1mM_Micit','10mM_Micit','5FU_1mM_Micit','5FU_10mM_Micit')]

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names
pval_mat = pval_mat[,c('Control','5FU','1mM_Micit','10mM_Micit','5FU_1mM_Micit','5FU_10mM_Micit')]

Heatmap(heat_mat, 
        name = "z-score",
        # column_km = 3,
        # row_km = 3,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 8),
        column_names_rot =30, 
        column_names_side = "top",
        # cell_fun = function(j, i, x, y, width, height, fill) {
        #   grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 5))},
        column_names_gp = gpar(fontsize = 6)
)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_p53.pdf'),
             width = 5, height = 4, useDingbats = FALSE)


heat_mat[107,]


### cell cycle genes ####

prots = data %>% 
  filter(str_detect(KEGG_name,'cell cycle') | str_detect(GOBP_name, 'cell cycle') |  
           str_detect(GOMF_name, 'cell cycle')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()



# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names


Heatmap(heat_mat, 
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 3),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        })

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_cell_cycle.pdf'),
             width = 5, height = 20, useDingbats = FALSE)





### mitochondrion genes ####

prots = data %>% 
  filter(str_detect(KEGG_name,'mitochondrion') | str_detect(GOBP_name, 'mitochondrion') |  
           str_detect(GOMF_name, 'mitochondrion')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()



# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names


Heatmap(heat_mat, 
        name = "z-score",
        # column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        })

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_mitochondrion.pdf'),
             width = 5, height = 20, useDingbats = FALSE)








### Pyr genes (log2FC) ####

prots = data %>% 
  filter(str_detect(KEGG_name,'Pyrimidine')) %>% 
  distinct(Gene_names) %>% pull(Gene_names)


prots = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  filter(FDR <= 0.05) %>% 
  distinct(Gene_names) %>%  pull(Gene_names)

# a further step of filtering (optional)
prots = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  filter(FDR <= 0.05) %>% 
  mutate(selected = case_when(Target == '10mM_Micit' & estimate < 0 ~ 'selected',
                              Target == '5FU' & estimate < 0 ~ 'selected',
                              TRUE ~ 'remove')) %>% 
  distinct(Gene_names) %>% pull(Gene_names)




heat = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups[2:6]) %>% 
  select(Gene_names, Target, estimate) %>% 
  pivot_wider(names_from = Target, values_from = estimate)


pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:6])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names

pval_mat = pval_mat[,1:5]

Heatmap(heat_mat, 
        name = "log2FC",
        # column_km = 3,
        # row_km = 3,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 10),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        }
        )

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_pyrimidine_log2FC.pdf'),
             width = 5, height = 10, useDingbats = FALSE)





### Pur genes (log2FC) ####

prots = data %>% 
  filter(str_detect(KEGG_name,'Purine')) %>% 
  distinct(Gene_names) %>% pull(Gene_names)


prots = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  filter(FDR <= 0.05) %>% 
  distinct(Gene_names) %>%  pull(Gene_names)



heat = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups[2:6]) %>% 
  select(Gene_names, Target, estimate) %>% 
  pivot_wider(names_from = Target, values_from = estimate)


pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:6])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names

pval_mat = pval_mat[,1:5]

Heatmap(heat_mat, 
        name = "log2FC",
        # column_km = 3,
        # row_km = 3,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 10),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        }
)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_purine_log2FC.pdf'),
             width = 5, height = 12, useDingbats = FALSE)



### mitochondrion genes (log2FC) ####


prots = data %>% 
  filter(str_detect(KEGG_name,'mitochondrion') | str_detect(GOBP_name, 'mitochondrion') |  
           str_detect(GOMF_name, 'mitochondrion'))  %>% 
  distinct(Gene_names) %>% pull(Gene_names)


prots = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups[2:6]) %>% 
  filter(FDR <= 0.05) %>% 
  distinct(Gene_names) %>%  pull(Gene_names)



heat = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups[2:6]) %>% 
  select(Gene_names, Target, estimate) %>% 
  pivot_wider(names_from = Target, values_from = estimate)


pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:6])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names

pval_mat = pval_mat[,1:5]

Heatmap(heat_mat, 
        name = "log2FC",
        # column_km = 3,
        # row_km = 3,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 10),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        }
)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_mitochondrion_log2FC.pdf'),
             width = 5, height = 17, useDingbats = FALSE)







### category selection genes ####

kegg_select = c('Cell cycle','p53 signaling pathway', 'Citrate cycle (TCA cycle)',
                'mTOR signaling pathway', 'RNA transport', 'Purine metabolism',
                'Pyrimidine metabolism', 'Ribosome', 'Fatty acid metabolism',
                'Focal adhesion', 'Glycolysis','Gap junction','ErbB signaling pathway',
                'MAPK signaling pathway')


prots = data %>% 
  separate_rows(KEGG_name, sep = ';') %>% 
  filter(KEGG_name %in% kegg_select) %>% 
  # select(Gene_names) %>% 
  distinct(Gene_names) %>%
  t %>% as.character()


# cluster = new_cluster(8)
# calculate means and zscore and put it wider
heat = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  # partition(cluster) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  # collect() %>%
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway)

# pval_mtor = statsR %>% 
#   filter(!(Gene_names %in% removals)) %>% 
#   filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
#   # create dummy variable for Control
#   mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
#                                TRUE ~ FDR_stars),
#          Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
#                             TRUE ~ Target)) %>% 
#   select(Gene_names,Target,FDR_stars) %>% 
#   pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,3:8])
row.names(heat_mat) = heat$Gene_names

# pval_mat = as.matrix(pval_mtor[,2:7])
# row.names(pval_mat) = pval_mtor$Gene_names

library(circlize)
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

ha = rowAnnotation(Pathway = heat$Pathway)

h1 = Heatmap(heat_mat, 
        name = "z-score",
        # column_km = 3,
        # row_km = 3,
        cluster_rows = FALSE, # turns off row clustering
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 4),
        right_annotation = ha,
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7))

heat$Pathway

# plot every pathway as a subheatmap
# LONG AND UGLY AND LONG AND UGLY CODE
# BUT IT WORKS, OK? 

cell_ht = Heatmap(heat_mat[1:76,], 
        name = "Z-score",
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 4),
        right_annotation = rowAnnotation(
          Pathway = heat$Pathway[1:76],
          show_annotation_name = F),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7))

tca_ht = Heatmap(heat_mat[77:103,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[77:103],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

erbB_ht = Heatmap(heat_mat[104:139,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[104:139],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

fa_ht = Heatmap(heat_mat[140:168,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[140:168],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

focal_ht = Heatmap(heat_mat[169:243,], 
                name = "Z-score",
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 4),
                right_annotation = rowAnnotation(
                  Pathway = heat$Pathway[169:243],
                  show_annotation_name = F),
                column_names_rot =30, 
                column_names_side = "top",
                column_names_gp = gpar(fontsize = 7))

gap_ht = Heatmap(heat_mat[244:280,], 
                   name = "Z-score",
                   show_row_names = FALSE,
                   row_names_gp = gpar(fontsize = 4),
                   right_annotation = rowAnnotation(
                     Pathway = heat$Pathway[244:280],
                     show_annotation_name = F),
                   column_names_rot =30, 
                   column_names_side = "top",
                   column_names_gp = gpar(fontsize = 7))

mapk_ht = Heatmap(heat_mat[281:361,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[281:361],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

mtor_ht = Heatmap(heat_mat[362:379,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[362:379],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

p53_ht = Heatmap(heat_mat[380:405,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[380:405],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

purine_ht = Heatmap(heat_mat[406:487,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[406:487],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

pyr_ht = Heatmap(heat_mat[488:555,], 
                    name = "Z-score",
                    show_row_names = FALSE,
                    row_names_gp = gpar(fontsize = 4),
                    right_annotation = rowAnnotation(
                      Pathway = heat$Pathway[488:555],
                      show_annotation_name = F),
                    column_names_rot =30, 
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 7))

ribosome_ht  = Heatmap(heat_mat[556:635,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[556:635],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

rnatransport_ht = Heatmap(heat_mat[636:757,], 
                      name = "Z-score",
                      show_row_names = FALSE,
                      row_names_gp = gpar(fontsize = 4),
                      right_annotation = rowAnnotation(
                        Pathway = heat$Pathway[636:757],
                        show_annotation_name = F),
                      column_names_rot =30, 
                      column_names_side = "top",
                      column_names_gp = gpar(fontsize = 7))



ht_list = cell_ht %v% tca_ht %v% erbB_ht %v% 
  fa_ht %v% focal_ht %v% gap_ht %v% mapk_ht %v% 
  mtor_ht %v% p53_ht %v% purine_ht %v% pyr_ht %v%  
  ribosome_ht %v% rnatransport_ht

lgd = Legend(labels = unique(heat$Pathway), title = "", 
             legend_gp = gpar(fill = 1:13),
             grid_height = unit(0, "cm"), grid_width = unit(0, "mm"),
             labels_gp = gpar(col = "white", fontsize = 0))

draw(ht_list, annotation_legend_list = lgd) 



dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_KEGG.pdf'),
             width = 8, height = 14, useDingbats = FALSE)






###  sig category selection genes ####

prots = data %>% 
  separate_rows(KEGG_name, sep = ';') %>% 
  filter(KEGG_name %in% kegg_select) %>% 
  # select(Gene_names) %>% 
  distinct(Gene_names) %>%
  t %>% as.character()

# take the proteins and see which ones are at least significant
# in one of the groups in respect to the control
prots_sig = statsR %>%
  filter(!(Gene_names %in% removals)) %>%
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  filter(FDR<=0.05) %>% 
  distinct(Gene_names) %>% t %>% as.character
  


# calculate means and zscore and put it wider
heat = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway)



# generate matrices
heat_mat = as.matrix(heat[,3:8])
row.names(heat_mat) = heat$Gene_names

# pval_mat = as.matrix(pval_mtor[,2:7])
# row.names(pval_mat) = pval_mtor$Gene_names


ha = rowAnnotation(Pathway = heat$Pathway)

h1 = Heatmap(heat_mat, 
             name = "z-score",
             # column_km = 3,
             # row_km = 3,
             cluster_rows = FALSE, # turns off row clustering
             show_row_names = FALSE,
             row_names_gp = gpar(fontsize = 4),
             right_annotation = ha,
             column_names_rot =30, 
             column_names_side = "top",
             column_names_gp = gpar(fontsize = 7))

h1


# get indexes to subset 
require(rlist)
cat_list = list()
for (elm in unique(heat$Pathway)){
  min_index = min(str_which(heat$Pathway,as.character(elm)))
  max_index = max(str_which(heat$Pathway,as.character(elm)))
  cat_list = list.append(cat_list, elm = c(min_index,max_index))
}
# fix a couple things
names(cat_list) = unique(heat$Pathway)
cat_list$`Citrate cycle (TCA cycle)` = c(39,64)


# this loop will create many variables for each pathway
# each one of them will be a heatmap
paths = unique(heat$Pathway)
list_of_vars = c()
for (i in 1:length(paths)){
  var_name = paste0(str_replace_all(paths[i],' ','_'),'_ht')
  print(var_name)
  list_of_vars = c(var_name,list_of_vars)
  
  min_index = cat_list[[i]][1]
  max_index = cat_list[[i]][2]
  print(min_index)
  print(max_index)
  
  var_name_heat = Heatmap(heat_mat[min_index:max_index,], 
                    name = "Z-score",
                    show_row_names = FALSE,
                    row_names_gp = gpar(fontsize = 4),
                    right_annotation = rowAnnotation(
                      Pathway = heat$Pathway[min_index:max_index],
                      show_annotation_name = F),
                    column_names_rot =30, 
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 7))
  
  assign(var_name, var_name_heat)
  
}

dummy = c()
for (i in list_of_vars){
  cosa = paste(i,'%v% ')
  dummy = c(cosa,dummy)
}
str_c(dummy, collapse = '')

ht_list = Cell_cycle_ht %v% `Citrate_cycle_(TCA_cycle)_ht` %v% 
  ErbB_signaling_pathway_ht %v% Fatty_acid_metabolism_ht %v% 
  Focal_adhesion_ht %v% Gap_junction_ht %v% MAPK_signaling_pathway_ht %v%
  mTOR_signaling_pathway_ht %v% 
  p53_signaling_pathway_ht %v% Purine_metabolism_ht %v% 
  Pyrimidine_metabolism_ht %v% Ribosome_ht %v% RNA_transport_ht

lgd = Legend(labels = unique(heat$Pathway), title = "", 
             legend_gp = gpar(fill = 1:13),
             grid_height = unit(0, "cm"), grid_width = unit(0, "mm"),
             labels_gp = gpar(col = "white", fontsize = 0))

draw(ht_list, annotation_legend_list = lgd) 



dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_KEGG_sigProts.pdf'),
             width = 8, height = 14, useDingbats = FALSE)






# radar plot zscores --------------------------------------------------------------

### barplots of selected pathways ####

zscores_paths = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway) %>% 
  group_by(Pathway, Sample) %>% 
  summarise(zMean = mean(z_score,na.rm = TRUE),
            zSD = sd(z_score, na.rm = TRUE),
            zSEM = zSD/sqrt(n())) %>% 
  ungroup




zscores_paths %>% 
  ungroup %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('Control', 
                                    '5FU',
                                    '1mM_Micit',
                                    '5FU_1mM_Micit',
                                    '10mM_Micit',
                                    '5FU_10mM_Micit'))) %>% 
  ggplot(aes(Pathway, zMean, fill = Sample)) +
  geom_histogram(stat='identity', position = 'dodge2') +
  geom_errorbar(aes(ymin = zMean - zSEM, ymax = zMean + zSEM), 
                position = position_dodge(0.9), width = 0) +
  scale_fill_manual(values = c('grey50',
                               '#FA2713', # 5FU
                               '#F0C600', # 1 micit
                               '#E05C10', # 5fu + 1 micit
                               '#F79E05', # 10 micit
                               '#AD4413'  # 5fu + 10 micit
  )) +
  labs(
    x = 'KEGG Pathway',
    y = 'z-score average (+- SEM)',
    fill = 'Condition'
  ) +
  facet_wrap(~Pathway, scales = 'free')


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'zScore_means_sigProts.pdf'),
             width = 10, height = 9, useDingbats = FALSE)




### calculate zscore mean for everything ####

# take the proteins and see which ones are at least significant
# in one of the groups in respect to the control
prots_sig = statsR %>%
  filter(!(Gene_names %in% removals)) %>%
  filter(contrast %in% contr_groups) %>% 
  filter(FDR<=0.05) %>% 
  distinct(Gene_names) %>% t %>% as.character


zscores_paths_all = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  # filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway) %>% 
  group_by(Pathway, Sample) %>% 
  summarise(zMean = mean(z_score,na.rm = TRUE),
            zSD = sd(z_score, na.rm = TRUE),
            zSEM = zSD/sqrt(n())) %>% 
  ungroup





zscores_paths_all %>% 
  ungroup %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('Control', 
                                    '5FU', 
                                    '1mM_Micit', '10mM_Micit',
                                    '5FU_1mM_Micit', '5FU_10mM_Micit'))) %>% 
  ggplot(aes(Pathway, zMean, fill = Sample)) +
  geom_histogram(stat='identity', position = 'dodge2') +
  geom_errorbar(aes(ymin = zMean - zSEM, ymax = zMean + zSEM), 
                position = position_dodge(0.9), width = 0) +
  scale_fill_manual(values = c('grey50',
                               '#FA2713', # 5FU
                               '#F0C600', # 1 micit
                               '#F79E05', # 10 micit
                               '#E05C10', # 5fu + 1 micit
                               '#AD4413'  # 5fu + 10 micit
  )) +
  labs(
    x = 'KEGG Pathway',
    y = 'z-score average (+- SEM)',
    fill = 'Condition'
  ) +
  facet_wrap(~Pathway, scales = 'free')


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'zScore_means_sigProts_ALL.pdf'),
             width = 70, height = 60, useDingbats = FALSE)







## create the spider plot #####


mins = c(min(zscores_paths$zMean))
maxs = c(max(zscores_paths$zMean))

# define values that are symmetric
mins = -1.5
maxs = 1.5

zscores_wide = zscores_paths %>% 
  select(Pathway:zMean) %>% 
  pivot_wider(names_from = Sample, values_from = zMean) %>%
  mutate(Max = maxs, Min = mins, Limit = 0, .before = `10mM_Micit`)

zscores_mat = as.matrix(zscores_wide[,2:10])

rownames(zscores_mat) = zscores_wide$Pathway

zscores_mat = t(zscores_mat)



zscores_mat = as.data.frame(zscores_mat[1:9,])


#~~~~~~~~~~~~~~~~~

## This piece of code will plot all Cell samples
# per separate in a single plot, might be good 
# for comparison

opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(0.8,4))
par(mfrow = c(3,2))
# Produce a radar-chart for each student
for (i in 4:nrow(zscores_mat)) {
  radarchart(
    zscores_mat[c(1:3, i), ],
    pfcol = c("#99999980",NA),
    pcol= c(NA,2), plty = 1, plwd = 2,
    title = row.names(zscores_mat)[i]
  )
}
# Restore the standard par() settings
par <- par(opar) 


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_ALL_samples.pdf'),
             width = 15, height = 10, useDingbats = FALSE)


#~~~~~~~~~~~~

select_zscores = zscores_mat[c(1:3,4,6,7),] 

names(select_zscores) = c("Cell cycle", "Citrate cycle\n(TCA cycle)", 
                          "ErbB signaling\npathway", "Fatty acid\nmetabolism", 
                          "Focal adhesion", "Gap junction", 
                          "MAPK signaling\npathway" ,
                          "mTOR signaling\npathway", "p53 signaling\npathway", 
                          "Purine\nmetabolism", "Pyrimidine\nmetabolism", 
                          "Ribosome", "RNA\ntransport" )

opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(1,4))

radarchart(
  select_zscores,
  caxislabels = c(-1.5, -0.75, 0, 0.75, 1.5),
  # vlcex = 1.1,
  pfcol = c("#99999980",NA,NA,NA),
  pcol= c(NA,
          '#F2AB0C',
          '#E60B1A',
          '#1700F2'), plty = 1, plwd = 2
)

# legend(
#   # x = "bottom",
#   x = -1,y = -1,
#   legend = rownames(select_zscores[-c(1,2,3),]), horiz = TRUE,
#   bty = "n", pch = 20 ,col= c('#F2AB0C',
#                           '#E60B1A',
#                           '#1700F2'),
#   text.col = "black", cex = 1, pt.cex = 2.5
# )
# 

row_nams = c('10mM Micit',
             '5-FU',
             '5-FU + Micit')
legend(
  # x = "bottom",
  x = 0.85,y = 1.1,
  legend = row_nams, horiz = F,
  bty = "n", pch = 20 ,col= c('#F2AB0C',
                              '#E60B1A',
                              '#1700F2'),
  text.col = "black", cex = 1, pt.cex = 2.5
)

par <- par(opar) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_5FU_Micit_comparison.pdf'),
             width = 7, height = 6.5, useDingbats = FALSE)




#~~~~~~~~~~~~~~
### ggradar version ####

library(ggradar)

zscores_wide = zscores_paths %>% 
  select(Pathway:zMean) %>% 
  pivot_wider(names_from = Sample, values_from = zMean) 

zscores_mat = as.matrix(zscores_wide[,2:7])

rownames(zscores_mat) = zscores_wide$Pathway

zscores_mat = t(zscores_mat)

zscores_mat = zscores_mat %>% as_tibble() %>% 
  mutate(Sample = unique(zscores_paths$Sample), .before = `Cell cycle`)


zscores_mat %>% 
  # filter(Pathway == 'Cell cycle') %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('Control', 
                                    '5FU', 
                                    '1mM_Micit', '10mM_Micit',
                                    '5FU_1mM_Micit', '5FU_10mM_Micit'))) %>% 
  ggradar(
    values.radar = c("-1.5", "0", "1.5"),
    grid.min = -1.5, grid.mid = 0, grid.max = 1.5,
    group.line.width = 1, 
    group.colours = c('grey50',
                      '#FA2713', # 5FU
                      '#F0C600', # 1 micit
                      '#F79E05', # 10 micit
                      '#E05C10', # 5fu + 1 micit
                      '#AD4413'),
    group.point.size = 3,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "bottom"
  )



names(zscores_mat)

zscores_mat %>% 
  add_row(Sample = 'DUMMY', `Cell cycle` = 0, `Citrate cycle (TCA cycle)` = 0,
          `ErbB signaling pathway` = 0, `Fatty acid metabolism` = 0,
          `Focal adhesion` = 0, `Gap junction`=0,`MAPK signaling pathway`=0,
          `mTOR signaling pathway`=0,`p53 signaling pathway`=0,
          `Purine metabolism`=0,`Pyrimidine metabolism`=0,
          Ribosome = 0, `RNA transport`=0) %>% 
  filter(Sample %in% c('DUMMY','5FU','10mM_Micit','5FU_10mM_Micit')) %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('DUMMY',
                                    '5FU', 
                                    '10mM_Micit', 
                                    '5FU_10mM_Micit'
                                    ))) %>% 
  ggradar(
    values.radar = c("-1.5", "0", "1"),
    grid.min = -1.5, grid.mid = 0, grid.max = 1,
    group.line.width = 1, 
    group.colours = c('grey70',
                      '#FA2713', # 5FU
                      '#F79E05', # 10 micit
                      '#AD4413'),
    group.point.size = 3,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "bottom"
  ) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_categories_v1.pdf'),
             width = 12, height = 11, useDingbats = FALSE)





# revisiting radar plots --------------------------------------------------


### barplots of selected pathways ####
# zscores_paths = 
# zscores_paths = data_long %>% 
#   separate_rows(KEGG_name, sep=';') %>% 
#   filter(!(Gene_names %in% removals)) %>% 
#   filter(Gene_names %in% prots_sig) %>% 
#   arrange(Gene_names) %>% 
#   filter(KEGG_name %in%  kegg_select) %>% 
#   unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
#   group_by(ID, Sample) %>% 
#   summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
#   mutate(z_score = scale(Mean)) %>%
#   select(-Mean) %>% 
#   separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
#   arrange(Pathway) %>% 
#   group_by(Pathway, Sample) %>% 
#   summarise(zMean = mean(z_score,na.rm = TRUE),
#             zSD = sd(z_score, na.rm = TRUE),
#             zSEM = zSD/sqrt(n())) %>% 
#   ungroup



selected_stats = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  group_by(Gene_names, KEGG_name) %>% 
  t_test(Intensity ~ Sample, 
         ref.group = 'Control', 
         detailed = TRUE, 
         p.adjust.method = 'fdr')

selected_stats %>% 
  filter(KEGG_name == "p53 signaling pathway") %>% 
  filter(group2 %in% c("10mM_Micit")) %>% 
  write_csv("summary/p53_proteins_micit.csv")


# selected_stats_sig = data_long %>% 
#   filter(Gene_names %in% prots_sig) %>% 
#   separate_rows(KEGG_name, sep=';') %>% 
#   filter(!(Gene_names %in% removals)) %>% 
#   filter(KEGG_name %in%  kegg_select) %>% 
#   group_by(Gene_names, KEGG_name) %>% 
#   t_test(Intensity ~ Sample, ref.group = 'Control', detailed = TRUE)



estimate_paths = selected_stats %>% 
  filter(p.adj < 0.05) %>%
  # group_by(group2) %>% 
  # mutate(z_score = scale(estimate), .before = estimate) %>% 
  mutate(estimate = -estimate) %>% 
  group_by(group2, KEGG_name) %>% 
  summarise(Mean_score = mean(estimate),
            SD_score = sd(estimate),
            SEM = SD_score/sqrt(n())) %>% 
  rename(Sample = group2)
  


estimate_paths %>% 
  write_csv("summary/FIGURES_paper/log2FC_pathways.csv")

estimate_paths %>% 
  ungroup %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('Control', 
                                    '5FU',
                                    '1mM_Micit',
                                    '5FU_1mM_Micit',
                                    '10mM_Micit',
                                    '5FU_10mM_Micit'))) %>% 
  ggplot(aes(KEGG_name, Mean_score, fill = Sample)) +
  geom_histogram(stat='identity', position = 'dodge2') +
  geom_errorbar(aes(ymin = Mean_score - SEM, ymax = Mean_score + SEM), 
                position = position_dodge(0.9), width = 0) +
  scale_fill_manual(values = c(
                               '#FA2713', # 5FU
                               '#F0C600', # 1 micit
                               '#E05C10', # 5fu + 1 micit
                               '#F79E05', # 10 micit
                               '#AD4413'  # 5fu + 10 micit
  )) +
  labs(
    x = NULL,
    y = 'Fold Change (+- SEM)',
    fill = 'Condition'
  ) +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 7),
    strip.background = element_rect(fill = 'white', color = 'black')
  ) +
  facet_wrap(~KEGG_name, scales = 'free')


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'FC_paths_estimates.pdf'),
             width = 10, height = 9, useDingbats = FALSE)




# spider plot with FC -----------------------------------------------------


## create the spider plot #####


mins = c(min(estimate_paths$Mean_score))
maxs = c(max(estimate_paths$Mean_score))

# define values that are symmetric
mins = -0.6
maxs = 0.6

estimate_wide = estimate_paths %>% 
  select(Pathway = KEGG_name, Sample, Mean_score) %>% 
  pivot_wider(names_from = Sample, values_from = Mean_score) %>% 
  mutate(Max = maxs, Min = mins, Limit = 0, .before = `10mM_Micit`)

# zscores_wide = zscores_paths %>% 
#   select(Pathway:zMean) %>% 
#   pivot_wider(names_from = Sample, values_from = zMean) %>%
#   mutate(Max = maxs, Min = mins, Limit = 0, .before = `10mM_Micit`)

zscores_mat = as.matrix(estimate_wide[,2:9])

rownames(zscores_mat) = estimate_wide$Pathway

zscores_mat = t(zscores_mat)



zscores_mat = as.data.frame(zscores_mat[1:8,])

# HERE ------------

zscores_mat[c(1:3,4,6,7),] 

select_zscores = zscores_mat[c(1:3,4,6,7),] 

names(select_zscores) = c("Cell cycle", "Citrate cycle\n(TCA cycle)", 
                          "ErbB signaling\npathway", "Fatty acid\nmetabolism", 
                          "Focal adhesion", "Gap junction", 
                          "MAPK signaling\npathway" ,
                          "mTOR signaling\npathway", "p53 signaling\npathway", 
                          "Purine\nmetabolism", "Pyrimidine\nmetabolism", 
                          "Ribosome", "RNA\ntransport" )

opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(1,4))

radarchart(
  select_zscores,
  caxislabels = c(-0.5, -0.25, 0, 0.25, 0.5),
  vlcex = 1.1,
  pfcol = c("#99999980",NA,NA,NA),
  pcol= c(NA,
          '#F2AB0C',
          '#E60B1A',
          '#1700F2'), plty = 1, plwd = 2
)


row_nams = c('10mM Micit',
             '5-FU',
             '5-FU + Micit')
legend(
  # x = "bottom",
  x = 0.85,y = 1.1,
  legend = row_nams, horiz = F,
  bty = "n", pch = 20 ,col= c('#F2AB0C',
                              '#E60B1A',
                              '#1700F2'),
  text.col = "black", cex = 1, pt.cex = 2.5
)

par <- par(opar) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_5FU_Micit_FC_comparison_sigprots.pdf'),
             width = 7, height = 6.5, useDingbats = FALSE)


# protein means comparison ------------------------------------------------

# Here I'll do the same comparison as for the radar plots, but 
# using the ratio of Treatment over control, might be more visual
# It's important to notice that this will NOT be a different analysis
# only a different way to represent things


ratios_paths = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  group_by(ID) %>% 
  mutate(ratio = Mean - Mean[Sample == 'Control']) %>% 
  filter(Sample != 'Control') %>% 
  select(-Mean) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway) %>% 
  group_by(Pathway, Sample) %>% 
  summarise(Mean = mean(ratio ,na.rm = TRUE),
            SD = sd(ratio, na.rm = TRUE),
            SEM = SD/sqrt(n())) %>% 
  ungroup



ratios_paths %>% 
  ungroup %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('5FU',
                                    '1mM_Micit',
                                    '10mM_Micit',
                                    '5FU_1mM_Micit', '5FU_10mM_Micit'))) %>% 
  ggplot(aes(Pathway, Mean, fill = Sample)) +
  geom_histogram(stat='identity', position = 'dodge2') +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), 
                position = position_dodge(0.9), width = 0) +
  scale_fill_manual(values = c(
                               '#FA2713', # 5FU
                               '#F0C600', # 1 micit
                               '#F79E05', # 10 micit
                               '#E05C10', # 5fu + 1 micit
                               '#AD4413'  # 5fu + 10 micit
  )) +
  labs(
    x = 'KEGG Pathway',
    y = 'Ratio agains control average (+- SEM)',
    fill = 'Condition'
  ) +
  facet_wrap(~Pathway, scales = 'free')


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'ratios_means_sigProts.pdf'),
             width = 10, height = 9, useDingbats = FALSE)






## create the spider plot #####


mins = c(min(ratios_paths$Mean))
maxs = c(max(ratios_paths$Mean))

# define values that are symmetric
mins = -0.5
maxs = 0.25

ratios_wide = ratios_paths %>% 
  select(Pathway:Mean) %>% 
  pivot_wider(names_from = Sample, values_from = Mean) %>%
  mutate(Max = maxs, Min = mins, Limit = 0, .before = `10mM_Micit`)

ratios_mat = as.matrix(ratios_wide[,2:9])

rownames(ratios_mat) = ratios_wide$Pathway

ratios_mat = t(ratios_mat)



ratios_mat = as.data.frame(ratios_mat[1:8,])


#~~~~~~~~~~~~~~~~~

## This piece of code will plot all Cell samples
# per separate in a single plot, might be good 
# for comparison

opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(0.8,4))
par(mfrow = c(3,2))
# Produce a radar-chart for each student
for (i in 4:nrow(ratios_mat)) {
  radarchart(
    ratios_mat[c(1:3, i), ],
    pfcol = c("#99999980",NA),
    pcol= c(NA,2), plty = 1, plwd = 2,
    title = row.names(ratios_mat)[i]
  )
}
# Restore the standard par() settings
par <- par(opar) 


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_ratios_ALL_samples.pdf'),
             width = 15, height = 10, useDingbats = FALSE)


#~~~~~~~~~~~~

ratios_mat_select = ratios_mat[c(1:3,4,6,7),] 


opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(1,4))

radarchart(
  ratios_mat_select,
  caxislabels = c(-1.5, -0.75, 0, 0.75, 1.5),
  pfcol = c("#99999980",NA,NA,NA),
  pcol= c(NA,
          '#F2AB0C',
          '#E60B1A',
          '#1700F2'), plty = 1, plwd = 2
)

legend(
  # x = "bottom",
  x = -1,y = -1.2,
  legend = rownames(ratios_mat_select[-c(1,2,3),]), horiz = TRUE,
  bty = "n", pch = 20 ,col= c('#F2AB0C',
                              '#E60B1A',
                              '#1700F2'),
  text.col = "black", cex = 1, pt.cex = 2.5
)

par <- par(opar) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_ratios_5FU_Micit_comparison.pdf'),
             width = 13, height = 11.5, useDingbats = FALSE)



# Tanara's presentation ---------------------------------------------------

### TYMS plot ####

# plot examples of genes
gene = c('TYMS')
data_long %>% filter(Gene_names %in%  gene) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  # facet_wrap(~Gene_names) +
  labs(y = 'Protein expression (log2)',
       x = 'Condition') +
  scale_x_discrete(labels = c("5FU" = "5-FU",
                              "1mM_Micit" = "1mM",
                              "10mM_Micit" = "10mM",
                              "5FU_1mM_Micit" = "1mM",
                              "5FU_10mM_Micit" = "10mM")) +
  geom_vline(xintercept = 1.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 2, y = 29.2, label = 'atop(bold("5-FU"))', size = 7,
           color = 'orange', parse = T) +
  geom_vline(xintercept = 2.5, linetype="dashed", color = 'grey60') +
  geom_vline(xintercept = 4.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 3.5, y = 29.2, label = 'atop(bold("Micit"))', size = 7, 
           color = 'blue', parse = T) +
  geom_vline(xintercept = 4.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 5.5, y = 29.2, label = 'atop(bold("Micit"))', size = 7, 
           color = 'blue', parse = T) +
  annotate("text", x = 5.5, y = 29.02, label = 'atop(bold(" + 5-FU"))', size = 7, 
           color = 'orange', parse = T) +
  scale_fill_manual(
    values = c('grey60',
               '#F5EB44',
               '#87CFFF',
               '#2733FF',
               '#8DEB99',
               '#00C91B'), 
    labels = c('Control',
               '5-FU',
               '1mM',
               '10mM',
               '1mM + 5-FU',
               '10mM + 5-FU')) +
  theme(axis.text.x = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=15, color = 'black'),
        axis.text.y = element_text(size=15, color = 'black'),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

ggsave(here('presentation','TYMS_boxplot.pdf'), height = 5, width = 6.3)


### TP53 plot ####

# plot examples of genes
gene = c('TP53')
data_long %>% filter(Gene_names %in%  gene) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  # facet_wrap(~Gene_names) +
  labs(y = 'Protein expression (log2)',
       x = 'Condition') +
  scale_x_discrete(labels = c("5FU" = "5-FU",
                              "1mM_Micit" = "1mM",
                              "10mM_Micit" = "10mM",
                              "5FU_1mM_Micit" = "1mM",
                              "5FU_10mM_Micit" = "10mM")) +
  geom_vline(xintercept = 1.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 2, y = 22.7, label = 'atop(bold("5-FU"))', size = 6,
           color = 'orange', parse = T) +
  geom_vline(xintercept = 2.5, linetype="dashed", color = 'grey60') +
  geom_vline(xintercept = 4.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 3.5, y = 22.7, label = 'atop(bold("Micit"))', size = 7,
           color = 'blue', parse = T) +
  geom_vline(xintercept = 4.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 5.5, y = 22.7, label = 'atop(bold("Micit"))', size = 7,
           color = 'blue', parse = T) +
  annotate("text", x = 5.5, y = 22.6, label = 'atop(bold(" + 5-FU"))', size = 7,
           color = 'orange', parse = T) +
  scale_fill_manual(
    values = c('grey60',
               '#F5EB44',
               '#87CFFF',
               '#2733FF',
               '#8DEB99',
               '#00C91B'), 
    labels = c('Control',
               '5-FU',
               '1mM',
               '10mM',
               '1mM + 5-FU',
               '10mM + 5-FU')) +
  # theme(axis.text.x = element_text()) +
  theme(axis.text.x = element_text(size=13, color = 'black', hjust = 0.5),
        axis.text.y = element_text(size=13, color = 'black'),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

ggsave(here('presentation','TP53_boxplot.pdf'), height = 5, width = 6.3)



### TP21 plot ####

# plot examples of genes
gene = c('CDKN1A')
data_long %>% filter(Gene_names %in%  gene) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  # facet_wrap(~Gene_names) +
  labs(y = 'Protein expression (log2)',
       x = 'Condition') +
  scale_x_discrete(labels = c("5FU" = "5-FU",
                              "1mM_Micit" = "1mM",
                              "10mM_Micit" = "10mM",
                              "5FU_1mM_Micit" = "1mM",
                              "5FU_10mM_Micit" = "10mM")) +
  geom_vline(xintercept = 1.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 2, y = 26.2, label = 'atop(bold("5-FU"))', size = 7,
           color = 'orange', parse = T) +
  geom_vline(xintercept = 2.5, linetype="dashed", color = 'grey60') +
  geom_vline(xintercept = 4.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 3.5, y = 26.2, label = 'atop(bold("Micit"))', size = 7,
           color = 'blue', parse = T) +
  geom_vline(xintercept = 4.5, linetype="dashed", color = 'grey60') +
  annotate("text", x = 5.5, y = 26.2, label = 'atop(bold("Micit"))', size = 7,
           color = 'blue', parse = T) +
  annotate("text", x = 5.5, y = 25.94, label = 'atop(bold(" + 5-FU"))', size = 7,
           color = 'orange', parse = T) +
  scale_fill_manual(
    values = c('grey60',
               '#F5EB44',
               '#87CFFF',
               '#2733FF',
               '#8DEB99',
               '#00C91B'), 
    labels = c('Control',
               '5-FU',
               '1mM',
               '10mM',
               '1mM + 5-FU',
               '10mM + 5-FU')) +
  theme(axis.text.x = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=15, color = 'black'),
        axis.text.y = element_text(size=15, color = 'black'),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

ggsave(here('presentation','TP21_boxplot.pdf'),height = 5, width = 6.3)



### heatmaps ####
# this code is a repetition of the above code

heat_pres = heat_mat
# colnames(heat_pres) = c('10mM Metabolite', "1mM Metabolite", '5-FU',
#                         '5FU + 10mM Metabolite', '5FU + 1mM Metabolite', 
#                         'Control')

# simpler names
colnames(heat_pres) = c('10mM', "1mM", '5-FU',
                        '5FU +\n10mM', '5FU +\n1mM', 
                        'Control')


heat_pres = heat_pres[,c(6,3,2,1,5,4)]

ha = HeatmapAnnotation(
  Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
              '5FU + Metabolite', '5FU + Metabolite'),
  col = list(Samples = c("Control" = "grey60", 
                         '5-FU' = 'red',
                         "Metabolite" = "orange", 
                         "5FU + Metabolite" = "blue")
                       ),
                       border = TRUE)



# cell cycle
Heatmap(heat_pres[1:76,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  column_names_rot =0, 
                  cluster_columns = F,
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7),
        column_names_centered = TRUE,
             top_annotation = HeatmapAnnotation( 
               Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                           '5FU + Metabolite', '5FU + Metabolite'),
               col = list(Samples = c("Control" = "grey60", 
                                      '5-FU' = 'red',
                                      "Metabolite" = "orange", 
                                      "5FU + Metabolite" = "blue")
               ),border = TRUE,
               show_legend = FALSE
               )
        )



dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_cell_cycle.pdf'),
             width = 5, height = 5, useDingbats = FALSE)

# TCA cycle
Heatmap(heat_pres[77:103,], 
                 name = "Z-score",
        column_names_centered = TRUE,
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                cluster_columns = F,
                column_names_rot =0, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7),
        top_annotation = HeatmapAnnotation( 
          Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                      '5FU + Metabolite', '5FU + Metabolite'),
          col = list(Samples = c("Control" = "grey60", 
                                 '5-FU' = 'red',
                                 "Metabolite" = "orange", 
                                 "5FU + Metabolite" = "blue")
          ),border = TRUE,
          show_legend = FALSE
        )
        )

dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_TCA.pdf'),
             width = 5, height = 5, useDingbats = FALSE)


# mTOR pathway
Heatmap(heat_pres[362:379,], 
                  name = "Z-score",
                  show_row_names = FALSE,
        cluster_columns = F,
                  row_names_gp = gpar(fontsize = 4),
                  column_names_rot =0, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7),
        column_names_centered = TRUE,
        top_annotation = HeatmapAnnotation( 
          Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                      '5FU + Metabolite', '5FU + Metabolite'),
          col = list(Samples = c("Control" = "grey60", 
                                 '5-FU' = 'red',
                                 "Metabolite" = "orange", 
                                 "5FU + Metabolite" = "blue")
          ),border = TRUE,
          show_legend = FALSE
        )
        )

dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_mTOR_poster.pdf'),
             width = 5, height = 2, useDingbats = FALSE)



# p53
pushViewport(viewport(gp = gpar(fontfamily = "Arial")))
ht = Heatmap(heat_pres[380:405,], 
                 name = "Z-score",
                 show_row_names = T,
                 row_names_gp = gpar(fontsize = 4),
                column_names_rot =0, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        column_names_centered = TRUE,
        top_annotation = HeatmapAnnotation( 
          Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                      '5FU + Metabolite', '5FU + Metabolite'),
          col = list(Samples = c("Control" = "grey60", 
                                 '5-FU' = 'red',
                                 "Metabolite" = "orange", 
                                 "5FU + Metabolite" = "blue")
          ),border = TRUE,
          show_legend = FALSE
        )
        )
draw(ht, newpage = FALSE)
popViewport()
dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_p53_poster.pdf'),
             width = 5, height = 2.4, useDingbats = FALSE)



# purine
Heatmap(heat_pres[406:487,], 
                    name = "Z-score",
                    show_row_names = FALSE,
                    row_names_gp = gpar(fontsize = 4),
                   column_names_rot =0, 
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        column_names_centered = TRUE,
        top_annotation = HeatmapAnnotation( 
          Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                      '5FU + Metabolite', '5FU + Metabolite'),
          col = list(Samples = c("Control" = "grey60", 
                                 '5-FU' = 'red',
                                 "Metabolite" = "orange", 
                                 "5FU + Metabolite" = "blue")
          ),border = TRUE,
          show_legend = FALSE
        )
        )


dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_purine.pdf'),
             width = 5, height = 5, useDingbats = FALSE)





# pyrimidine
Heatmap(heat_pres[488:555,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 column_names_rot =0, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        column_names_centered = TRUE,
        top_annotation = HeatmapAnnotation( 
          Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                      '5FU + Metabolite', '5FU + Metabolite'),
          col = list(Samples = c("Control" = "grey60", 
                                 '5-FU' = 'red',
                                 "Metabolite" = "orange", 
                                 "5FU + Metabolite" = "blue")
          ),border = TRUE,
          show_legend = FALSE
        )
        )

dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_pyrimidines.pdf'),
             width = 5, height = 5, useDingbats = FALSE)



# ribosome
Heatmap(heat_pres[556:635,], 
                       name = "Z-score",
                       show_row_names = FALSE,
                       row_names_gp = gpar(fontsize = 4),
                       column_names_rot =0, 
                       column_names_side = "top",
                       column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        column_names_centered = TRUE,
        top_annotation = HeatmapAnnotation( 
          Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                      '5FU + Metabolite', '5FU + Metabolite'),
          col = list(Samples = c("Control" = "grey60", 
                                 '5-FU' = 'red',
                                 "Metabolite" = "orange", 
                                 "5FU + Metabolite" = "blue")
          ),border = TRUE,
          show_legend = FALSE
        )
        )

dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_ribosome.pdf'),
             width = 5, height = 5, useDingbats = FALSE)




# RNA transport
Heatmap(heat_mat[636:757,], 
                          name = "Z-score",
                          show_row_names = FALSE,
                          row_names_gp = gpar(fontsize = 4),
                          column_names_rot =0, 
                          column_names_side = "top",
                          column_names_gp = gpar(fontsize = 7),
                          cluster_columns = F,
                          column_names_centered = TRUE,
                          top_annotation = HeatmapAnnotation( 
                            Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
                                        '5FU + Metabolite', '5FU + Metabolite'),
                            col = list(Samples = c("Control" = "grey60", 
                                                   '5-FU' = 'red',
                                                   "Metabolite" = "orange", 
                                                   "5FU + Metabolite" = "blue")
                            ),border = TRUE,
                            show_legend = FALSE
                          ))


dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_rna_transport.pdf'),
             width = 5, height = 5, useDingbats = FALSE)




# THE heatmap -------------------------------------------------------------



# Let's extract the most significant genes from the main enriched/significant
# pathways: TCA, ribosome, pyrimidine, purine, p53, mTOR

# folate
gene = c("MOCS2","FPGS","GGH","DHFR","DHFR","QDPR","SHMT1","SPR",
                "TYMS","MTHFD1L","MTHFD2","ATIC","MTHFD1","GART","SHMT2")

# pyrimidine
pyrimidine_gene = c("CDA","POLR3G","POLE4","POLR3H","UCKL1","POLR3E","POLR3B","POLD3",
         "POLR3A","POLE2","POLR2J;POLR2J3;POLR2J2","NME7","TYMP","POLR3F",
         "CANT1","POLR2L","RRM2B","ZNRD1","DCTD","TXNRD2","POLR3C","POLR2G",
         "POLR2K","POLR3D","NT5E","POLA2","NME3","NT5C","NT5C3A","DCK","POLR2I",
         "POLR2D","POLA1","UCK2","NT5C2","NUDT2","POLD2","DHODH","PRIM1","POLE3",
         "POLE","PRIM2","AK3","ITPA","DTYMK","POLR2C","POLR1E","POLR1C","POLR1B",
         "CMPK1","POLD1","POLR2H","POLR2E","UMPS","PNPT1","POLR1A","TYMS","DUT",
         "POLR2B","TXNRD1","POLR2A","TK1","CTPS1","PNP","RRM1","RRM2","CAD",
         "NME1","NME2;NME2P1")

gene = statsR %>% 
  filter(Gene_names %in% pyrimidine_gene, 
         Reference == 'Control',
         FDR < 0.05) %>% 
  filter(Target == '5FU_10mM_Micit') %>% 
  pull(Gene_names)

# purine
pur_gene = c("POLR3G","POLE4","GMPR","POLR3H","DGUOK","POLR3E","POLR3B","NTPCR",
             "POLD3","POLR3A","POLE2","POLR2J;POLR2J3;POLR2J2","GUK1","NME7",
             "POLR3F","CANT1","AMPD2","PAPSS2","POLR2L","RRM2B","ZNRD1","POLR3C",
             "POLR2G","POLR2K","POLR3D","NT5E","POLA2","NME3","NT5C","NUDT9",
             "NT5C3A","DCK","GMPR2","POLR2I","POLR2D","POLA1","NT5C2","NUDT2",
             "POLD2","PRIM1","POLE3","POLE","AK1","PRIM2","PGM1","ITPA","IMPDH1",
             "PRPS1","POLR2C","POLR1E","GDA","ADK","PGM2","POLR1C","ADSS","POLR1B",
             "POLD1","AK4","POLR2H","POLR2E","PNPT1","PPAT","APRT","POLR1A","ADSL",
             "PAPSS1","HPRT1","PFAS","NUDT5","POLR2B","POLR2A","PNP","PRPS2",
             "RRM1","RRM2","AK2","GMPS","ATIC","GART","PAICS","NME1","IMPDH2",
             "NME2;NME2P1","PKM")

gene = statsR %>% 
  filter(Gene_names %in% pur_gene, 
         Reference == 'Control',
         FDR < 0.05) %>% 
  pull(Gene_names)



# ribosome genes

ribo_gene = c("RPL26L1","RSL24D1","UBA52","RPS27","RPS29","RPL38","RPL39P5;RPL39",
              "RPL22L1","RPL29","RPS27L","RPL36AL","RPS21","RPS15","RPL22",
              "RPS15A","MRPL13","RPL14","RPL35","RPS25","RPS262;RPS26P11",
              "RPL19","RPS20","RPL34","RPL36","RPLP1","RPS24","RPL13A","RPL37A",
              "RPL27A","RPS17","RPL11","RPL31","RPS23","RPL9","RPS5","RPL24",
              "RPL35A","RPS28","RPL26","RPS10","RPS13","RPL23A","RPL21","RPL32",
              "RPL30","RPS14","RPL18A","RPS12","RPL8","RPL27","RPL23","RPL18",
              "RPL28","RPS19","RPS6","RPS18","RPS16","RPL12","RPL10","RPL17",
              "RPL15","RPL10A","RPS11","RPL7","RPS7","RPL13","RPL7A",
              "RPLP0;RPLP0P6","RPS2","RPS9","RPL5","RPS27A;UBB;UBC",
              "RPLP2","RPS8","RPS3A","RPS3","RPS4X","RPL6","RPSA","RPL3","RPL4")


gene = statsR %>% 
  filter(Gene_names %in% ribo_gene, 
         Reference == 'Control',
         FDR < 0.05) %>% 
  pull(Gene_names)



# p53
p53_gene = c("MDM4","TP53","ATR","TP53I3","FAS","CHEK2","CASP3","EI24","CDKN1A",
             "CHEK1","CCND1","CCNB2","RRM2B","BAX","DDB2","CDK4","GTSE1","CDK6",
             "CASP8","BID","CCNB1","SERPINB5","CDK2","RRM2","CDK1","CYCS","SFN")


gene = statsR %>% 
  filter(Gene_names %in% p53_gene, 
         # Reference == 'Control',
         FDR < 0.05) %>% 
  pull(Gene_names)




## TCA
tca_genes = c("SDHC","PC","PCK2","PDHX","IDH2","ACO1","IDH1","IDH3G","SUCLG1",
              "SUCLA2","SUCLG2","SDHB","DLST","IDH3B","OGDH","PDHB","DLAT",
              "PDHA1","ACO2","SDHA","IDH3A","MDH1","FH","DLD","ACLY","CS","MDH2")

gene = statsR %>% 
  filter(Gene_names %in% tca_genes, 
         # Reference == 'Control',
         FDR < 0.05) %>% 
  pull(Gene_names)


# gene = c('PPAT')

data_long %>% filter(Gene_names %in%  gene) %>% 
  filter(Sample %in% c('Control', '5FU', '10mM_Micit', '5FU_10mM_Micit')) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            # '1mM_Micit',
                                            '10mM_Micit',
                                            # '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  # filter(Sample != 'Control') %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  # labs(title = gene) +
  facet_wrap(~Gene_names, scales= 'free_y') +
  theme(axis.text.x = element_text(hjust = 1, angle=45))



# final gene list
gene = c(# folate
         'GART', 'MTHFD1', 'MTHFD1L', 'MTHFD1', 'SPR', 'DHFR','FPGS','TYMS',
         # pyrimidine
         'CAD','NME1','NME2;NME2P1','PNP','NT5C','POLE3','POLR1A','POLR1B',
         'POLR2B','POLR2D','POLR2G','POLR2H','PRIM1','RRM1','RRM2','UCK2','TK1','UMPS','TYMS',
         # purine
        'ATIC','GART','ADSS','GMPS','IMPDH1','PAICS','PFAS','POLA1','PPAT','PRIM1','PRPS2',
         # ribosome
         'RPL10','RPL11','RPL13','RPS14','RPS13','RPS10',
         # p53
         'TP53','CDKN1A','CYCS','DDB2','FAS','RRM2', 'SFN', 'SERPINB5','CCNB2','CDK2',
         # TCA
         'DLD','MDH2','IDH2','PC','SUCLG1','SUCLG2','PCK2'
         )




### Simple heatmap ####

conference_heat = heat_mat_all[unique(gene),]

# simpler names
colnames(conference_heat) = c('10mM', "1mM", '5-FU',
                        '5FU +\n10mM', '5FU +\n1mM', 
                        'Control')

conference_heat = conference_heat[,c(6,3,2,1,5,4)]

conference_heat

Heatmap(conference_heat[,c(1,2,4,6)],
        show_heatmap_legend = F,
        # name = "z-score",
        row_names_gp = gpar(fontsize = 5),
        column_names_rot =0, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 4),
        column_names_centered = TRUE,
        top_annotation = HeatmapAnnotation(
          Samples = c('Control','5-FU', 'Metabolite',
                      '5FU + Metabolite'),
          col = list(Samples = c("Control" = "grey60",
                                 '5-FU' = 'red',
                                 "Metabolite" = "orange",
                                 "5FU + Metabolite" = "blue")
          ),border = TRUE,
          show_legend = FALSE),
        cluster_columns = F
        
        )










### category selection genes ####


folate = c('GART', 'MTHFD1', 'MTHFD1L', 'MTHFD1', 'SPR', 'DHFR','FPGS','TYMS')
pyr = c('CAD','NME1','NME2;NME2P1','PNP','NT5C','POLE3','POLR1A','POLR1B',
        'POLR2B','POLR2D','POLR2G','POLR2H','PRIM1','RRM1','RRM2','UCK2','TK1','UMPS','TYMS')
pur = c('ATIC','GART','ADSS','GMPS','IMPDH1','PAICS','PFAS','POLA1','PPAT','PRIM1','PRPS2')
rib = c('RPL10','RPL11','RPL13','RPS14','RPS13','RPS10')
p53 = c('TP53','CDKN1A','CYCS','DDB2','FAS','RRM2', 'SFN', 'SERPINB5','CCNB2','CDK2')
tca = c('DLD','MDH2','IDH2','PC','SUCLG1','SUCLG2','PCK2')

## helping functions

heat_mat_gen = function(mat=heat_mat_all,genes) {
  temp_mat = mat[genes,]
  colnames(temp_mat) = c('10mM', "1mM", '5-FU',
                         '5FU +\n10mM', '5FU +\n1mM', 
                         'Control')
  temp_mat = temp_mat[,c(6,3,2,1,5,4)]
  return(temp_mat)
}

draw_HT_with_top = function(temp_mat = temp_mat, pathway, show_names = T) {
  Heatmap(temp_mat[,c(1,2,4,6)], 
          # name = "Z-score",
          show_heatmap_legend = F,
          show_row_names = show_names,
          row_names_gp = gpar(fontsize = 6),
          column_names_rot =0, 
          column_names_side = "top",
          cluster_columns = F,
          # annotations
          top_annotation = HeatmapAnnotation(
            Samples = c('Control','5-FU', 'Metabolite',
                        '5FU + Metabolite'),
            col = list(Samples = c("Control" = "grey60",
                                   '5-FU' = 'red',
                                   "Metabolite" = "orange",
                                   "5FU + Metabolite" = "blue")
            ),border = TRUE,
            show_legend = FALSE),
          right_annotation = rowAnnotation(
            Pathway = rep(pathway,dim(temp_mat)[1]),
            show_annotation_name = F),
          column_names_gp = gpar(fontsize = 7))
  
}

draw_HT = function(temp_mat = temp_mat, pathway, show_names = T) {
 Heatmap(temp_mat[,c(1,2,4,6)], 
          # name = "Z-score",
         show_heatmap_legend = F,
          row_names_gp = gpar(fontsize = 6),
         show_row_names = show_names,
          column_names_rot =0, 
          column_names_side = "top",
          cluster_columns = F,
          # annotations
         right_annotation = rowAnnotation(
            Pathway = rep(pathway,dim(temp_mat)[1]),
            show_annotation_name = F),
          column_names_gp = gpar(fontsize = 7))
  
}





# draw individual heatmaps
# folate
temp_mat = heat_mat_gen(genes = folate)
fol_ht = draw_HT_with_top(temp_mat , pathway = 'Folate', show_names = F)

# pyrimidines
temp_mat = heat_mat_gen(genes = pyr)
pyr_ht = draw_HT(temp_mat , pathway = 'Pyrimidines', show_names = F)

# purines
temp_mat = heat_mat_gen(genes = pur)
pur_ht = draw_HT(temp_mat , pathway = 'Purines', show_names = F)

# ribosome
temp_mat = heat_mat_gen(genes = rib)
rib_ht = draw_HT(temp_mat , pathway = 'Ribosome',show_names = F)

# p53
temp_mat = heat_mat_gen(genes = p53)
p53_ht = draw_HT(temp_mat , pathway = 'p53',show_names = F)

# tca
temp_mat = heat_mat_gen(genes = tca)
tca_ht = draw_HT(temp_mat , pathway = 'TCA cycle',show_names = F)


pathways = c('Folate', 'Pyrimidines','Purines','Ribosome','p53','TCA cycle')

ht_list = fol_ht %v% pyr_ht %v% pur_ht %v% 
  rib_ht %v% p53_ht %v% tca_ht 

lgd = Legend(labels = pathways, title = "", 
             legend_gp = gpar(fill = 1:6),
             grid_height = unit(0, "cm"), grid_width = unit(0, "mm"),
             labels_gp = gpar(col = "white", fontsize = 0))

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
zlgd = Legend(col_fun = col_fun, title = "Z-score")

pd = packLegend(lgd, zlgd)

draw(ht_list, annotation_legend_list = pd, annotation_legend_side = "right") 



dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_KEGG.pdf'),
             width = 6, height = 8, useDingbats = FALSE)

dev.copy2pdf(device = cairo_pdf,
             file = here('presentation', 'heatmap_KEGG_poster.pdf'),
             width = 6, height = 3.4, useDingbats = FALSE)




# e2f targets -------------------------------------------------------------





e2f = c('CCND1','CCND3','JUN','MYC','MYCN','CCNE1','CCNE2','CDC25A','CDK2','E2F1','E2F2','E2F3',
        'NPAT','MYB','MYBL2','TFDP1','AURKB','CCNA1','CCNA2','CDC25A','CDC20','CKS1B','CKS2','NDC80',
        'MKI67','KIF4A','KIF22','PLK1','PRC1','SMC2','SMC4','AURKB','CDC14B','CDC45','CDC6',
        'CDC7','CDT1','DCK','DHFR','DUT','LIG1','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','ORC1',
        'PCNA','POLA1','POLA2','POLD1','PRIM2','RFC1','RFC2','RFC3','RFC4','RPA1','RPA2','RPA3','RRM1',
        'RRM2','TK1','TOP2A','TYMS','CDKN1C','CDKN2C','CDKN2D','E2F7','RB1','RBL1','BRCA1','BRCA2',
        'BUB1','BUB1B','BUB3','CENPE','CHEK1','MAD2L1','TP53','TTK','BARD1','CSTF1','FEN1','MGMT',
        'MLH1','MSH2','MSH6','PMS2','PRKDC','RAD51','RAD54L','UNG','UNG2','APAF1','BAD','BAK1',
        'BCL2','BID','BOK','CASP3','CASP7','CASP8','MAP3K14','MAP3K5','TP73','NKX3-2','EED','EN2',
        'EZH2','FOS','HEY1','HOXA4','HOXA5','HOXA7','HOXA9','HOXA10','HOXA11','HOXB9','HOXD8',
        'PITX1','SIX1','SUZ12','BMP2','FST','TGFA','PPARGC1A','JUNB','TEAD4')


e2f_res = statsR %>% 
  filter(contrast == '10mMmicit - Ctrl') %>% 
  filter(Gene_names %in% e2f) %>% 
  mutate(Direction = case_when(estimate < 0 ~ 'Down',
                               estimate > 0 ~ 'Up')) %>% 
  select(estimate, FDR, FDR_stars, Direction, Gene_names) 

e2f_res 


# arrange table by vector
e2f_res[match(e2f, e2f_res$Gene_names),] %>% 
  write_csv(here('summary','e2f_targets_hct116_prots.csv'))








# multiomics --------------------------------------------------------------



names(data_tidy)

data_tidy %>% 
  select(Gene_names, Control_1, Control_2, Control_3, Control_4,
         `10mM_Micit_1`, `10mM_Micit_2`, `10mM_Micit_3`, `10mM_Micit_4`) %>% 
  # replace_na(list(Control_1 = 1, Control_2 = 1, Control_3 = 1, Control_4 = 1,
  #                 `10mM_Micit_1` = 1, `10mM_Micit_2` = 1, `10mM_Micit_3` = 1, `10mM_Micit_4` = 1)) %>% 
  write_csv(here('summary', 'prot_micit_multiomics.csv'))





# PAPER heatmaps ----------------------------------------------------------



### p53 genes ####

p53 =  data %>% 
  filter(str_detect(KEGG_name, 'p53 signaling pathway')) %>% 
  select(Gene_names) %>% pull


# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)


data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)[,1]) %>% 
  mutate(Pathway = 'P53', .before = "Mean") %>% 
  write_csv("summary/FIGURES_paper/p53_zscore.csv")



# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

# heat_mat = heat_mat[,c('Control','5FU','1mM_Micit','10mM_Micit','5FU_1mM_Micit','5FU_10mM_Micit')]
heat_mat = heat_mat[,c('Control','5FU','10mM_Micit','5FU_10mM_Micit')]

# colnames(heat_mat) = c('Control', '5-FU', "1mM", '10mM', 
#                         '5FU +\n1mM', '5FU +\n10mM')

colnames(heat_mat) = c('Control', '5-FU', '10mM', '5FU +\n10mM')

# 
# ha = HeatmapAnnotation(
#   Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
#               '5FU + Metabolite', '5FU + Metabolite'),
#   col = list(Samples = c("Control" = "grey60", 
#                          '5-FU' = 'red',
#                          "Metabolite" = "orange", 
#                          "5FU + Metabolite" = "blue")
#   ),
#   border = TRUE)

ha = HeatmapAnnotation(
  Samples = c('Control','5-FU','Metabolite', '5FU + Metabolite'),
  col = list(Samples = c("Control" = "grey60", 
                         '5-FU' = 'red',
                         "Metabolite" = "orange", 
                         "5FU + Metabolite" = "blue")
  ),
  border = TRUE)


ht = Heatmap(heat_mat, 
        name = "z-score",
        # column_km = 3,
        # row_km = 3,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 8),
        column_names_rot = 0, 
        column_names_side = "top",
        top_annotation = ha,
        column_names_centered = TRUE,
        border = TRUE,
        show_heatmap_legend = T,
        heatmap_legend_param = list(direction = "horizontal", 
                                    title_position = "topcenter"),
        # cell_fun = function(j, i, x, y, width, height, fill) {
        #   grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 5))},
        column_names_gp = gpar(fontsize = 6)
)

draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


dev.copy2pdf(device = cairo_pdf,
             file = here('summary/FIGURES_paper', 'heatmap_p53.pdf'),
             width = 5, height = 5, useDingbats = FALSE)



## purine-pyrimidine metab ------------



p53 =  data %>% 
  filter(str_detect(KEGG_name, 'Purine metabolism') | str_detect(KEGG_name, "Pyrimidine metabolism")) %>% 
  select(Gene_names) %>% pull


# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

# heat_mat = heat_mat[,c('Control','5FU','1mM_Micit','10mM_Micit','5FU_1mM_Micit','5FU_10mM_Micit')]
heat_mat = heat_mat[,c('Control','5FU','10mM_Micit','5FU_10mM_Micit')]

# colnames(heat_mat) = c('Control', '5-FU', "1mM", '10mM', 
#                         '5FU +\n1mM', '5FU +\n10mM')

colnames(heat_mat) = c('Control', '5-FU', '10mM', '5FU +\n10mM')

# 
# ha = HeatmapAnnotation(
#   Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
#               '5FU + Metabolite', '5FU + Metabolite'),
#   col = list(Samples = c("Control" = "grey60", 
#                          '5-FU' = 'red',
#                          "Metabolite" = "orange", 
#                          "5FU + Metabolite" = "blue")
#   ),
#   border = TRUE)

ha = HeatmapAnnotation(
  Samples = c('Control','5-FU','Metabolite', '5FU + Metabolite'),
  col = list(Samples = c("Control" = "grey60", 
                         '5-FU' = 'red',
                         "Metabolite" = "orange", 
                         "5FU + Metabolite" = "blue")
  ),
  border = TRUE)


ht = Heatmap(heat_mat, 
             name = "z-score",
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 8),
             column_names_rot = 0, 
             column_names_side = "top",
             top_annotation = ha,
             column_names_centered = TRUE,
             border = TRUE,
             show_heatmap_legend = T,
             heatmap_legend_param = list(direction = "horizontal", 
                                         title_position = "topcenter"),
             column_names_gp = gpar(fontsize = 6)
)

draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


dev.copy2pdf(device = cairo_pdf,
             file = here('summary/FIGURES_paper', 'heatmap_nucleotide_metab.pdf'),
             width = 6, height = 12, useDingbats = FALSE)





## pyrimidine metab ------------



p53 =  data %>% 
  filter(str_detect(KEGG_name, "Pyrimidine metabolism")) %>% 
  select(Gene_names) %>% pull

# refined list

p53 = c("AK3","CAD","CDA", "CMPK1", "CTPS1", "DCTD","DHODH", "DTYMK", 
        "DUT", "ITPA", "NME1","NME2;NME2P1","NME3","NME7","NT5C",
        "NT5C2","NT5C3A","NT5E","PNP","RRM1","RRM2","RRM2B","TK1",
        "TXNRD1","TXNRD2","TYMP","TYMS", "UCK2", "UCKL1","UMPS")


# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)


data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)[,1]) %>% 
  mutate(Pathway = "Pyrimidine metabolism", .before = "Mean") %>% 
  write_csv("summary/FIGURES_paper/pyr_metab_zscore.csv")

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

# heat_mat = heat_mat[,c('Control','5FU','1mM_Micit','10mM_Micit','5FU_1mM_Micit','5FU_10mM_Micit')]
heat_mat = heat_mat[,c('Control','5FU','10mM_Micit','5FU_10mM_Micit')]

# colnames(heat_mat) = c('Control', '5-FU', "1mM", '10mM', 
#                         '5FU +\n1mM', '5FU +\n10mM')

colnames(heat_mat) = c('Control', '5-FU', '10mM', '5FU +\n10mM')

# 
# ha = HeatmapAnnotation(
#   Samples = c('Control','5-FU', "Metabolite" ,'Metabolite',  
#               '5FU + Metabolite', '5FU + Metabolite'),
#   col = list(Samples = c("Control" = "grey60", 
#                          '5-FU' = 'red',
#                          "Metabolite" = "orange", 
#                          "5FU + Metabolite" = "blue")
#   ),
#   border = TRUE)

ha = HeatmapAnnotation(
  Samples = c('Control','5-FU','Metabolite', '5FU + Metabolite'),
  col = list(Samples = c("Control" = "grey60", 
                         '5-FU' = 'red',
                         "Metabolite" = "orange", 
                         "5FU + Metabolite" = "blue")
  ),
  border = TRUE)


ht = Heatmap(heat_mat, 
             name = "z-score",
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 8),
             column_names_rot = 0, 
             column_names_side = "top",
             top_annotation = ha,
             column_names_centered = TRUE,
             border = TRUE,
             show_heatmap_legend = T,
             heatmap_legend_param = list(direction = "horizontal", 
                                         title_position = "topcenter"),
             column_names_gp = gpar(fontsize = 6)
)

draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


dev.copy2pdf(device = cairo_pdf,
             file = here('summary/FIGURES_paper', 'heatmap_pyrimidine_metab.pdf'),
             width = 5, height = 5, useDingbats = FALSE)
