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

library(broom)
library(plotly)
library(cluster)

theme_set(theme_cowplot(14))




drugs = read_excel("supplementary_tables/Table S7.xlsx", 
           sheet = "drug_fingerprints")


drugs

set.seed(468)

tsne_fit = drugs %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  Rtsne(check_duplicates = FALSE, pca = FALSE,
        num_threads = 8,
        normalize = FALSE, 
        max_iter = 2000, 
        perplexity = 20, theta = 0.5, dims = 3)

# generate data frame from tnse results
tsne.df = data.frame(tsne_fit$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df = tsne.df %>% 
  mutate(Drug = drugs$Drug)

tsne.df %>%
  as_tibble() %>% 
  left_join(drugs %>% 
              select(Drug, Database)) %>% 
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color = ~Database) %>% 
  add_markers()

t1 = list(family = 'Arial', size = 22)
tsne.df %>%
  mutate(DB = drugs$Database) %>% 
  as_tibble() %>% 
  mutate(opacity = case_when(DB == 'DrugBank' ~ 1,
                             TRUE ~ 0.2)) %>% 
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~DB, 
          opacity= ~opacity,
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



### k-clust of tSNE_biolog ####

### get biolog compounds from the previous tSNE ####

# explore the distribution
# keep only the compounds from biolog
tsne.df.biolog = tsne.df %>%
  mutate(DB = drugs$Database) %>% 
  as_tibble() %>% 
  filter(DB == 'Biolog') 


### silhouette plot ####
tsne_mat_biolog = tsne.df.biolog %>% select(Dim1:Dim3) %>% data.frame

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


#### silhouette analysis ####

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


tsne.df.biolog = tsne.df.biolog %>%
  mutate(cluster = res.km$cluster, .before = Drug,
         cluster = as.factor(cluster))

tsne.df.biolog %>%
  plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()






