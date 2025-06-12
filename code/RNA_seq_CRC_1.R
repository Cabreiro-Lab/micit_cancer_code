### RNA seq analysis of human cell lines from Tanara

# In this script we will analyse the RNA seq with the following cell lines:
# 	- HCT116
# 	- LoVo
# 	- SW948 
# 	- DLD-1

# useful links:
# 	- https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
# 	- https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# 	- https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

# maybe use ulimit -s 16384 before start R console

### libraries ####
# library(tximport)
library(tidyverse)
# library(DESeq2)
# notice that DESeq2 library masks 'rename' function from dplyr 
library(ensembldb) # use only if you are going to deal with db
here::set_here()
library(here)
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
# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2



# Sequencing quality ------------------------------------------------------

library(fastqcr)

# Demo QC directory containing zipped FASTQC reports

fastqc_dir = paste(getwd(),'/FastQC',sep='')

qc = qc_aggregate(fastqc_dir)
qc

summary(qc)

qc_stats(qc) %>% 
  mutate(tot.seq = as.numeric(tot.seq)) %>% 
  ggplot(aes(x = tot.seq)) +
  geom_histogram(fill = '#3640D6', color = 'black', bins = 25) +
  # xlim(13694173,24023731) +
  scale_x_continuous(labels = scales::comma_format(big.mark = ".",
                                           decimal.mark = ","),
                     limits = c(14694173,24023731)) +
  labs(x = 'Total seqs',
       y = 'Count')

ggsave(here('summary', 'sequence_length_histogram.pdf'))


qc_stats(qc) %>% 
  mutate(tot.seq = as.numeric(tot.seq)) %>% 
  summarise(Mean = mean(tot.seq),
            SD = sd(tot.seq))


# # # # # # # # # # #
# RNA seq analysis --------------------------------------------------------
# # # # # # # # # # #


# Get sample info ---------------------------------------------------------



samples = read.delim("sampleInfo.txt") 

dir = getwd()
rownames(samples) = samples$Name

quants_dir = 'quants_103'

# load kegg tables from wormenrichr
# kegg = read.delim("KEGG_2019.txt", header = FALSE) 
# kegg = kegg[,-2]
# 
# 
# rownames(kegg) = kegg[,1] ; kegg = kegg[,-1]

# prepare a list with file names
files = file.path(dir,quants_dir, samples$Name, "quant.sf")
names(files) = samples$Name
all(file.exists(files)) # check that files exist




read.samples = function(samp_file = "sampleInfo.txt", quants = 'quants_103') {
  
  samples = read.delim(samp_file) 
  
  dir = getwd()
  rownames(samples) = samples$Name
  
  quants_dir = quants

  # prepare a list with file names
  files = file.path(dir,quants_dir, samples$Name, "quant.sf")
  names(files) = samples$Name

  return(files)
}


read.samples(samp_file = 'sampleInfo_HCT116.txt')

# Get gene labels from databases ------------------------------------------




#### manual id handling  #####

# 
# gtf_file = 'Homo_sapiens.GRCh38.104.gtf.gz'
# gtf_file = 'Homo_sapiens.GRCh38.103.gtf.gz'
# txdb = GenomicFeatures::makeTxDbFromGFF(gtf_file, organism='Homo sapiens')
# 
# 
# k = keys(txdb, keytype="TXNAME")
# tx_map = ensembldb::select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")
# 
# head(tx_map)
# 
# tx2gene = tx_map

# write.csv(tx2gene,file="tx2gene.csv",row.names = FALSE,quote=FALSE)
#
# quants = read_tsv(files[1])
# 
# quants <- separate(quants, Name, c("TXNAME","Number"),remove = FALSE)
# head(quants)
# 
# quants <- left_join(quants, tx_map, by="TXNAME")
# head(quants)
# 
# tx2gene <- dplyr:::select(quants, Name, GENEID)
# head(tx2gene)
# tx2gene <- filter(tx2gene, !is.na(GENEID))




#### using ensembl genomes ####

# let's make our database from ensembldb
ah = AnnotationHub::AnnotationHub(proxy='127.0.0.1:10801')
ahDb = AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb", 103))
ahEdb = ahDb[[1]]
# generate the database
tx2gene.complete = transcripts(ahEdb, return.type = "DataFrame")

# fetch descriptions of genes
info = genes(ahEdb) %>%
  as_tibble() %>%
  dplyr::select(width, gene_id, gene_name, gene_biotype, description, entrezid) %>%
  unnest(cols = c(entrezid))

# join transcription info with gene ids and entrezids
info.join = tx2gene.complete %>%
  as_tibble() %>%
  dplyr::select(tx_id, tx_biotype, gene_id, tx_name) %>%
  left_join(info)

write.csv(info, here('summary','gene_ids_mapping.csv'))

# subset to have tx_id in first column, and gene_id in second
# tx2gene = data.frame(tx2gene.complete[,c(1,7)])


# subset to have tx_id in first column, and gene_id in second
tx2gene = data.frame(tx2gene.complete[,c(9,7)])

colnames(tx2gene) = c('tx_id','gene_id')

head(tx2gene)





### TX import + DESeq ####

# import quantification data 
txi = tximport::tximport(files, type = "salmon", tx2gene = tx2gene)


### starting analysis with DESeq2
# create DESeq data type to be analysed


samples.batch = samples %>% 
  mutate(Batch = as.factor(Batch),
         Condition = as.factor(Condition),
         Sample = as.factor(Sample))

ddsTxi = DESeqDataSetFromTximport(txi, colData = samples.batch, 
                                  design = ~Batch + Sample)



#### filter by rowSums ####
# prefilter, but that might not be necessary
keep = rowSums(counts(ddsTxi)) >= 100

ddsTxi = ddsTxi[keep,]


#### 0s filter ####
# filter by presence of 0s in the samples
# STRATEGY: get sample columns per separate, calculate how many
# columns have 0s, select the genes that have 6 (or more) columns with
# 0s and change the rest of the values to 0


temp = counts(ddsTxi)

cells = unique(samples$Cell_line)

threshold = 6

# loop to cycle for every sample subset and fix weird values
for (cell in cells) {
samp_names = samples %>% 
  filter(Cell_line == cell) %>% 
  pull(Name)
flawed = rownames(temp[,samp_names][rowSums(temp[,samp_names] == 0) >= threshold,])
temp[flawed,samp_names] = 0
}

 
temp %>% view




# Run DESeq2 --------------------------------------------------------------




ddsTxi.filt = DESeqDataSetFromMatrix(temp ,
                       colData = samples.batch,
                       design = ~Batch + Sample)


# 
# temp[flawed,samp_names]





ddsTxi.filt$Sample = relevel(ddsTxi$Sample, ref = "HCT116_C")


# run the Differential Expression Analysis
dds = DESeq(ddsTxi.filt)



# save files for GEO submission -------------------------------------------

samp_names = samples %>% unite(Sample, Sample, Replicate) %>% pull(Sample)
  
samp_names = c("Gene", samp_names)

norm_counts = counts(dds, normalized = TRUE) %>% as_tibble(rownames = "Gene")
raw_counts = counts(dds, normalized = FALSE) %>% as_tibble(rownames = "Gene")


colnames(norm_counts) = samp_names
colnames(raw_counts) = samp_names

norm_counts %>% write_csv("DATA SUBMISSION/norm_counts_deseq2.csv")
raw_counts %>% write_csv("DATA SUBMISSION/raw_counts_deseq2.csv")

### tidy results ####---------
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% cbind(gene_list,.) %>% as_tibble()

gene_counts = gene_counts %>% 
  gather(Name, counts, TVP_1:TVP_32) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(as_tibble(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) 

gene_counts %>% 
  write_csv(here('summary','gene_counts_norm.csv'))

# TYMS, CDKN1A, TP53

gene = 'TYMS'
gene_counts%>% 
  dplyr::filter(gene_name == gene) %>% 
  # dplyr::filter(Cell_line == 'HCT116') %>%
  ggplot(aes(x = Sample, y = counts, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, position = position_jitterdodge())+
  facet_wrap(~gene_id*Cell_line, scales = 'free') +
  scale_fill_manual(values = c("#1C86EE", "#EEC900")) +
  labs(x = 'Sample',
       y = 'Counts (normalised)')

ggsave(here('summary', glue('boxplot_{gene}.pdf')), height = 8, width = 10)

# transofrm data

vsd = vst(dds, blind = FALSE)
rld = rlog(dds, blind = FALSE)


# plot differences between different transformation data methods
df = bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("mean", "sd")  

ggplot(df, aes(x = mean, y = sd)) + 
  geom_hex(bins = 100) +
  # coord_fixed() + 
  facet_grid( . ~ transformation) +
  theme_light()

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'transformation_comparison.pdf'),
             height = 4.5, width = 12, useDingbats = FALSE)


###
# Sample distances

sampleDists = dist(t(assay(rld)))
sampleDists = as.matrix(sampleDists)


# USE THIS TO CHANGE THE ROW/COLUMN NAMES
# names = colnames(sampleDists) %>%
#   str_split('_', simplify = T) %>%
#   data.frame %>% tbl_df() %>%
#   unite(sample, X1, X2, sep = " - ") %>%
#   dplyr::select(sample) %>%
#   t %>% as.vector

names = gene_counts %>% 
  unite(ID, Sample, Replicate) %>% 
  distinct(ID) %>% pull(ID)


colnames(sampleDists) = names; rownames(sampleDists) = names

col_fun = colorRamp2(c(0, 100), c("white", "blue"))
# col_fun(seq(0, 100))
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

Heatmap(sampleDists, name = 'Euclidean \ndistances', 
        col = colors)



dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'Euclidean_distances_samples.pdf'),
             height = 8, width = 9, useDingbats = FALSE)





# PCA ---------------------------------------------------------------------




pcaData = plotPCA(rld, intgroup = c("Cell_line", "Condition"), returnData = TRUE)
pcaData


# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                 scale = c(sd(x) * sc, sd(y) * sc),
                                 centre = c(mean(x), mean(y))))
}



# get info for the ellipses
ell = pcaData %>% group_by(Condition, Cell_line) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
ggplot(pcaData, aes(x = PC1, y = PC2, color = Cell_line, group = interaction(Condition, Cell_line))) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), linetype = Condition), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), 
                               linetype = Condition, fill = Cell_line), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black')) 


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_main_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)


# 
# 
# 
# KEGG databases
# get databases for genes and pathways from KEGG
kegg.links.entrez = limma::getGeneKEGGLinks('hsa', convert = TRUE)
kegg.links.ids = limma::getGeneKEGGLinks('hsa')
path.ids = limma::getKEGGPathwayNames('hsa', remove.qualifier = TRUE)
kegg.links = cbind(kegg.links.entrez, kegg.links.ids[,1])
colnames(kegg.links) = c('entrezid', 'PathwayID', 'KEGG_genes')

kegg.links = kegg.links %>%
  as_tibble %>%
  mutate(entrezid = as.integer(entrezid)) %>%
  left_join(path.ids) %>%
  mutate(PathwayID = str_replace(PathwayID, 'path:hsa', ''))




# General stats -----------------------------------------------------------


# get results and tidy it
res = results(dds) 


# results with different shape of contrasts, tidy
res.hct = results(dds,   contrast = c("Sample", "HCT116_M" , "HCT116_C"))  
res.hct = lfcShrink(dds, contrast = c("Sample", "HCT116_M" , "HCT116_C"), res = res.hct, type = 'ashr')

res.dld = results(dds,  contrast = c("Sample",  "DLD1_M", "DLD1_C")) 
res.dld = lfcShrink(dds, contrast = c("Sample",  "DLD1_M", "DLD1_C"), res = res.dld, type = 'ashr')

res.lovo = results(dds,  contrast = c("Sample",  "LoVo_M", "LoVo_C"))   
res.lovo = lfcShrink(dds, contrast = c("Sample", "LoVo_M", "LoVo_C"), res = res.lovo, type = 'ashr')

res.sw = results(dds, contrast = c("Sample",   "SW948_M", "SW948_C")) 
res.sw = lfcShrink(dds, contrast = c("Sample", "SW948_M", "SW948_C"), res = res.sw, type = 'ashr')

#### tidying results ####
res.hct.tidy = as_tibble(res.hct, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'HCT116',
  Contrast_description = 'Comparison of HCT116 Micit vs HCT116 Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.dld.tidy = as_tibble(res.dld, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'DLD-1',
  Contrast_description = 'Comparison of DLD-1 Micit vs DLD-1 Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.lovo.tidy = as_tibble(res.lovo, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'LoVo',
  Contrast_description = 'Comparison of LoVo Micit vs LoVo Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.sw.tidy = as_tibble(res.sw, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'SW948',
  Contrast_description = 'Comparison of SW948 Micit vs SW948 Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

results.complete = res.hct.tidy %>% rbind(res.dld.tidy, res.lovo.tidy, res.sw.tidy)

results.complete.kegg = results.complete %>% left_join(kegg.links)

# write results in excel files
list_of_datasets = list('HCT116' = res.hct.tidy, 
                        'DLD-1' = res.dld.tidy, 
                        'LoVo' = res.lovo.tidy,
                        'SW948' = res.sw.tidy)

write.xlsx(list_of_datasets, here('summary', 'complete_stats.xlsx'), 
           colNames = T, rowNames = F) 


#### p-val dist ####


results.complete %>% 
  drop_na(pvalue) %>% 
  ggplot(aes(x = pvalue)) +
  geom_histogram(alpha=.8, position = 'identity',
                 fill = '#4A5CE6',
                 bins = 80, color='black') +
  labs(title='Histogram of unadjusted p-values') +
  xlab('Unadjusted p-values') +
  facet_wrap(~Contrast)

ggsave(here('summary', 'unadjusted_pval_dist.pdf'), height = 8, width = 10)

results.complete %>% 
  drop_na(padj) %>% 
  ggplot(aes(x = padj)) +
  geom_histogram(alpha=.8, color='black',
                 fill = '#EE3B3B',
                 position = 'identity', bins = 80) +
  labs(title='Histogram of adjusted p-values (FDR)') +
  xlab('Adjusted p-values (FDR)') +
  facet_wrap(~Contrast)

ggsave(here('summary', 'FDR_pval_dist.pdf'), height = 8, width = 10)



# MA plots ----------------------------------------------------------------



### MA plots for every comparison
plotMA(res.hct,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_HCT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)

plotMA(res.dld,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_DLD.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.lovo,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_LoVo.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.sw,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_SW948.pdf'),
             height = 8, width = 11, useDingbats = FALSE)







# PCA (manual) ------------------------------------------------------------


#### old method ####

###
# PCA data manual way

# Worm = Condition
# Bacteria =  Cell.line

# # pca_data = t(counts(rld, normalized = TRUE))
pca_data = t(assay(rld))

# HCT 116 samples
pca_data = pca_data[c(1,2,9,10,17,18,25,26),]
# DLD_1
pca_data = pca_data[c(3,4,11,12,19,20,27,28),]
# Lovo
pca_data = pca_data[c(3,4,11,12,19,20,27,28)+2,]
# SW948
pca_data = pca_data[c(3,4,11,12,19,20,27,28)+4,]

# # lets compute the PCA
res.pca = PCA(pca_data, scale.unit = FALSE, ncp = 5, graph = F)

# # metadata 
meta_var = samples

# HCT samples
meta_var = samples[c(1,2,9,10,17,18,25,26),]

# DLD_1
meta_var = samples[c(3,4,11,12,19,20,27,28),]

# LoVo
meta_var = samples[c(3,4,11,12,19,20,27,28)+2,]

# SW948
meta_var = samples[c(3,4,11,12,19,20,27,28)+4,]

# # extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], meta_var$Condition,
  meta_var$Cell_line)

colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Condition', 'Cell_line')


# # make a data frame from ellipses
ell = ind_df %>% group_by(Condition, Cell_line) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
line = 'HCT116'
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Condition, group = interaction(Condition, Cell_line))) +
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), linetype = Condition), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), 
                               linetype = Condition, fill = Condition), size = 1, alpha = 0.3) +
  # xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  # ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black'))

ggsave(here('summary',glue('PCA_{line}.pdf')), height = 7, width = 8)







# PCA and UMAP (tidymodels) --------------------------------------------------------------------


#### PCA tidymodels ####

library(tidymodels)
library(ggrepel)



rld_table = t(assay(rld)) %>% as_tibble(rownames = 'Name')

rld_table = samples %>% select(Name, Cell_line, Replicate, Sample) %>% 
  left_join(rld_table) %>% as_tibble() %>% 
  select(-Name) 

# remove SW948 cell line
rld_table = rld_table %>% 
  filter(Cell_line != "LoVo")

pca_rec = recipe(~., data = rld_table) %>%
  update_role(Replicate, Sample, Cell_line, new_role = "id") %>%
  # step_normalize(all_predictors()) %>%
  step_zv(all_predictors()) %>% 
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_pca(all_predictors(), num_comp = 2)

pca_prep = prep(pca_rec)

names(pca_prep)

pca_prep

sdev = pca_prep$steps[[4]]$res$sdev
percent_variation =  100 * round(sdev^2 / sum(sdev^2),4)


tidied_pca = tidy(pca_prep,1)

selection = tidied_pca %>%
  arrange(desc(abs(value))) %>% 
  head(15) %>% 
  pull(terms)

tidied_pca %>% 
  filter(terms %in% selection) %>% 
  filter(component %in% paste0("PC", 1:5)) %>%
  mutate(component = fct_inorder(component)) %>%
  ggplot(aes(value, terms, fill = terms)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~component, nrow = 1) +
  labs(y = NULL) +
  theme_classic()


juice(pca_prep) %>%
  select(-Cell_line) %>% 
  separate(Sample, into = c("Cell_line", "Condition"), sep = "_") %>% 
  ggplot(aes(PC1, PC2, fill = Cell_line)) +
  geom_point(aes(fill = Cell_line,
                 color = Cell_line,
                 shape = Condition), 
             alpha = 0.9, size = 3) +
  # geom_text(check_overlap = F, hjust = "inward") +
  # geom_text_repel(max.overlaps = 14) +
  # labs(color = NULL) +
  # guides(color = NULL) +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_manual(name = "Cell line",
                    values = c("#F0822E",
                               "#5D3EF0",
                               "#3D9B42")) +
  scale_color_manual(name = "Cell line",
                    values = c("#F0822E",
                               "#5D3EF0",
                               "#3D9B42")) +
  theme_cowplot(17, font_family = 'Arial') +
  theme()

ggsave(here('summary', 'PCA_tidymodels_noLoVo.pdf'), height = 6, width = 7)


juice(pca_prep) %>% 
  write_csv("summary/PCA_data_paper.csv")



# calculate by hand -------------------------------------------------------

rld_table = t(assay(rld)) %>% as_tibble(rownames = 'Name')

rld_table = samples %>% select(Name, Cell_line, Replicate, Sample) %>% 
  left_join(rld_table) %>% as_tibble() %>% 
  select(-Name) 

# # remove SW948 cell line
# rld_table = rld_table %>% 
#   filter(Cell_line != "SW948")



rld_table

pca_fit = rld_table %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F) # do PCA on scaled data

eigens = pca_fit %>%
  tidy(matrix = "eigenvalues") %>% 
  mutate(percent = percent * 100)


pca_fit %>%
  augment(rld_table) %>% # add original dataset back in
  select(-Cell_line) %>% 
  separate(Sample, into = c("Cell_line", "Condition"), sep = "_") %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, 
             color = Cell_line,
             shape = Condition)) + 
  geom_point(size = 2) +
  labs(
    x = glue("PC1 ({eigens[1,]$percent}%)"),
    y = glue("PC2 ({eigens[2,]$percent}%)")
  ) +
  theme_half_open(12, font_family = "Arial") 
# +
#   background_grid()

# ggsave(here('summary', 'PCA_noSW.pdf'), height = 6, width = 7)
ggsave(here('summary', 'PCA_REVISION.pdf'), height = 6, width = 7)


#### UMAP ####

library(embed)

umap_rec = recipe(~., data = rld_table) %>%
  update_role(Replicate, Sample, Cell_line, new_role = "id") %>%
  step_umap(all_predictors())



umap_prep = prep(umap_rec)

umap_prep

juice(umap_prep) %>%
  ggplot(aes(umap_1, umap_2, label = Sample)) +
  geom_point(aes(color = Sample), alpha = 0.7, size = 4) +
  geom_text_repel() +
  labs(color = NULL,
       x = 'UMAP 1',
       y = 'UMAP 2') +
  theme_cowplot(17) +
  theme(legend.position = "none")

ggsave(here('summary', 'UMAP_tidymodels_all.pdf'), height = 7, width = 8)

# GSEA --------------------------------------------------------------------



res.hct.tidy %>% 
  filter(str_detect(gene_name,'CPT'))

gene = ''
gene_counts%>% 
  dplyr::filter(gene_name == gene) %>% 
  filter(Cell_line == 'HCT116') %>%
  ggplot(aes(x = Sample, y = counts, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, position = position_jitterdodge(), aes(color = Replicate))+
  facet_wrap(~gene_id, scales = 'free_y')




# DE genes (up/down) to analyse via StringDB

# helper function that extracts DE genes with upper and lower thresholds
DEgenes = function(dataset = results.complete, 
                   cell = 'HCT116', 
                   min_thrs = 0.5,
                   max_thrs = 10) {
  
  up = dataset %>%  
    filter(Contrast == cell) %>% 
    filter(log2FoldChange < max_thrs & log2FoldChange > min_thrs) %>% 
    filter(padj <= 0.05) %>%
    pull(gene_id) %>% unique
  
  
  down = dataset %>% 
    filter(Contrast == cell) %>% 
    filter(log2FoldChange > -max_thrs & log2FoldChange < -min_thrs) %>% 
    filter(padj <= 0.05) %>%
    pull(gene_id) %>% unique
  
  
  up = c('genes', up)
  down = c('genes', down)
  
  up_name = glue('{cell}_UP')
  down_name = glue('{cell}_DOWN')
  
  list_of_datasets = list(
    up_name = up,
    down_name = down
  )
  
  names(list_of_datasets) = c(up_name, down_name)
  return(list_of_datasets)
  # write.xlsx(list_of_datasets, here('summary', glue('{cell}_genes_updown.xlsx')))
}

# HCT
hct_de = DEgenes(results.complete, cell = 'HCT116', min_thrs = 0.3, max_thrs = 10)
# DLD-1
dld_de = DEgenes(results.complete, cell = 'DLD-1', min_thrs = 0, max_thrs = 10)
# LoVo
lovo_de = DEgenes(results.complete, cell = 'LoVo', min_thrs = 0, max_thrs = 10)
# SW948
sw_de = DEgenes(results.complete, cell = 'SW948', min_thrs = 0, max_thrs = 10)

list_of_datasets = list('HCT116_UP' = hct_de[[1]],
                        'HCT116_DOWN'= hct_de[[2]],
                        'DLD_UP' = dld_de[[1]],
                        'DLD_DOWN'= dld_de[[2]],
                        'LoVo_UP' = lovo_de[[1]],
                        'LoVo_DOWN'= lovo_de[[2]]
                        )

write.xlsx(list_of_datasets, here('summary', 'total_genes_updown.xlsx'),overwrite = TRUE)



# HTML report -------------------------------------------------------------

library(ReportingTools)
library(hwriter)


colData(dds)$Sample

colData(dds[,c(1,2,9,10,17,18,25,26)])



## RUN THIS ONLY ONCE
library(biomaRt)
ens.mart = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")

# list the available datasets (species)
listDatasets(ens.mart) %>% 
  filter(str_detect(description, "Human"))

ensembl = useDataset("hsapiens_gene_ensembl", mart=ens.mart)

listAttributes(ensembl) %>% 
  filter(str_detect(name, "go"))

## annotate using biomaRt
## note this is slightly different from what Mike pointed you to, as we
## are calling the IDs 'Ensembl', and are using mgi_symbol instead of hgnc_symbol
add.anns <- function(df, mart, ...) {
  nm <- rownames(df)
  anns <- getBM( attributes = c("ensembl_gene_id", "external_gene_name","description"),
                 filters = "ensembl_gene_id", values = nm, mart = mart)
  anns <- anns[match(nm, anns[, 1]), ]
  colnames(anns) <- c("Ensembl", "Gene Symbol", "Gene Description")
  df <- cbind(anns, df[, 2:ncol(df)])
  rownames(df) <- nm
  df
}

#Add links to Ensembl.org, because that's how we roll.
ensemblLinks <- function(df, ...){
  naind <- is.na(df$Ensembl)
  df$Ensembl <- hwrite(as.character(df$Ensembl), 
                       link = paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
                       as.character(df$Ensembl)), table = FALSE)
  df$Ensembl[naind] <- ""
  return(df)
}



#### HCT report ####

# extract the data for hct cell line
hct.dds = dds[,c(1,2,9,10,17,18,25,26)]


des2Report = HTMLReport(shortName = 'HCT',
                         title = 'RNA-seq analysis of differential expression for HCT116 cells',
                         reportDirectory = "./report")

publish(res.hct,
        des2Report,
        DataSet = hct.dds,
        # pvalueCutoff = 0.05,
        pvalueCutoff = 1,
        make.plots = T,
        # annotation.db = "org.Hs.eg.db",
        .modifyDF = list(add.anns, modifyReportDF,ensemblLinks),
        mart = ensembl,
        factor = colData(hct.dds)$Sample
        )

finish(des2Report)



#### DLD report ####

# extract the data for hct cell line
dld.dds = dds[,c(3,4,11,12,19,20,27,28)]


des2Report = HTMLReport(shortName = 'DLD',
                        title = 'RNA-seq analysis of differential expression for DLD-1 cells',
                        reportDirectory = "./report")

publish(res.dld,
        des2Report,
        DataSet = dld.dds,
        # pvalueCutoff = 0.05,
        pvalueCutoff = 1,
        make.plots = T,
        # annotation.db = "org.Hs.eg.db",
        .modifyDF = list(add.anns, modifyReportDF,ensemblLinks),
        mart = ensembl,
        factor = colData(hct.dds)$Sample
)

finish(des2Report)

#### LoVo report ####

# extract the data for hct cell line
lovo.dds = dds[,c(3,4,11,12,19,20,27,28)+2]


des2Report = HTMLReport(shortName = 'LoVo',
                        title = 'RNA-seq analysis of differential expression for LoVo cells',
                        reportDirectory = "./report")

publish(res.lovo,
        des2Report,
        DataSet = lovo.dds,
        # pvalueCutoff = 0.05,
        pvalueCutoff = 1,
        make.plots = T,
        # annotation.db = "org.Hs.eg.db",
        .modifyDF = list(add.anns, modifyReportDF,ensemblLinks),
        mart = ensembl,
        factor = colData(hct.dds)$Sample
)

finish(des2Report)



# pathway exploration -----------------------------------------------------------

# first, filter out SW948 as it is uninformative for us

results.filt = results.complete %>% 
  filter(Contrast != 'SW948') %>% 
  filter(gene_name != 'TAP2') %>% 
  mutate(Contrast = str_replace(Contrast, 'DLD-1', 'DLD')) %>% 
  select(Contrast, gene_id, gene_name, log2FoldChange, padj, p_adj_stars, Direction) %>% 
  group_by(Contrast) %>% 
  distinct(gene_id, .keep_all = T) %>% 
  ungroup

# filter genes that are, at least in 1 condition, significative

sig.genes = results.filt %>% 
  filter(padj <= 0.05) %>% 
  distinct(gene_id) %>% pull(gene_id)

# generate matrix
results.filt  %>% 
  filter(gene_id %in% sig.genes) %>% 
  select(gene_id, Contrast, log2FoldChange, padj) %>% 
  pivot_wider(names_from = Contrast, values_from = c(log2FoldChange, padj)) 


# see the kegg pathways we have
kegg.links %>% 
  distinct(Description, .keep_all = T) %>% view




pathways = c('Autophagy - other', 'Autophagy - animal',
             'Citrate cycle (TCA cycle)', 'Purine metabolism',
             'Pyrimidine metabolism', 'Folate biosynthesis',
             'mTOR signaling pathway','p53 signaling pathway',
             'Ribosome','RNA transport','FoxO signaling pathway',
             'Cell cycle','Oxidative phosphorylation',
             'Hippo signaling pathway','MicroRNAs in cancer',
             'ErbB signaling pathway','MAPK signaling pathway',
             'Spliceosome','Protein processing in endoplasmic reticulum')

for (path in pathways) {
  
  autoph.genes = kegg.links %>% 
    filter(Description == path) %>% 
    left_join(info) %>% 
    distinct(gene_id) %>% 
    pull(gene_id)
  
  path_res = results.filt  %>% 
    filter(gene_id %in% sig.genes) %>%
    filter(gene_id %in% autoph.genes)
  
  
  path_res %>% 
    ungroup %>% 
    left_join(df) %>% 
    mutate(Contrast = factor(Contrast, levels = c('HCT116', 'DLD', 'LoVo'))) %>%
    ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = 0) +
    geom_tile() +
    geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
    labs(x = 'Cell line',
         y = 'Gene name',
         fill = 'Fold Change\n(log2)') +
    theme(axis.text.y = NULL)
  
  ggsave(here('summary/heatmaps', glue('{path}_heatmap.pdf')), height = 9, width = 12)


}


# HEATMAP PAPER VERSION ---------------------------------------------------



autoph.genes = kegg.links %>% 
  filter(Description %in% c('Citrate cycle (TCA cycle)',
                            'Oxidative phosphorylation')) %>% 
  left_join(info) %>% 
  distinct(gene_id) %>% 
  pull(gene_id)

path_res = results.filt  %>% 
  filter(gene_id %in% sig.genes) %>%
  filter(gene_id %in% autoph.genes)


path_res %>% 
  ungroup %>% 
  left_join(df) %>% 
  mutate(Contrast = factor(Contrast, levels = c('HCT116', 'DLD', 'LoVo'))) %>%
  ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  geom_tile() +
  geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
  labs(x = 'Cell line',
       y = 'Gene name',
       fill = 'Fold Change\n(log2)') +
  theme(axis.text.y = NULL) +
  theme_cowplot(15, font_family = "Arial")

ggsave(here('summary/heatmaps', "TCA_oxphos_heatmap_PAPER.pdf"), height = 12, width = 8)


### plot gene boxplots ####


# merge counts with stats info
results.merge = gene_counts %>% 
  left_join(results.complete %>% 
              mutate(Contrast = case_when(Contrast == 'DLD-1' ~ 'DLD1',
                                          TRUE ~ Contrast)) %>% 
              select(gene_id, Cell_line = Contrast, log2FoldChange, 
                     lfcSE, pvalue, padj, p_adj_stars)) 

## generate gene list from pathways selected above
gene_list = c()
for (path in pathways) {
  
  gene_names = kegg.links %>% 
    filter(Description == path) %>% 
    left_join(info) %>% 
    distinct(gene_name) %>% 
    drop_na(gene_name) %>% 
    pull(gene_name)
  
  gene_list = c(gene_list, gene_names)
  
}

gene_list = unique(gene_list)


# filter genes that are not present
gene_list = results.complete %>% 
  filter(gene_name %in% gene_list) %>% 
  distinct(gene_name) %>% 
  pull(gene_name)

# arrange genes alphabetically
gene_list = gene_list[order(gene_list)]


# generates a boxplot per gene, with P-value and log2FC annotated 
for (gene in gene_list){

  print(glue('Plotting gene {gene}'))
  
  results.merge %>% 
    dplyr::filter(gene_name == gene) %>% 
    # filter(Cell_line == 'HCT116') %>%
    group_by(Cell_line) %>% 
    mutate(position_y = max(counts) + (max(counts) - min(counts))*0.04) %>% 
    ungroup %>% 
    ggplot(aes(x = Sample, y = counts, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 2.5, position = position_jitterdodge())+
    facet_wrap(~gene_id*Cell_line, scales = 'free', ncol = 4) +
    # geom_text(x = 1, y = 5000, size = 10,aes(label = p_adj_stars)) +
        geom_text(x = 1.5, size = 4,aes(y = position_y+(position_y*0.025), 
                                    label = paste('log2FC = ',round(log2FoldChange,3)))) +
    geom_text(x = 1.5, size = 10, color = '#E03636',
              aes(y = position_y, label = p_adj_stars)) +
    scale_fill_manual(values = c("#1C86EE", "#EEC900")) +
    labs(x = 'Sample',
         y = 'Counts (normalised)') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(here('summary/gene_boxplots', glue('boxplot_{gene}.pdf')), height = 8, width = 13)

}



#### plot ALL genes ####

# filter genes that are not present
gene_list = results.complete %>% 
  distinct(gene_name) %>% 
  pull(gene_name)

# arrange genes alphabetically
gene_list = gene_list[order(gene_list)]

# generates a boxplot per gene, with P-value and log2FC annotated 
for (gene in gene_list){
  
  print(glue('Plotting gene {gene}'))
  
  results.merge %>% 
    dplyr::filter(gene_name == gene) %>% 
    # filter(Cell_line == 'HCT116') %>%
    group_by(Cell_line) %>% 
    mutate(position_y = max(counts) + (max(counts) - min(counts))*0.04) %>% 
    ungroup %>% 
    ggplot(aes(x = Sample, y = counts, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 2.5, position = position_jitterdodge())+
    facet_wrap(~gene_id*Cell_line, scales = 'free', ncol = 4) +
    # geom_text(x = 1, y = 5000, size = 10,aes(label = p_adj_stars)) +
    geom_text(x = 1.5, size = 4,aes(y = position_y+(position_y*0.025), 
                                    label = paste('log2FC = ',round(log2FoldChange,3)))) +
    geom_text(x = 1.5, size = 10, color = '#E03636',
              aes(y = position_y, label = p_adj_stars)) +
    scale_fill_manual(values = c("#1C86EE", "#EEC900")) +
    labs(x = 'Sample',
         y = 'Counts (normalised)') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(here('summary/ALL_gene_boxplots', glue('boxplot_{gene}.pdf')), height = 8, width = 13)
  
}


# PCA of sig genes --------------------------------------------------------

library(tidymodels)

results.wide = results.filt  %>% 
  filter(gene_id %in% sig.genes) %>% 
  select(gene_id, Contrast, log2FoldChange) %>% 
  pivot_wider(names_from = Contrast, values_from = c(log2FoldChange)) %>% 
  select(HCT116:LoVo)

t.results =  t(results.wide)

rownames(t.results)

colnames(t.results) = sig.genes

res.pca = PCA(t.results, ncp = 5, scale.unit = T, graph = F)

# # extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2])

colnames(ind_df) = c('Dim1', 'Dim2')

ind_df['Cell'] = rownames(ind_df)

ind_df %>% 
  ggplot(aes(x = Dim1, y = Dim2, color = Cell)) +
  geom_point(size = 9)





# Gene enrichment -------------------------------------------------

# Enrichment Browser

library(EnrichmentBrowser)


## Run this once!

# obtaining gene sets
kegg.gs = getGenesets(org = "hsa", db = "kegg")
go.gs = getGenesets(org = "hsa", db = "go", onto = "BP")





#### HCT ####
# need to do again DESeq2 for hct per separate
files.hct = read.samples(samp_file = 'sampleInfo_HCT116.txt')
# import quantification data 
txi.hct = tximport::tximport(files.hct, type = "salmon", tx2gene = tx2gene)

samples.red = samples.batch %>% filter(Cell_line == 'HCT116')
ddsTxi.hct = DESeqDataSetFromTximport(txi.hct, 
                                      colData = samples.red, 
                                      design = ~Sample)

# filtering by min number of sequences
keep = rowSums(counts(ddsTxi.hct)) >= 40
ddsTxi.hct = ddsTxi.hct[keep,]

ddsTxi.hct$Sample = relevel(ddsTxi.hct$Sample, ref = "HCT116_C")


# run the Differential Expression Analysis
dds.hct = DESeq(ddsTxi.hct)

res.hct.pure = results(dds.hct)

# import for gene enrichment
hct.SE = import(dds.hct, res.hct.pure, from = c('DESeq2'), anno = 'hsa')

# map IDs
hct.SE = idMap(hct.SE, org = "hsa", from = "ENSEMBL", to = "ENTREZID")


head(rownames(hct.SE))

# set based enrichment analysis
sbeaMethods()

# normalize counts
hct.SE = normalize(hct.SE, norm.method = "vst")

# run GSEA enrichment
# kegg 
alpha = 0.3
hct.gsea = sbea(method = "gsea", se = hct.SE, gs = kegg.gs, alpha = alpha)
hct.ora = sbea(method = "ora", se = hct.SE, gs = kegg.gs, alpha = alpha)
hct.padog = sbea(method = "padog", se = hct.SE, gs = kegg.gs, alpha = alpha)

hct.go.gsea = sbea(method = "gsea", se = hct.SE, gs = go.gs, alpha = 0.05)
gsRanking(hct.go.gsea)

# 
# gsRanking(hct.gsea)
# gsRanking(hct.ora)

eaBrowse(hct.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/hct.gsea', 
         report.name = 'hct.gsea')

eaBrowse(hct.ora, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/hct.ora', 
         report.name = 'hct.ora')

eaBrowse(hct.padog, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/hct.padog', 
         report.name = 'hct.padog')


eaBrowse(hct.go.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/GO_BP/hct.go.gsea', 
         report.name = 'hct.go.gsea')

# network regulation analysis

hsa.grn = compileGRN(org="hsa", db="kegg")


nbeaMethods()

nbea.res = nbea(method="ggea", se=hct.gsea, gs=kegg.gs, grn=hsa.grn)

gsRanking(nbea.res)








#### DLD ####

# need to do again DESeq2 for hct per separate
files.dld = read.samples(samp_file = 'sampleInfo_DLD.txt')
# import quantification data 
txi.dld = tximport::tximport(files.dld, type = "salmon", tx2gene = tx2gene)


samples.red = samples.batch %>% filter(Cell_line == 'DLD1')

ddsTxi.dld = DESeqDataSetFromTximport(txi.dld, 
                                      colData = samples.red, 
                                      design = ~Batch + Sample)
# filtering by min number of sequences
keep = rowSums(counts(ddsTxi.dld)) >= 40
ddsTxi.dld = ddsTxi.dld[keep,]


ddsTxi.dld$Sample = relevel(ddsTxi.dld$Sample, ref = "DLD1_C")


# run the Differential Expression Analysis
dds.dld = DESeq(ddsTxi.dld)

res.dld.pure = results(dds.dld)


dld.SE = import(dds.dld, res.dld.pure, from = c('DESeq2'), anno = 'hsa')

dld.SE = idMap(dld.SE, org = "hsa", from = "ENSEMBL", to = "ENTREZID")


head(rownames(dld.SE))

# normalize counts
dld.SE = normalize(dld.SE, norm.method = "vst")

# run GSEA enrichment
# kegg 
alpha = 0.3
dld.gsea = sbea(method = "gsea", se = dld.SE, gs = kegg.gs, alpha = alpha)
dld.ora = sbea(method = "ora", se = dld.SE, gs = kegg.gs, alpha = alpha)
dld.padog = sbea(method = "padog", se = dld.SE, gs = kegg.gs, alpha = alpha)

dld.go.gsea = sbea(method = "gsea", se = dld.SE, gs = go.gs, alpha = 0.05)
 
# gsRanking(hct.gsea)
# gsRanking(hct.ora)

eaBrowse(dld.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/dld.gsea', 
         report.name = 'dld.gsea')

eaBrowse(dld.ora, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/dld.ora', 
         report.name = 'dld.ora')

eaBrowse(dld.padog, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/dld.padog', 
         report.name = 'dld.padog')


eaBrowse(dld.go.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/GO_BP/dld.go.gsea', 
         report.name = 'dld.go.gsea')


#### LoVo ####

# need to do again DESeq2 for hct per separate
files.lovo = read.samples(samp_file = 'sampleInfo_LoVo.txt')
# import quantification data 
txi.lovo = tximport::tximport(files.lovo, type = "salmon", tx2gene = tx2gene)


samples.red = samples.batch %>% filter(Cell_line == 'LoVo')

ddsTxi.lovo = DESeqDataSetFromTximport(txi.lovo, 
                                      colData = samples.red, 
                                      design = ~Batch + Sample)
# filtering by min number of sequences
keep = rowSums(counts(ddsTxi.lovo)) >= 40
ddsTxi.lovo = ddsTxi.lovo[keep,]

ddsTxi.lovo$Sample = relevel(ddsTxi.lovo$Sample, ref = "LoVo_C")


# run the Differential Expression Analysis
dds.lovo = DESeq(ddsTxi.lovo)

res.lovo.pure = results(dds.lovo)




lovo.SE = import(dds.lovo, res.lovo.pure, from = c('DESeq2'), anno = 'hsa')


lovo.SE = idMap(lovo.SE, org = "hsa", from = "ENSEMBL", to = "ENTREZID")


head(rownames(lovo.SE))

# normalize counts
lovo.SE = normalize(lovo.SE, norm.method = "vst")

# run GSEA enrichment
# kegg 
alpha = 0.3
lovo.gsea = sbea(method = "gsea", se = lovo.SE, gs = kegg.gs, alpha = alpha)
lovo.ora = sbea(method = "ora", se = lovo.SE, gs = kegg.gs, alpha = alpha)
lovo.padog = sbea(method = "padog", se = lovo.SE, gs = kegg.gs, alpha = alpha)


lovo.go.gsea = sbea(method = "gsea", se = lovo.SE, gs = go.gs, alpha = 0.05)


# 
# gsRanking(hct.gsea)
# gsRanking(hct.ora)

eaBrowse(lovo.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/lovo.gsea', 
         report.name = 'lovo.gsea')

eaBrowse(lovo.ora, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/lovo.ora', 
         report.name = 'lovo.ora')

eaBrowse(lovo.padog, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/lovo.padog', 
         report.name = 'lovo.padog')

eaBrowse(lovo.go.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/GO_BP/lovo.go.gsea', 
         report.name = 'lovo.go.gsea')



#### SW948 ####

# need to do again DESeq2 for hct per separate
files.sw = read.samples(samp_file = 'sampleInfo_SW.txt')
# import quantification data 
txi.sw = tximport::tximport(files.sw, type = "salmon", tx2gene = tx2gene)


samples.red = samples.batch %>% filter(Cell_line == 'SW948')

ddsTxi.sw = DESeqDataSetFromTximport(txi.sw, 
                                     colData = samples.red, 
                                     design = ~Batch + Sample)
# filtering by min number of sequences
keep = rowSums(counts(ddsTxi.sw)) >= 40
ddsTxi.sw = ddsTxi.sw[keep,]

ddsTxi.sw$Sample = relevel(ddsTxi.sw$Sample, ref = "SW948_C")


# run the Differential Expression Analysis
dds.sw = DESeq(ddsTxi.sw)

res.sw.pure = results(dds.sw)


sw.SE = import(dds.sw, res.sw.pure, from = c('DESeq2'), anno = 'hsa')

sw.SE = idMap(sw.SE, org = "hsa", from = "ENSEMBL", to = "ENTREZID")


head(rownames(sw.SE))

# normalize counts
sw.SE = normalize(sw.SE, norm.method = "vst")

# run GSEA enrichment
# kegg 
alpha = 0.3
sw.gsea = sbea(method = "gsea", se = sw.SE, gs = kegg.gs, alpha = alpha)
sw.ora = sbea(method = "ora", se = sw.SE, gs = kegg.gs, alpha = alpha)
sw.padog = sbea(method = "padog", se = sw.SE, gs = kegg.gs, alpha = alpha)

sw.go.gsea = sbea(method = "gsea", se = sw.SE, gs = go.gs, alpha = 0.05)


# 
# gsRanking(hct.gsea)
# gsRanking(hct.ora)

eaBrowse(sw.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/sw.gsea', 
         report.name = 'sw.gsea')

eaBrowse(sw.ora, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/sw.ora', 
         report.name = 'sw.ora')

eaBrowse(sw.padog, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/sw.padog', 
         report.name = 'sw.padog')

eaBrowse(sw.go.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/GO_BP/sw.go.gsea', 
         report.name = 'sw.go.gsea')






# E2F targets -------------------------------------------------------------


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

e2f_res = res.hct.tidy %>% 
  filter(gene_name %in% e2f) %>% 
  select(gene_name, log2FoldChange, padj, p_adj_stars, Direction) 


# arrange table by vector
e2f_res[match(e2f, e2f_res$gene_name),] %>% 
  write_csv(here('summary','e2f_targets_hct116.csv'))






# Dorothea analysis -----------------------------------------------------------



## We load the required packages
library(dorothea)
library(bcellViper)
library(viper)


# Running VIPER with DoRothEA regulons
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# get the gene counts from HCT116 cell line
hct_counts = gene_counts %>% 
  filter(Cell_line == 'HCT116') %>% 
  select(gene_name, Sample, Replicate, counts) %>% 
  unite(Sample, Sample, Replicate) %>% 
  pivot_wider(names_from = Sample, values_from = counts,
              values_fn = {mean}) %>% # there are a few repeated elements
  data.frame()

# rownames
row.names(hct_counts) = hct_counts[,1]
hct_counts[,1] = NULL

# run viper wrapper with dorothea
tf_hct = run_viper(hct_counts, regulons, tidy = TRUE,
          options =  list(method = "scale", minsize = 4, 
                          eset.filter = FALSE, cores = 1, 
                          verbose = T))

# tidy the results
tidy_tf = tf_hct %>% 
  separate(sample , into = c('Cell', 'Condition', 'Replicate'), sep = '_') %>% 
  as_tibble() %>% 
  group_by(tf, Cell, Condition, confidence) %>% 
  summarise(mean_activity = mean(activity, na.rm = TRUE),
            sd_activity = sd(activity, na.rm = TRUE)) %>% 
  arrange(desc(abs(mean_activity)))  %>% 
  ungroup()


# calculate the stats: Treatment vs Control
tf_hct_stats = tf_hct %>% 
  separate(sample , into = c('Cell', 'Condition', 'Replicate'), sep = '_') %>% 
  as_tibble() %>% 
  group_by(tf, Cell, confidence) %>% 
  rstatix::t_test(activity ~ Condition, 
                  p.adjust.method = "fdr",
                  ref.group = 'C') %>% 
  mutate(p.stars = gtools::stars.pval(p),
         Direction = case_when(statistic < 0 ~ 'Up',
                                      statistic > 0 ~ 'Down',
                                      TRUE ~ 'Neutral'))

# get significant tf
sig_tf = tf_hct_stats  %>% 
  filter(p < 0.05) %>% 
  pull(tf)



#### dot plots #####

# print selection of TFs

tidy_tf %>% 
  left_join(tf_hct_stats) %>% 
  filter(tf %in% sig_tf) %>% 
  arrange(desc(abs(mean_activity))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down'))) %>% 
  filter(tf %in% unique(head(tidy_tf$tf, 25))) %>%
  # drop_na() %>% 
  ggplot(aes(x = fct_reorder(tf, p), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 10.5) +
  theme_cowplot(15) +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity.pdf'),
       height = 8, width = 10)



tidy_tf %>% 
  left_join(tf_hct_stats) %>% 
  filter(tf %in% sig_tf) %>% 
  arrange(desc(abs(mean_activity))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down'))) %>% 
  # filter(tf %in% unique(head(tidy_tf$tf, 25))) %>%
  # drop_na() %>% 
  ggplot(aes(x = fct_reorder(tf, p), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 10.5) +
  theme_cowplot(15) +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity_complete_HCT116.pdf'),
       height = 8, width = 30)



tidy_tf %>% 
  left_join(tf_hct_stats) %>% 
  filter(tf %in% sig_tf) %>% 
  filter(Cell == 'HCT116') %>% 
  arrange(mean_activity) %>% 
  write_csv("summary/tf_activity_data.csv")


# print out the complete set of significant TFs
tidy_tf %>% 
  left_join(tf_hct_stats) %>% 
  filter(tf %in% sig_tf) %>%
  arrange((abs(p))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down')),
         tf = factor(tf)) %>% 
  ggplot(aes(x = fct_reorder(tf, (p)), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 10, angle = 90) +
  theme_cowplot(15) +
  labs(x = 'Transcription factors',
       y = 'Mean activity (+- std)') +
  # facet_wrap(~direction) +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(here('summary', 'tf_activity_complete.pdf'),
       height = 8, width = 30)





#### heatmap ####


tf_mat = tf_hct %>% 
  separate(sample , into = c('Cell', 'Condition', 'Replicate'), 
           sep = '_', remove = FALSE) %>% 
  as_tibble()  %>% 
  left_join(tf_hct_stats) %>% 
  filter(tf %in% sig_tf) %>%
  arrange(p) %>% 
  select(tf, sample, activity) %>% 
  pivot_wider(names_from = sample, values_from = activity) %>% 
  select(-tf) %>% 
  as.matrix()

tf_rownames = tf_hct %>% 
  separate(sample , into = c('Cell', 'Condition', 'Replicate'), 
           sep = '_', remove = FALSE) %>% 
  as_tibble()  %>% 
  left_join(tf_hct_stats) %>% 
  filter(tf %in% sig_tf) %>%
  arrange(p) %>% 
  select(tf, sample, activity) %>% 
  pivot_wider(names_from = sample, values_from = activity) %>% 
  pull(tf)


rownames(tf_mat) = tf_rownames

colnames(tf_mat) = c('Control 1','Micit 1',
                     'Control 2','Micit 2',
                     'Control 3','Micit 3',
                     'Control 4','Micit 4')


subset = c(rownames(head(tf_mat, 20)), 'SMAD3')


Heatmap(tf_mat[subset,],
        name = 'TF activity',
        column_title = "Samples",
        column_title_side = "bottom",
        show_column_dend = FALSE,
        row_dend_side = "right",
        row_dend_width = unit(2, "cm"),
        row_km = 2,
        column_km = 2)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'heatmap_subset.pdf'),
             height = 8, width = 7, useDingbats = FALSE)





Heatmap(tf_mat,
        name = 'TF activity',
        column_title = "Samples",
        column_title_side = "bottom",
        show_column_dend = FALSE,
        row_dend_side = "right",
        row_dend_width = unit(2, "cm"),
        row_km = 2,
        column_km = 2)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'heatmap_complete.pdf'),
             height = 18, width = 7, useDingbats = FALSE)





tf_hct_stats %>% 
  filter(p < 0.05) %>% 
  write_csv(here('summary','tf_stats_sig.csv'))



# save the tables to be analysed with String DB
# down TFs
down_tf = c('genes', tf_hct_stats %>% 
            filter(statistic < 0, p <= 0.05) %>% 
            pull(tf)
)

# up TFs
up_tf = c('genes', tf_hct_stats %>% 
              filter(statistic > 0, p <= 0.05) %>% 
              pull(tf)
)

tf_list = list(
  'TF_UP' = up_tf,
  'TF_DOWN' = down_tf
)


write.xlsx(tf_list, 
           here('summary', 'tf_updown.xlsx'),
           overwrite = TRUE)


#### heatmap with cancer hits ####

# from the enrichment done in String Db, I can read the file and get the genes
# that are related to the cancer categories in KEGG

library(readxl)
TF_output = read_excel("summary/TF_enrich/TF/TF_output.xlsx", 
                        sheet = "KEGG")

# clean a bit this messy file
TF_output = TF_output %>% 
  mutate(inputGenes = str_replace_all(inputGenes, pattern = '\\[', ''),
         inputGenes = str_replace_all(inputGenes, pattern = '\\]', ''),
         inputGenes = str_replace_all(inputGenes, pattern = '\'', ''),
         inputGenes = str_replace_all(inputGenes, pattern = ' ', ''))


cancer_tf = TF_output %>% 
  filter(str_detect(description,'cancer')) %>% 
  filter(fdr < 0.05) %>% 
  filter(direction == 'UP') %>% 
  separate_rows(inputGenes, sep = ',') %>% 
  distinct(inputGenes) %>% pull(inputGenes)
  

Heatmap(tf_mat[cancer_tf,],
        name = 'TF activity',
        column_title = "Samples",
        column_title_side = "bottom",
        show_column_dend = FALSE,
        row_dend_side = "right",
        row_dend_width = unit(2, "cm"),
        # row_km = 2,
        column_km = 2)






# dorothea - others -------------------------------------------------------


#### DLD1 ####

# ad hoc function to get summary stuff

get_tf = function(cell_line = 'HCT116') {
  cat('Formating the data \n')
  # get the gene counts from HCT116 cell line
  hct_counts = gene_counts %>% 
    filter(Cell_line == cell_line) %>% 
    select(gene_name, Sample, Replicate, counts) %>% 
    unite(Sample, Sample, Replicate) %>% 
    pivot_wider(names_from = Sample, values_from = counts,
                values_fn = {mean}) %>% # there are a few repeated elements
    data.frame()
  
  # rownames
  row.names(hct_counts) = hct_counts[,1]
  hct_counts[,1] = NULL
  
  cat('Running viper with the regulons \n')
  # run viper wrapper with dorothea
  tf_hct = run_viper(hct_counts, regulons, tidy = TRUE,
                     options =  list(method = "scale", minsize = 4, 
                                     eset.filter = FALSE, cores = 1, 
                                     verbose = T))
  
  # tidy the results
  tidy_tf = tf_hct %>% 
    separate(sample , into = c('Cell', 'Condition', 'Replicate'), sep = '_') %>% 
    as_tibble() %>% 
    group_by(tf, Cell, Condition, confidence) %>% 
    summarise(mean_activity = mean(activity, na.rm = TRUE),
              sd_activity = sd(activity, na.rm = TRUE)) %>% 
    arrange(desc(abs(mean_activity)))  %>% 
    ungroup()
  
  cat('Computing the stats \n')
  # calculate the stats: Treatment vs Control
  tf_hct_stats = tf_hct %>% 
    separate(sample , into = c('Cell', 'Condition', 'Replicate'), sep = '_') %>% 
    as_tibble() %>% 
    group_by(tf, Cell, confidence) %>% 
    rstatix::t_test(activity ~ Condition, 
                    p.adjust.method = "fdr",
                    ref.group = 'C') %>% 
    mutate(p.stars = gtools::stars.pval(p),
           Direction = case_when(statistic < 0 ~ 'Up',
                                 statistic > 0 ~ 'Down',
                                 TRUE ~ 'Neutral'))
  
  res_list = list(tf_hct, tidy_tf, tf_hct_stats)
  
  return(res_list)
  
}



#### calculate stuff ####

dld_tf_res = get_tf(cell_line = 'DLD1')

dld_tf = dld_tf_res[[1]]
dld_tf_tidy = dld_tf_res[[2]]
dld_tf_stats = dld_tf_res[[3]]



# get significant tf
sig_tf = dld_tf_stats  %>% 
  filter(p < 0.05) %>% 
  pull(tf)



##### plots ####

# print selection of TFs

dld_tf_tidy %>% 
  left_join(dld_tf_stats) %>% 
  filter(tf %in% sig_tf) %>% 
  arrange(desc(abs(mean_activity))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down'))) %>% 
  filter(tf %in% unique(head(tidy_tf$tf, 25))) %>%
  ggplot(aes(x = fct_reorder(tf, p), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 7.5) +
  theme_cowplot(15) +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity_DLD1.pdf'),
       height = 8, width = 10)



# print out the complete set of significant TFs
dld_tf_tidy %>% 
  left_join(dld_tf_stats) %>% 
  filter(tf %in% sig_tf) %>%
  arrange((abs(p))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down')),
         tf = factor(tf)) %>% 
  ggplot(aes(x = fct_reorder(tf, (p)), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 8, angle = 90) +
  theme_cowplot(15) +
  labs(x = 'Transcription factors',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(here('summary', 'tf_activity_complete_DLD1.pdf'),
       height = 8, width = 30)



#### LoVo ####
# "LoVo"   "SW948" 

lovo_tf_res = get_tf(cell_line = 'LoVo')
lovo_tf = lovo_tf_res[[1]]
lovo_tf_tidy = lovo_tf_res[[2]]
lovo_tf_stats = lovo_tf_res[[3]]


# get significant tf
sig_tf = lovo_tf_stats  %>% 
  filter(p < 0.05) %>% 
  pull(tf)


##### plots ####

# print selection of TFs

lovo_tf_tidy %>% 
  left_join(lovo_tf_stats) %>% 
  filter(tf %in% sig_tf) %>% 
  arrange(desc(abs(mean_activity))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down'))) %>% 
  filter(tf %in% unique(head(tidy_tf$tf, 25))) %>%
  ggplot(aes(x = fct_reorder(tf, p), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 10.5) +
  theme_cowplot(15) +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity_LoVo.pdf'),
       height = 8, width = 10)



# print out the complete set of significant TFs
lovo_tf_tidy %>% 
  left_join(lovo_tf_stats) %>% 
  filter(tf %in% sig_tf) %>%
  arrange((abs(p))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down')),
         tf = factor(tf)) %>% 
  ggplot(aes(x = fct_reorder(tf, (p)), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 10.5, angle = 90) +
  theme_cowplot(15) +
  labs(x = 'Transcription factors',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(here('summary', 'tf_activity_complete_LoVo.pdf'),
       height = 8, width = 30)



#### SW948 ####
# "LoVo"   "SW948" 

sw_tf_res = get_tf(cell_line = 'SW948')
sw_tf = sw_tf_res[[1]]
sw_tf_tidy = sw_tf_res[[2]]
sw_tf_stats = sw_tf_res[[3]]


# get significant tf
sig_tf = sw_tf_stats  %>% 
  filter(p < 0.05) %>% 
  pull(tf)


##### plots #####

# print selection of TFs

sw_tf_tidy %>% 
  left_join(sw_tf_stats) %>% 
  filter(tf %in% sig_tf) %>% 
  arrange(desc(abs(mean_activity))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down'))) %>% 
  filter(tf %in% unique(head(tidy_tf$tf, 25))) %>%
  ggplot(aes(x = fct_reorder(tf, p), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 4.4) +
  theme_cowplot(15) +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity_SW948.pdf'),
       height = 8, width = 10)



# print out the complete set of significant TFs
sw_tf_tidy %>% 
  left_join(sw_tf_stats) %>% 
  filter(tf %in% sig_tf) %>%
  arrange((abs(p))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down')),
         tf = factor(tf)) %>% 
  ggplot(aes(x = fct_reorder(tf, (p)), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y =6, angle = 90) +
  theme_cowplot(15) +
  labs(x = 'Transcription factors',
       y = 'Mean activity (+- std)') +
  facet_grid(~Direction, scales = 'free_x', space = 'free') +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(here('summary', 'tf_activity_complete_SW948.pdf'),
       height = 8, width = 20)




#### all together ####


# Running VIPER with DoRothEA regulons
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# get the gene counts from HCT116 cell line
total_counts = gene_counts %>% 
  # filter(Cell_line == 'HCT116') %>% 
  select(gene_name, Sample, Replicate, counts) %>% 
  unite(Sample, Sample, Replicate) %>% 
  drop_na(counts) %>% 
  pivot_wider(names_from = Sample, values_from = counts,
              values_fn = {mean}) %>% # there are a few repeated elements
  data.frame()

# rownames
row.names(total_counts) = total_counts[,1]
total_counts[,1] = NULL

# total_counts['source'] = rowMeans(total_counts)


# run viper wrapper with dorothea
tf_total = run_viper(total_counts, regulons, tidy = TRUE,
                   options = list(method = "scale", minsize = 4, 
                                  eset.filter = TRUE, cores = 1, 
                                  verbose = F))

# tidy the results
tidy_tf = tf_total %>% 
  separate(sample , into = c('Cell', 'Condition', 'Replicate'), sep = '_') %>% 
  as_tibble() %>% 
  group_by(tf, Cell, Condition, confidence) %>% 
  summarise(mean_activity = mean(activity, na.rm = TRUE),
            sd_activity = sd(activity, na.rm = TRUE)) %>% 
  arrange(desc(abs(mean_activity)))  %>% 
  ungroup()


# calculate the stats: Treatment vs Control
tf_stats = tf_total %>% 
  separate(sample , into = c('Cell', 'Condition', 'Replicate'), sep = '_') %>% 
  as_tibble() %>% 
  group_by(tf, Cell, confidence) %>% 
  rstatix::t_test(activity ~ Condition, 
                  p.adjust.method = "fdr",
                  detailed = TRUE) %>% 
  mutate(p.stars = gtools::stars.pval(p),
         Direction = case_when(estimate < 0 ~ 'Up',
                               estimate > 0 ~ 'Down',
                               TRUE ~ 'Neutral'))

# get significant tf
sig_tf = tf_stats  %>% 
  filter(p < 0.05) %>% 
  pull(tf)


tf_stats %>% 
  filter(p < 0.05) %>% 
  write_csv("summary/tf_stats_sig_ALL.csv")




#### dot plots #####

# print selection of TFs

tidy_tf_w_stats = tidy_tf %>% 
  left_join(tf_stats) %>% 
  # filter(tf %in% sig_tf) %>% 
  arrange(desc(abs(mean_activity))) 

tf_selection = tidy_tf_w_stats %>% slice_head(n = 40) %>% distinct(tf) %>% pull(tf)

length(tf_selection)

tidy_tf_w_stats %>% 
  filter(tf %in% tf_selection ) %>%
  ggplot(aes(x = fct_reorder(tf, p), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 5, angle = 90) +
  theme_cowplot(15) +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  # facet_grid(~Direction*Cell, scales = 'free_x', space = 'free') +
  facet_wrap(~Direction*Cell, scales='free_x', ncol = 4) +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity_ALL_top30.pdf'),
       height = 8, width = 23)



# print out the complete set of significant TFs
tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(tf %in% sig_tf) %>% 
  arrange(desc(abs(mean_activity))) %>% 
  mutate(Direction = factor(Direction, levels = c('Up', 'Down'))) %>% 
  ggplot(aes(x = fct_reorder(tf, (p)), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 10, angle = 90) +
  theme_cowplot(15) +
  labs(x = 'Transcription factors',
       y = 'Mean activity (+- std)') +
  facet_wrap(~Direction*Cell, scales='free_x', ncol = 4) +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(here('summary', 'tf_activity_complete_ALL.pdf'),
       height = 8, width = 30)



#### shared TFs in HCT, LoVo and DLD1 -----------

tf_shared_2_down = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell != "SW948") %>% 
  filter(p < 0.05) %>%
  filter(Direction == "Down") %>%
  count(tf) %>% 
  mutate(n = n/2) %>% 
  filter(n >= 2) %>% 
  pull(tf)

tidy_tf_w_stats %>% 
  filter(Cell != "SW948") %>% 
  filter(tf %in% tf_shared_2_down ) %>%
  filter(Direction == "Down") %>% 
  ggplot(aes(x = fct_reorder(tf, estimate), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 5, angle = 90) +
  theme_cowplot(15) +
  background_grid() +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  # facet_grid(~Direction*Cell, scales = 'free_x', space = 'free') +
  facet_wrap(~Direction*Cell, scales='free_x', ncol = 1) +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity_shared_DOWN.pdf'),
       height = 13, width = 15)


##### shared UP -------------------

tf_shared_2_up = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell != "SW948") %>% 
  filter(p < 0.05) %>%
  filter(Direction == "Up") %>%
  count(tf) %>% 
  mutate(n = n/2) %>% 
  filter(n >= 2) %>% 
  pull(tf)

tidy_tf_w_stats %>% 
  filter(Cell != "SW948") %>% 
  filter(tf %in% tf_shared_2_up ) %>%
  filter(Direction == "Up") %>% 
  ggplot(aes(x = fct_reorder(tf, estimate), y = mean_activity)) +
  geom_errorbar(aes(ymax = mean_activity + sd_activity, 
                    ymin = mean_activity - sd_activity),
                width = 0.1, color = 'grey70') +
  geom_point(aes(color = Condition, shape = Condition), 
             size = 5) +
  geom_text(aes(label = p.stars),
            y = 5, angle = 90) +
  theme_cowplot(15) +
  background_grid() +
  labs(x = 'Regulons',
       y = 'Mean activity (+- std)') +
  # facet_grid(~Direction*Cell, scales = 'free_x', space = 'free') +
  facet_wrap(~Direction*Cell, scales='free_x', ncol = 1) +
  scale_color_manual(values = c('#E6454E', '#1F993D')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 

ggsave(here('summary', 'tf_activity_shared_UP.pdf'),
       height = 12, width = 15)




# venn diagram of TF ------------------------------------------------------

#### UP -----------
library(VennDiagram)

tf_hct_up = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell == "HCT116") %>% 
  filter(p < 0.05, Direction == "Up") %>% 
  pull(tf)

tf_lovo_up = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell == "LoVo") %>% 
  filter(p < 0.05, Direction == "Up") %>% 
  pull(tf)

tf_dld_up = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell == "DLD1") %>% 
  filter(p < 0.05, Direction == "Up") %>% 
  pull(tf)

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(tf_hct_up, tf_lovo_up, tf_dld_up),
  category.names = c("HCT116" , "LoVo" , "DLD-1"),
  filename = 'TF_Venn_UP.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 780 , 
  width = 780 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# create a tbl with tfs that are shared

# Combine all TFs and their source
all_tfs <- unique(c(tf_hct_up, tf_lovo_up, tf_dld_up))

tf_table <- tibble(
  TF = all_tfs,
  HCT116 = all_tfs %in% tf_hct_up,
  LoVo = all_tfs %in% tf_lovo_up,
  DLD1 = all_tfs %in% tf_dld_up
)

# Convert logical to presence/absence indicator (optional)
tf_table <- tf_table %>%
  mutate(across(HCT116:DLD1, ~ ifelse(.x, "Yes", "No")))  # Or 1/0


# Arrange for better readability (optional)
tf_table <- tf_table %>%
  arrange(TF)

tf_table %>% 
  write.xlsx("summary/TF_shared_UP.xlsx")




#### DOWN ---------


tf_hct_down = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell == "HCT116") %>% 
  filter(p < 0.05, Direction == "Down") %>% 
  pull(tf)

tf_lovo_down = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell == "LoVo") %>% 
  filter(p < 0.05, Direction == "Down") %>% 
  pull(tf)

tf_dld_down = tidy_tf %>% 
  left_join(tf_stats) %>% 
  filter(Cell == "DLD1") %>% 
  filter(p < 0.05, Direction == "Down") %>% 
  pull(tf)

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(tf_hct_down, tf_lovo_down, tf_dld_down),
  category.names = c("HCT116" , "LoVo" , "DLD-1"),
  filename = 'TF_Venn_DOWN.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 780 , 
  width = 780 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# create a tbl with tfs that are shared

# Combine all TFs and their source
all_tfs <- unique(c(tf_hct_down, tf_lovo_down, tf_dld_down))

tf_table <- tibble(
  TF = all_tfs,
  HCT116 = all_tfs %in% tf_hct_down,
  LoVo = all_tfs %in% tf_lovo_down,
  DLD1 = all_tfs %in% tf_dld_down
)

# Convert logical to presence/absence indicator (optional)
tf_table <- tf_table %>%
  mutate(across(HCT116:DLD1, ~ ifelse(.x, "Yes", "No")))  # Or 1/0


# Arrange for better readability (optional)
tf_table <- tf_table %>%
  arrange(TF)

tf_table %>% 
  write.xlsx("summary/TF_shared_DOWN.xlsx")




# multiomics --------------------------------------------------------------


gene_counts %>% 
  unite(Sample, Sample, Replicate, sep = '_', remove=T) %>% 
  filter(Cell_line == 'HCT116') %>% 
  arrange(Sample) %>% 
  select(gene_name, Sample, counts) %>% 
  # distinct(gene_name, .keep_all = TRUE) %>% 
  pivot_wider(names_from = Sample, values_from = counts, values_fn = {mean}) %>% 
  write_csv(here('summary', 'rnaseq_micit_multiomics.csv'))
  





# Enrich plots presentation -----------------------------------------------

# these plots are for Filipe's (and mine?) presentation
# 


library(readxl)

#### HCT116 ####

HCT116_output = read_excel("summary/GSEA/HCT116/HCT116_output.xlsx", 
                            sheet = "KEGG")



removals = c('Antigen processing and presentation',
             'Epstein-Barr virus infection',
             'Pathogenic Escherichia coli infection',
             'Protein processing in endoplasmic reticulum')

HCT116_output %>% 
  filter(!(description %in% removals)) %>% 
  mutate(direction=factor(direction, levels = c('UP', 'DOWN'))) %>% 
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
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits=rev)

ggsave('summary/GSEA/HCT116/HCT116_heatmap_enrich_KEGG.pdf',
       width = 6, height = 6.5)  


#### DLD ####

DLD_output = read_excel("summary/GSEA/DLD/DLD_output.xlsx", 
                           sheet = "KEGG")



include_dld = c(  # UP
                'Autophagy - animal', 'Apoptosis',
                'Mitophagy - animal',
                'FoxO signaling pathway',
                'Transcriptional misregulation in cancer',
                'ErbB signaling pathway', 'Endocytosis',
                'Colorectal cancer','mTOR signaling pathway',
                'Insulin resistance','p53 signaling pathway',
                'TNF signaling pathway',
                  # DOWN
                'Ribosome', 'Spliceosome', 'Pyrimidine metabolism',
                'RNA transport','Proteasome','Purine metabolism',
                'Cell cycle', 'Glycolysis / Gluconeogenesis',
                'DNA replication','Mismatch repair','RNA polymerase',
                'Citrate cycle (TCA cycle)', 'Pyruvate metabolism',
                'RNA degradation','Carbon metabolism'
                
                )

DLD_output %>% 
  filter((description %in% include_dld)) %>% 
  mutate(direction=factor(direction, levels = c('UP', 'DOWN'))) %>% 
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
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits=rev)

ggsave('summary/GSEA/DLD/DLD_heatmap_enrich_KEGG.pdf',
       width = 6, height = 9.5)  


#### LoVo ####

lovo_output = read_excel("summary/GSEA/LoVo/LoVo_output.xlsx", 
                        sheet = "KEGG")



removals_lovo = c(
  'Chemical carcinogenesis','Steroid hormone biosynthesis',
  'Human papillomavirus infection',
  'Protein processing in endoplasmic reticulum',
  'Antigen processing and presentation','Epstein-Barr virus infection'
)
  

lovo_output %>% 
  filter(!(description %in% removals_lovo)) %>% 
  mutate(direction=factor(direction, levels = c('UP', 'DOWN'))) %>% 
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
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits=rev)

ggsave('summary/GSEA/LoVo/LoVo_heatmap_enrich_KEGG.pdf',
       width = 6, height = 8.5)  






# poster version ----------------------------------------------------------

# join 3 datasets

hct_cats = HCT116_output %>% 
  mutate(direction=factor(direction, levels = c('UP', 'DOWN'))) %>% 
  mutate(FDR = cut(fdr, 
                   labels = c('0.001', '0.01', '0.05', 'ns'),
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf)),
         .before = term) %>% 
  mutate(description = str_wrap(description, width = 25),
         cell = 'HCT 116') %>% 
  filter(description %in% categories)


dld_cats =DLD_output %>% 
  mutate(direction=factor(direction, levels = c('UP', 'DOWN'))) %>% 
  mutate(FDR = cut(fdr, 
                   labels = c('0.001', '0.01', '0.05', 'ns'),
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf)),
         .before = term) %>% 
  mutate(description = str_wrap(description, width = 25),
         cell = 'DLD-1') %>% 
  filter(description %in% categories)
  
  
lovo_cats = lovo_output %>% 
  mutate(direction=factor(direction, levels = c('UP', 'DOWN'))) %>% 
  mutate(FDR = cut(fdr, 
                   labels = c('0.001', '0.01', '0.05', 'ns'),
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf)),
         .before = term) %>% 
  mutate(description = str_wrap(description, width = 25),
         cell = 'LoVo') %>% 
  filter(description %in% categories)


merge_cats = hct_cats %>% 
  bind_rows(dld_cats) %>% 
  bind_rows(lovo_cats) 

# categories = unique(merge_cats$description)



highlight = function(x, pat, color="black") {
  ifelse(grepl(pat, x), 
         glue("<b style='font-size:18pt; color:{color}'>{x}</b>"), x)
}

merge_cats %>% 
  ggplot(aes(y = description, x = direction ,fill = FDR)) + 
  geom_tile() +
  scale_fill_manual(values = c('#2432FF',
                               '#616BFF',
                               '#A3A9FF',
                               '#FFFFFF')) +
  labs(
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits=rev,
                   labels= function(x) highlight(x, "p53 signaling pathway", "red")) +
  facet_wrap(~cell) +
  theme(
    strip.text = element_text(size = 15),
    axis.text.y = element_markdown(size = 12, lineheight = 0.6)
  ) +
  panel_border(color = 'black', size = 0.5)
  

merge_cats %>% 
  ggplot(aes(y = description, x = direction ,fill = FDR)) + 
  geom_tile() +
  scale_fill_manual(values = c('#2432FF',
                               '#616BFF',
                               '#A3A9FF',
                               '#FFFFFF')) +
  labs(
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~cell) +
  theme_cowplot(19) +
  panel_border() +
  theme(
    strip.text = element_text(size = 15),
    axis.text.y = element_markdown(size = 15, lineheight = 0.6),
    text = element_text(family="Arial")
  ) +
  panel_border(color = 'black', size = 0.5)

ggsave('summary/GSEA/poster_heatmap_enrich_KEGG.pdf',
       width = 12, height = 10)  

dev.copy2pdf(device = cairo_pdf,
             file = 'summary/GSEA/poster_heatmap_enrich_KEGG.pdf',
             width = 12, height = 10, useDingbats = FALSE)

merge_cats %>% 
  write_csv("summary/GSEA/enrich_data_paper.csv")

# paper version -----------------------------------------------------------

## filter the categories that show up at least twice in the 3 cell lines

library(extrafont)
# font_import()

selected_descriptions = merge_cats %>% 
  group_by(description) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(description)


merge_cats %>% 
  filter(description %in% selected_descriptions) %>% 
  ggplot(aes(y = description, x = direction ,fill = FDR)) + 
  geom_tile() +
  scale_fill_manual(values = c('#2432FF',
                               '#616BFF',
                               '#A3A9FF',
                               '#FFFFFF')) +
  labs(
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~cell) +
  theme_cowplot(15) +
  panel_border() +
  theme(
    strip.text = element_text(size = 15),
    axis.text.y = element_markdown(size = 12, lineheight = 0.6)
  ) +
  panel_border(color = 'black', size = 0.5)


ggsave('summary/GSEA/paper_heatmap_enrich_KEGG.pdf',
       width = 8, height = 6)  




# Including SW948 ---------------------------------------------------------



selected_descriptions = merge_cats %>%
  filter(cell %in% c("DLD-1", "HCT 116")) %>% 
  group_by(description) %>% 
  count() %>% 
  filter(n == 2) %>% 
  pull(description)


small_merge_cats = merge_cats %>% 
  filter(description %in% selected_descriptions) %>% 
  filter(cell %in% c("DLD-1", "HCT 116"))

small_merge_cats %>% 
  write_csv("summary/GSEA/small_enrich_paper.csv")

small_merge_cats = read_csv("summary/GSEA/small_enrich_paper.csv")

small_merge_cats %>% 
  mutate(cell = factor(cell, levels = c("HCT 116", "DLD-1", "SW948"))) %>% 
  mutate(FDR = as.factor(FDR)) %>% 
  ggplot(aes(y = description, x = direction ,fill = FDR)) + 
  geom_tile() +
  scale_fill_manual(values = c('#2432FF',
                               '#616BFF',
                               '#A3A9FF',
                               '#FFFFFF')) +
  labs(
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits = rev) +
  facet_wrap(~cell) +
  theme_cowplot(15, font_family = "Arial") +
  panel_border() +
  theme(
    strip.text = element_text(size = 15),
    axis.text.y = element_markdown(size = 12, lineheight = 0.6)
  ) +
  panel_border(color = 'black', size = 0.5) 


ggsave('summary/GSEA/paper_heatmap_enrich_KEGG_SW.pdf',
       width = 8, height = 6)  


# ENRICHMENT REVISION VARIANT --------------------------------------------------------



# selected_descriptions = merge_cats %>%
#   filter(cell %in% c("DLD-1", "HCT 116")) %>% 
#   group_by(description) %>% 
#   count() %>% 
#   filter(n == 2) %>% 
#   pull(description)


small_merge_cats = merge_cats %>% 
  filter(description %in% selected_descriptions) %>% 
  filter(cell %in% c("DLD-1", "HCT 116", "LoVo", "SW948"))


small_merge_cats %>% 
  write_csv("summary/GSEA/small_enrich_paper.csv")

small_merge_cats = read_csv("summary/GSEA/small_enrich_paper.csv")

small_merge_cats %>% 
  mutate(cell = factor(cell, levels = c("HCT 116", "DLD-1", "LoVo", "SW948"))) %>% 
  mutate(FDR = as.factor(FDR)) %>% 
  ggplot(aes(y = description, x = direction ,fill = FDR)) + 
  geom_tile() +
  scale_fill_manual(values = c('#2432FF',
                               '#616BFF',
                               '#A3A9FF',
                               '#FFFFFF')) +
  labs(
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits = rev) +
  facet_wrap(~cell, ncol = 4) +
  theme_cowplot(15, font_family = "Arial") +
  panel_border() +
  theme(
    strip.text = element_text(size = 15),
    axis.text.y = element_markdown(size = 12, lineheight = 0.6)
  ) +
  panel_border(color = 'black', size = 0.5) 


ggsave('summary/GSEA/paper_heatmap_enrich_KEGG_SW.pdf',
       width = 8, height = 6)  
ggsave('summary/GSEA/paper_heatmap_enrich_KEGG_REVISION.pdf',
       width = 9, height = 6)  



# complete version with SW ------------------------------------------------

merge_cats %>%
  filter(cell %in% c("DLD-1", "HCT 116")) %>% 
  filter(fdr < 0.05) %>% 
  write_csv("summary/GSEA/enrich_data_paper.csv")

merge_cats_alt = read_csv("summary/GSEA/enrich_data_paper.csv") %>% 
  mutate(FDR = as.factor(FDR))

merge_cats_alt %>% 
  mutate(cell = factor(cell, levels = c("HCT 116", "DLD-1", "SW948"))) %>% 
  ggplot(aes(y = description, x = direction ,fill = FDR)) + 
  geom_tile() +
  scale_fill_manual(values = c('#2432FF',
                               '#616BFF',
                               '#A3A9FF',
                               '#FFFFFF')) +
  labs(
    x = 'Regulation',
    y = 'KEGG Terms'
  ) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~cell) +
  theme_cowplot(15, font_family = "Arial") +
  panel_border() +
  theme(
    strip.text = element_text(size = 15),
    axis.text.y = element_markdown(size = 12, lineheight = 0.6)
  ) +
  panel_border(color = 'black', size = 0.5)


ggsave('summary/GSEA/paper_heatmap_enrich_KEGG_SW_complete.pdf',
       width = 8, height = 8)  



# P53 pathway -------------------------------------------------------------


p53_genes = c("CDK2","CDK4","CDK6","CDKN1A","CDKN2A","GADD45G","CHEK1",
              "CHEK2","SESN3","DDB2","GADD45A","RCHY1","BBC3","SESN1",
              "SFN","APAF1","IGF1","IGFBP3","FAS","CD82","MDM2","MDM4",
              "GADD45B","ATM","RRM2B","SERPINE1","SHISA5","GTSE1",
              "SERPINB5","PMAIP1","CYCS","ATR","STEAP3","PIDD1","RPRM",
              "PTEN","ADGRB1","BAX","CCND1","RRM2","BID","TP53AIP1","PERP",
              "COP1","ZMAT3","SIAH1",",THBS1","TP53","TP73","TSC2","CASP3",
              "SESN2","CASP8","CASP9","PPM1D","CCNB3","TNFRSF10B","CCNB1",
              "CCND2","CCND3","CCNE1","CCNG1","CCNG2","CCNB2","CCNE2","EI24",
              "TP53I3","CDK1")

#### Volcano plots ========

gene_counts %>% 
  filter(gene_name %in% p53_genes) %>% 
  filter(Cell_line != "SW948")


p53_stats = results.complete.kegg %>% 
  filter(Contrast != "SW948") %>% 
  filter(gene_name %in% p53_genes) %>% 
  distinct(gene_id, Contrast, .keep_all = T) %>% 
  mutate(Gene = case_when(padj < 0.05 ~ "Significant",
                          TRUE ~ "Not significant"))

cell = "DLD-1"
p53_stats %>% 
  filter(Contrast == cell) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = Gene)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not significant" = "black"
    )
  ) +
  ggrepel::geom_text_repel(aes(label = gene_name)) +
  facet_wrap(~Contrast) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot()

ggsave(glue("summary/P53/{cell}_volcano_p53path.pdf"), height = 7, width = 9)



p53_stats %>% 
  filter(Contrast == "HCT116") %>% 
  write_csv("summary/P53/HCT116_p53_stats.csv")


#### pairwise correlation =====

p53_stats %>% 
  filter(Contrast %in% c("HCT116", "DLD-1")) %>%
  select(Contrast, gene_name, log2FoldChange) %>%
  distinct(Contrast, gene_name, .keep_all = T) %>% 
  pivot_wider(names_from = Contrast, 
              values_from = log2FoldChange) %>% 
  ggplot(aes(x = HCT116, y = `DLD-1`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor() +
  labs(x = "HCT116 (log2FC)",
       y = "DLD-1 (log2FC)") +
  theme_cowplot()

ggsave('summary/P53/corr_HCT_DLD1.pdf', height = 7, width = 8)

p53_stats %>% 
  filter(Contrast %in% c("HCT116", "LoVo")) %>%
  select(Contrast, gene_name, log2FoldChange) %>%
  distinct(Contrast, gene_name, .keep_all = T) %>% 
  pivot_wider(names_from = Contrast, 
              values_from = log2FoldChange) %>% 
  ggplot(aes(x = HCT116, y = `LoVo`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor() +
  labs(x = "HCT116 (log2FC)",
       y = "LoVo (log2FC)") +
  theme_cowplot()

ggsave('summary/P53/corr_HCT_LoVo.pdf', height = 7, width = 8)

p53_stats %>% 
  filter(Contrast %in% c("LoVo", "DLD-1")) %>%
  select(Contrast, gene_name, log2FoldChange) %>%
  distinct(Contrast, gene_name, .keep_all = T) %>% 
  pivot_wider(names_from = Contrast, 
              values_from = log2FoldChange) %>% 
  filter(!(gene_name %in% c("SESN3", "CCNG2"))) %>% 
  ggplot(aes(x = LoVo, y = `DLD-1`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor() +
  labs(x = "LoVo (log2FC)",
       y = "DLD-1 (log2FC)") +
  theme_cowplot()

ggsave('summary/P53/corr_LoVo_DLD1.pdf', height = 7, width = 8)






#### PROTEOMICS VS RNASEQ -------------------------------------------------


p53_rna = p53_stats %>% 
  filter(Contrast == "HCT116")

p53_prot = read_csv("p53_proteins_micit.csv")


p53_rna = p53_rna %>% 
  select(gene_name, rna_log2fc = log2FoldChange) 

p53_prot = p53_prot %>% 
  select(gene_name = Gene_names, prot_log2fc = estimate) %>% 
  mutate(prot_log2fc = -1 * prot_log2fc) # fix directionality of comparison


p53_prot_rna = p53_rna %>% 
  left_join(p53_prot)



p53_prot_rna %>% 
  drop_na() %>% 
  ggplot(aes(x = rna_log2fc, y = prot_log2fc)) +
  geom_smooth(method = 'lm') +
  geom_point() +
  ggpubr::stat_cor() +
  labs(x = "RNAseq log2FC",
       y = "Proteomics log2FC") +
  theme_cowplot(font_family = "Arial")

ggsave("summary/P53/rna_prot_correlation.pdf", width = 6, height = 5)


p53_prot_rna %>% 
  mutate(prot_log2fc = replace_na(prot_log2fc , 0)) %>% 
  ggplot(aes(x = rna_log2fc, y = prot_log2fc)) +
  geom_smooth(method = 'lm') +
  geom_point() +
  ggpubr::stat_cor() +
  labs(x = "RNAseq log2FC",
       y = "Proteomics log2FC") +
  theme_cowplot(font_family = "Arial")



p53_prot_rna %>% 
  mutate(prot_log2fc = replace_na(prot_log2fc , 0)) %>% 
  ggplot(aes(x = rna_log2fc, 
             y =fct_reorder(gene_name, rna_log2fc), fill = prot_log2fc)) +
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = 'grey50', alpha = 0.7) +
  geom_point(shape=21, size = 3) +
  theme_cowplot(font_family = "Arial") +
  labs(y = "Genes/Proteins",
       x = "RNA-seq Log2FC",
       fill = "Proteomics\nLog2FC") +
  scale_fill_gradient2()

ggsave("summary/P53/rna_prot_dotplot.pdf", width = 6, height = 9)


# Volcano plots general ---------------------------------------------------



results.complete %>% 
  filter(Contrast == "HCT116") %>% 
  mutate(Gene = case_when(padj < 0.05 ~ "Significant",
                          TRUE ~ "Not significant"),
         labels = case_when(Gene == "Significant" ~ gene_name)) %>% 
  filter(log2FoldChange < 3) %>% # filter an outlier
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = Gene)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not significant" = "black"
    )
  ) +
  xlim(-2, 2.5) + ylim(0, 180) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 10,
  #                          size = 3) +
  facet_wrap(~Contrast) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(14, font_family = "Arial")


ggsave(glue("summary/volcano_HCT116.pdf"), height = 7, width = 9)



results.complete %>% 
  filter(Contrast == "DLD-1") %>% 
  mutate(Gene = case_when(padj < 0.05 ~ "Significant",
                          TRUE ~ "Not significant"),
         labels = case_when(Gene == "Significant" ~ gene_name)) %>% 
  # filter(log2FoldChange < 3) %>% # filter an outlier
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = Gene)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not significant" = "black"
    )
  ) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 10,
  #                          size = 3) +
  facet_wrap(~Contrast) +
  xlim(-2, 2.5) + ylim(0, 180) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(12, font_family = "Arial")


ggsave(glue("summary/volcano_DLD1.pdf"), height = 7, width = 9)



results.complete %>% 
  filter(Contrast == "LoVo") %>% 
  mutate(Gene = case_when(padj < 0.05 ~ "Significant",
                          TRUE ~ "Not significant"),
         labels = case_when(Gene == "Significant" ~ gene_name)) %>% 
  # filter(log2FoldChange < 3) %>% # filter an outlier
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = Gene)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not significant" = "black"
    )
  ) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 15,
  #                          size = 3) +
  facet_wrap(~Contrast) +
  xlim(-2, 2.5) + ylim(0, 180) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(12, font_family = "Arial")


ggsave(glue("summary/volcano_LoVo.pdf"), height = 7, width = 9)




results.complete %>% 
  filter(Contrast == "SW948") %>% 
  mutate(Gene = case_when(padj < 0.05 ~ "Significant",
                          TRUE ~ "Not significant"),
         labels = case_when(Gene == "Significant" ~ gene_name)) %>% 
  # filter(log2FoldChange < 3) %>% # filter an outlier
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = Gene)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not significant" = "black"
    )
  ) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 15,
  #                          size = 3) +
  facet_wrap(~Contrast) +
  # xlim(-2, 4.5) + 
  # ylim(0, 180) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(12, font_family = "Arial")


ggsave(glue("summary/volcano_SW948_scale.pdf"), height = 7, width = 9)
