### RNA seq analysis of human cell lines from Tanara

# In this script we will analyse the RNA seq with the following cell lines:
# 	- HT29
# 	- SK-CO-1
# 	- SW1417

# useful links:
# 	- https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
# 	- https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# 	- https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

# maybe use ulimit -s 16384 before start R console

### libraries ####
library(tximport)
library(tidyverse)
library(DESeq2)
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
library(readxl)
# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2

extrafont::loadfonts()

theme_set(theme_cowplot(14, font_family = "Arial"))


# # # # # # # # # # #
# RNA seq analysis --------------------------------------------------------
# # # # # # # # # # #


# Get sample info ---------------------------------------------------------



samples = read.delim("sampleInfo.txt") 

dir = getwd()
rownames(samples) = samples$Name

quants_dir = 'quants_103'


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


# read.samples(samp_file = 'sampleInfo_HCT116.txt')

# Get gene labels from databases ------------------------------------------

#### using ensembl genomes ####

# let's make our database from ensembldb
ah = AnnotationHub::AnnotationHub()
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


# samples.batch = samples %>% 
#   mutate(Batch = as.factor(Batch),
#          Condition = as.factor(Condition),
#          Sample = as.factor(Sample))
# 
# ddsTxi = DESeqDataSetFromTximport(txi, colData = samples.batch, 
#                                   design = ~Batch + Sample)

samples = samples %>% 
  mutate(Sample = as.factor(Sample),
         Batch = as.factor(Batch),
         Condition = as.factor(Condition),
         Cell_line = as.factor(Cell_line))

ddsTxi = DESeqDataSetFromTximport(txi, colData = samples,
                                  design = ~Sample)



#### filter by rowSums ####
# prefilter, but that might not be necessary

# histogram of rowSums(counts(ddsTxi))
hist(rowSums(counts(ddsTxi)), 
     breaks = 1000, col = 'blue', 
     main = 'Histogram of rowSums(counts(ddsTxi))')

# histogram of rowSums(counts(ddsTxi)) below 10000
hist(rowSums(counts(ddsTxi)), 
     breaks = 100000, col = 'blue', 
     main = 'Histogram of rowSums(counts(ddsTxi)) below 10000',
     xlim = c(0, 5000))

keep = rowSums(counts(ddsTxi)) >= 100



ddsTxi = ddsTxi[keep,]


# Run DESeq2 --------------------------------------------------------------

# run the Differential Expression Analysis
dds = DESeq(ddsTxi)

detach("package:ensembldb", unload = TRUE)


dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)


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
  gather(Name, counts, HT29_C_1:SW1417_M_4) %>% 
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

gene = 'ENSG00000139517'
gene_counts%>% 
  dplyr::filter(gene_id == gene) %>% 
  # dplyr::filter(Cell_line == 'HCT116') %>%
  ggplot(aes(x = Sample, y = counts, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, position = position_jitterdodge())+
  facet_wrap(~gene_name*Cell_line, scales = 'free') +
  scale_fill_manual(values = c("#1C86EE", "#EEC900")) +
  labs(x = 'Sample',
       y = 'Counts (normalised)')




# transofrm data ------------------------

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



# PCA per cell line -------------------------------------------------------

#### HT29 -----------------------
files_ht29 = read.samples(samp_file = 'sampleInfo_HT29.txt')

# import quantification data
txi_ht29 = tximport::tximport(files_ht29, type = "salmon", tx2gene = tx2gene)

samples_ht29 = samples %>% filter(Cell_line == 'HT29')

ddsTxi_ht29 = DESeqDataSetFromTximport(txi_ht29, colData = samples_ht29, 
                                      design = ~Sample)

#filter by min number of sequences
keep = rowSums(counts(ddsTxi_ht29)) >= 40
ddsTxi_ht29 = ddsTxi_ht29[keep,]

dds_ht29 = DESeq(ddsTxi_ht29)

# rld
rld_ht29 = rlog(dds_ht29, blind = FALSE)

pcaData_ht29 = plotPCA(rld_ht29, intgroup = c("Condition"), returnData = TRUE)

# get ellipses
ell_ht29 = pcaData_ht29 %>% group_by(Condition) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar_ht29 = round(100 * attr(pcaData_ht29, "percentVar"))

# plot!
ggplot(pcaData_ht29, aes(x = PC1, y = PC2, color = Condition, group = Condition)) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell_ht29, aes(x = x, y = y, group = Condition, linetype = Condition), size = 1) +
  geom_polygon(data = ell_ht29, aes(x = x, y = y, group = Condition, 
                                    linetype = Condition, fill = Condition), size = 1, alpha = 0.1) +
  xlab(paste0("PC1: ", percentVar_ht29[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_ht29[2], "% variance")) +
  theme_classic() +
  ggrepel::geom_text_repel(aes(label = name), size = 3) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black'))

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_HT29_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)


#### SW1417 -----------------------
files_sw1417 = read.samples(samp_file = 'sampleInfo_SW1417.txt')

# import quantification data
txi_sw1417 = tximport::tximport(files_sw1417, type = "salmon", tx2gene = tx2gene)

samples_sw1417 = samples %>% filter(Cell_line == 'SW1417')

ddsTxi_sw1417 = DESeqDataSetFromTximport(txi_sw1417, colData = samples_sw1417, 
                                      design = ~Sample)

#filter by min number of sequences
keep = rowSums(counts(ddsTxi_sw1417)) >= 40
ddsTxi_sw1417 = ddsTxi_sw1417[keep,]

dds_sw1417 = DESeq(ddsTxi_sw1417)

# rld
rld_sw1417 = rlog(dds_sw1417, blind = FALSE)

pcaData_sw1417 = plotPCA(rld_sw1417, intgroup = c("Condition"), returnData = TRUE)

# get ellipses
ell_sw1417 = pcaData_sw1417 %>% group_by(Condition) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar_sw1417 = round(100 * attr(pcaData_sw1417, "percentVar"))

# plot!
ggplot(pcaData_sw1417, aes(x = PC1, y = PC2, color = Condition, group = Condition)) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell_sw1417, aes(x = x, y = y, group = Condition, linetype = Condition), size = 1) +
  geom_polygon(data = ell_sw1417, aes(x = x, y = y, group = Condition, 
                                    linetype = Condition, fill = Condition), size = 1, alpha = 0.1) +
  xlab(paste0("PC1: ", percentVar_sw1417[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_sw1417[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black'))

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_SW1417_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)


#### SK-CO-1 -----------------------

files_skco1 = read.samples(samp_file = 'sampleInfo_SKCO1.txt')

# import quantification data
txi_skco1 = tximport::tximport(files_skco1, type = "salmon", tx2gene = tx2gene)

samples_skco1 = samples %>% filter(Cell_line == 'SKCO1')

ddsTxi_skco1 = DESeqDataSetFromTximport(txi_skco1, colData = samples_skco1, 
                                      design = ~Sample)

#filter by min number of sequences
keep = rowSums(counts(ddsTxi_skco1)) >= 40
ddsTxi_skco1 = ddsTxi_skco1[keep,]

dds_skco1 = DESeq(ddsTxi_skco1)

# rld

rld_skco1 = rlog(dds_skco1, blind = FALSE)

pcaData_skco1 = plotPCA(rld_skco1, intgroup = c("Condition"), returnData = TRUE)

# get ellipses
ell_skco1 = pcaData_skco1 %>% group_by(Condition) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar_skco1 = round(100 * attr(pcaData_skco1, "percentVar"))

# plot!
ggplot(pcaData_skco1, aes(x = PC1, y = PC2, color = Condition, group = Condition)) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell_skco1, aes(x = x, y = y, group = Condition, linetype = Condition), size = 1) +
  geom_polygon(data = ell_skco1, aes(x = x, y = y, group = Condition, 
                                    linetype = Condition, fill = Condition), size = 1, alpha = 0.1) +
  xlab(paste0("PC1: ", percentVar_skco1[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_skco1[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black'))

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_SKCO1_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)

# KEGG databases -----------------------------------------------

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
res = results(dds, pAdjustMethod = 'fdr') 


# results with different shape of contrasts, tidy
res.ht29 = results(dds,   contrast = c("Sample", "HT29_M" , "HT29_C"),
                   pAdjustMethod = 'fdr')  
res.ht29 = lfcShrink(dds, contrast = c("Sample", "HT29_M" , "HT29_C"), 
                     res = res.ht29, type = 'ashr')

res.sw1417 = results(dds,  contrast = c("Sample",  "SW1417_M", "SW1417_C"),
                     pAdjustMethod = 'fdr') 
res.sw1417 = lfcShrink(dds, contrast = c("Sample",  "SW1417_M", "SW1417_C"), 
                       res = res.sw1417, type = 'ashr')

res.skco1 = results(dds,  contrast = c("Sample",  "SKCO1_M", "SKCO1_C"),
                    pAdjustMethod = 'fdr')   
res.skco1 = lfcShrink(dds, contrast = c("Sample", "SKCO1_M", "SKCO1_C"), 
                     res = res.skco1, type = 'ashr')

#### tidying results ####
res.ht29.tidy = as_tibble(res.ht29, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'HT29',
  Contrast_description = 'Comparison of HT29 Micit vs Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.sw1417.tidy = as_tibble(res.sw1417, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'SW1417',
  Contrast_description = 'Comparison of SW1417 Micit vs Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.skco1.tidy = as_tibble(res.skco1, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'SK-CO-1',
  Contrast_description = 'Comparison of SK-CO-1 Micit vs Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())



results.complete = res.ht29.tidy %>% rbind(res.sw1417.tidy, res.skco1.tidy)

results.complete.kegg = results.complete %>% left_join(kegg.links)

# write results in excel files
list_of_datasets = list('HT29' = res.ht29.tidy, 
                        'SW1417' = res.sw1417.tidy, 
                        'SK-CO-1' = res.skco1.tidy)

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
plotMA(res.ht29,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_HT29.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.sw1417,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_SW1417.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.skco1,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_SKCO1.pdf'),
             height = 8, width = 11, useDingbats = FALSE)



# GSEA --------------------------------------------------------------------



res.ht29.tidy %>% 
  dplyr::filter(str_detect(gene_name,'CPT'))

gene = 'ZMYM5'
gene_counts %>% 
  dplyr::filter(gene_name == gene) %>% 
  dplyr::filter(Cell_line == 'HT29') %>%
  ggplot(aes(x = Sample, y = counts, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, position = position_jitterdodge(), aes(color = Replicate))+
  facet_wrap(~gene_id, scales = 'free_y')




# DE genes (up/down) to analyse via StringDB

# helper function that extracts DE genes with upper and lower thresholds
DEgenes = function(dataset = results.complete, 
                   cell = 'HT29', 
                   min_thrs = 0,
                   max_thrs = 10,
                   padj_thrs = 0.05) {
  
  up = dataset %>%  
    filter(Contrast == cell) %>% 
    filter(log2FoldChange < max_thrs & log2FoldChange > min_thrs) %>% 
    filter(padj <= padj_thrs) %>%
    pull(gene_id) %>% unique
  
  
  down = dataset %>% 
    filter(Contrast == cell) %>% 
    filter(log2FoldChange > -max_thrs & log2FoldChange < -min_thrs) %>% 
    filter(padj <= padj_thrs) %>%
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

# HT29
ht29_de = DEgenes(results.complete, cell = 'HT29', 
                  min_thrs = 0, max_thrs = 10,
                  padj_thrs = 0.05)
# SW1417
sw1417_de = DEgenes(results.complete, cell = 'SW1417', 
                    min_thrs = 0, max_thrs = 10,
                    padj_thrs = 0.05)
# SKCO1
skco1_de = DEgenes(results.complete, cell = 'SK-CO-1', 
                   min_thrs = 0, max_thrs = 10,
                   padj_thrs = 0.05)


list_of_datasets = list('HT29_UP' = ht29_de[[1]],
                        'HT29_DOWN'= ht29_de[[2]],
                        'SW1417_UP' = sw1417_de[[1]],
                        'SW1417_DOWN'= sw1417_de[[2]],
                        'SKCO1_UP' = skco1_de[[1]],
                        'SKCO1_DOWN'= skco1_de[[2]]
)

write.xlsx(list_of_datasets, here('summary', 'total_genes_updown.xlsx'),overwrite = TRUE)




# boxplots ----------------------------------------------------------------

geneplot = function(gene = "ENSG00000259746"){
  gene_counts %>%
    # mutate(Replicate = as.factor(Replicate)) %>% 
    dplyr::filter(gene_id == gene) %>%
    # dplyr::filter(Cell_line == 'HCT116') %>%
    ggplot(aes(x = Sample, y = counts, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(
      aes(color = Replicate),
      size = 2.5, position = position_jitterdodge())+
    facet_wrap(~gene_name*Cell_line, scales = 'free') +
    scale_fill_manual(values = c("#1C86EE", "#EEC900")) +
    labs(x = 'Sample',
         y = 'Counts (normalised)')
}

geneplot("ENSG00000113722")


### plot gene boxplots ####


# merge counts with stats info
results.merge = gene_counts %>% 
  left_join(results.complete %>%
              mutate(Contrast = case_when(Contrast == "SK-CO-1" ~ "SKCO1",
                                          TRUE ~ Contrast)) %>% 
              dplyr::select(gene_id, Cell_line = Contrast, log2FoldChange, 
                     lfcSE, pvalue, padj, p_adj_stars)) 


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
  
  ggsave(here('summary/ALL_gene_boxplots', glue('{gene}_boxplot.pdf')), height = 8, width = 13)
  
}






# Volcano plots general ---------------------------------------------------



results.complete %>% 
  filter(Contrast == "HT29") %>% 
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
  xlim(-4, 3) +
  ylim(0, 180) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 10,
  #                          size = 3) +
  facet_wrap(~Contrast) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(14, font_family = "Arial")


ggsave(glue("summary/volcano_HT29.pdf"), height = 7, width = 9)

results.complete %>% 
  filter(Contrast == "SK-CO-1") %>% 
  drop_na(padj) %>% 
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
  xlim(-2, 3) +
  ylim(0, 180) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 10,
  #                          size = 3) +
  facet_wrap(~Contrast) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(14, font_family = "Arial")


ggsave(glue("summary/volcano_SKCO1.pdf"), height = 7, width = 9)



results.complete %>% 
  filter(Contrast == "SW1417") %>% 
  drop_na(padj) %>% 
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
  xlim(-2, 2) +
  ylim(0, 180) +
  # ggrepel::geom_text_repel(aes(label = labels),
  #                          box.padding = 0,
  #                          max.overlaps = 10,
  #                          size = 3) +
  facet_wrap(~Contrast) +
  labs(x = "log2(FC)",
       y = "-log10(P-value adj)") +
  theme_cowplot(14, font_family = "Arial")


ggsave(glue("summary/volcano_SW1417.pdf"), height = 7, width = 9)




# STRING enrichment plots -------------------------------------------------

# read data
# HT29
library(readxl)

ht29.enrich = read_excel("String_analysis/HT29/HT29_output.xlsx", sheet = "KEGG") %>% 
  dplyr::select(-`...1`) %>% mutate(Cell_line = "HT29")

sw1417.enrich = read_excel("String_analysis/SW1417/SW1417_output.xlsx", sheet = "KEGG") %>% 
  dplyr::select(-`...1`) %>% mutate(Cell_line = "SW1417")

skco1.enrich = read_excel("String_analysis/SKCO1/SKCO1_output.xlsx", sheet = "KEGG") %>% 
  dplyr::select(-`...1`) %>% mutate(Cell_line = "SKCO1")

# keep paths that show up at least twice
selected_kegg_enrich = ht29.enrich %>% 
  bind_rows(sw1417.enrich) %>%
  bind_rows(skco1.enrich) %>%
  dplyr::count(description) %>% 
  filter(n > 1) %>% 
  pull(description)


ht29.enrich %>% 
  filter((description %in% selected_kegg_enrich)) %>% 
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


ht29.enrich %>% 
  bind_rows(sw1417.enrich) %>%
  bind_rows(skco1.enrich) %>% 
  filter((description %in% selected_kegg_enrich)) %>%
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
    y = NULL
  ) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~Cell_line)

ggsave("summary/GSEA/KEGG_terms.pdf", height = 8, width = 10)


p53_ht29 = ht29.enrich %>% 
  filter(description == "p53 signaling pathway") %>% 
  mutate(preferredNames = str_replace(preferredNames, "\\[", ''),
         preferredNames = str_replace(preferredNames, "\\]", ""),
         preferredNames = str_replace_all(preferredNames, "'", "")) %>% 
  separate_rows(preferredNames, sep = ', ') %>% 
  pull(preferredNames)


p53_skco = skco1.enrich %>% 
  filter(description == "p53 signaling pathway") %>% 
  mutate(preferredNames = str_replace(preferredNames, "\\[", ''),
         preferredNames = str_replace(preferredNames, "\\]", ""),
         preferredNames = str_replace_all(preferredNames, "'", "")) %>% 
  separate_rows(preferredNames, sep = ', ') %>% 
  pull(preferredNames)

intersect(p53_ht29, p53_skco)


p53_ht29
p53_skco




# HTML Gene enrichment -------------------------------------------------

# Enrichment Browser

library(EnrichmentBrowser)


## Run this once!

# obtaining gene sets
kegg.gs = getGenesets(org = "hsa", db = "kegg")
go.gs = getGenesets(org = "hsa", db = "go", onto = "BP")





#### HT29 ####
# need to do again DESeq2 for hct per separate
files.ht = read.samples(samp_file = 'sampleInfo_HT29.txt')
# import quantification data 
txi.ht = tximport::tximport(files.ht, type = "salmon", tx2gene = tx2gene)

samples.red = samples %>% filter(Cell_line == 'HT29')
ddsTxi.ht = DESeqDataSetFromTximport(txi.ht, 
                                      colData = samples.red, 
                                      design = ~Sample)

# filtering by min number of sequences
keep = rowSums(counts(ddsTxi.ht)) >= 40
ddsTxi.ht = ddsTxi.ht[keep,]

ddsTxi.ht$Sample = relevel(ddsTxi.ht$Sample, ref = "HT29_C")


# run the Differential Expression Analysis
dds.ht = DESeq(ddsTxi.ht)

res.ht.pure = results(dds.ht)

# import for gene enrichment
ht.SE = import(dds.ht, res.ht.pure, from = c('DESeq2'), anno = 'hsa')

# map IDs
ht.SE = idMap(ht.SE, org = "hsa", from = "ENSEMBL", to = "ENTREZID")


head(rownames(ht.SE))

# set based enrichment analysis
sbeaMethods()

# normalize counts
ht.SE = normalize(ht.SE, norm.method = "vst")

# run GSEA enrichment
# kegg 
alpha = 0.05
ht.gsea = sbea(method = "gsea", se = ht.SE, gs = kegg.gs, alpha = alpha)
ht.ora = sbea(method = "ora", se = ht.SE, gs = kegg.gs, alpha = alpha)
ht.padog = sbea(method = "padog", se = ht.SE, gs = kegg.gs, alpha = alpha)

ht.go.gsea = sbea(method = "gsea", se = ht.SE, gs = go.gs, alpha = 0.05)
gsRanking(ht.go.gsea)

# 
# gsRanking(hct.gsea)
# gsRanking(hct.ora)

eaBrowse(ht.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/ht.gsea', 
         report.name = 'ht.gsea')

eaBrowse(ht.ora, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/ht.ora', 
         report.name = 'ht.ora')

eaBrowse(ht.padog, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/ht.padog', 
         report.name = 'ht.padog')


eaBrowse(ht.go.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/GO_BP/ht.go.gsea', 
         report.name = 'ht.go.gsea')

# network regulation analysis

hsa.grn = compileGRN(org="hsa", db="kegg")


nbeaMethods()

nbea.res = nbea(method="ggea", se=ht.gsea, gs=kegg.gs, grn=hsa.grn)

gsRanking(nbea.res)






# Dorothea analysis -----------------------------------------------------------



## We load the required packages
library(dorothea)
library(bcellViper)
library(viper)


# Running VIPER with DoRothEA regulons
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# get the gene counts from HT29 cell line
ht_counts = gene_counts %>% 
  filter(Cell_line == 'HT29') %>% 
  dplyr::select(gene_name, Sample, Replicate, counts) %>% 
  unite(Sample, Sample, Replicate) %>% 
  pivot_wider(names_from = Sample, values_from = counts,
              values_fn = {mean}) %>% # there are a few repeated elements
  data.frame()

# rownames
row.names(ht_counts) = ht_counts[,1]
ht_counts[,1] = NULL

# run viper wrapper with dorothea
tf_ht = run_viper(ht_counts, regulons, tidy = TRUE,
                   options =  list(method = "scale", minsize = 4, 
                                   eset.filter = FALSE, cores = 1, 
                                   verbose = T))

# tidy the results
tidy_tf = tf_ht %>% 
  separate(sample , into = c('Cell', 'Condition', 'Replicate'), sep = '_') %>% 
  as_tibble() %>% 
  group_by(tf, Cell, Condition, confidence) %>% 
  summarise(mean_activity = mean(activity, na.rm = TRUE),
            sd_activity = sd(activity, na.rm = TRUE)) %>% 
  arrange(desc(abs(mean_activity)))  %>% 
  ungroup()


# calculate the stats: Treatment vs Control
tf_ht_stats = tf_ht %>% 
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
sig_tf = tf_ht_stats  %>% 
  filter(p < 0.05) %>% 
  pull(tf)



#### dot plots #####

# print selection of TFs

tidy_tf %>% 
  left_join(tf_ht_stats) %>% 
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













# FULL DATASET ------------------------------------------------------------

samples_full = read.delim('sampleInfo_FULL.txt')


dir = getwd()
rownames(samples_full) = samples_full$Name

quants_dir = 'quants_103'


# prepare a list with file names
files = file.path(dir,quants_dir, samples_full$Name, "quant.sf")
names(files) = samples_full$Name
all(file.exists(files)) # check that files exist

### TX import + DESeq ####

# import quantification data 
txi = tximport::tximport(files, type = "salmon", tx2gene = tx2gene)


samples_full = samples_full %>% 
  mutate(Sample = as.factor(Sample),
         Batch = as.factor(Batch),
         Condition = as.factor(Condition),
         Cell_line = as.factor(Cell_line))

ddsTxi = DESeqDataSetFromTximport(txi, 
                                  colData = samples_full,
                                  design = ~Sample)



#### filter by rowSums ####
keep = rowSums(counts(ddsTxi)) >= 100


ddsTxi = ddsTxi[keep,]


# Run DESeq2 --------------------------------------------------------------

# run the Differential Expression Analysis
dds_full = DESeq(ddsTxi)

dds_full <- estimateSizeFactors(dds_full)
dds_full <- estimateDispersions(dds_full)
plotDispEsts(dds_full)


# transform the data
rld = rlog(dds_full, blind = FALSE)



# PCA ---------------------------------------------------------------------


pcaData = plotPCA(rld, intgroup = c("Cell_line", "Condition"), returnData = TRUE)
pcaData

pcaData %>% 
  as_tibble() %>% 
  select(-name) %>% 
  write.csv("summary/pca_data.csv")
 
# get info for the ellipses
ell = pcaData %>% group_by(Condition, Cell_line) %>% 
  do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame 

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
pcaData %>% 
  # filter(Cell_line != "LoVo") %>% 
  ggplot(aes(x = PC1, y = PC2, color = Cell_line, group = interaction(Condition, Cell_line))) + 
  geom_point(size = 4, show.legend = NA, alpha = 1,
             aes(shape = Condition)) + 
  # geom_path(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), linetype = Condition), size = 1) +
  # geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), 
  #                              linetype = Condition, fill = Cell_line), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black')) 


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_ALL_CELLS_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)





### tidy results ####---------
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts_full = counts(dds_full, normalized = TRUE)
gene_list = rownames(gene_counts_full)
gene_counts_full = gene_counts_full %>% cbind(gene_list,.) %>% as_tibble()

gene_counts_full = gene_counts_full %>% 
  gather(Name, counts, HT29_C_1:TVP_32) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(as_tibble(samples_full), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) 

gene_counts_full %>% 
  write_csv(here('summary', 'FULL_version', 'gene_counts_norm.csv'))



## STATS ------------
# load stats from the other cell lines and merge all into one table 

stats_Tanara_hct = read_excel("~/Documents/MRC_postdoc/Cancer_5FU/RNA-seq/summary/complete_stats.xlsx",
                              sheet = "HCT116")

stats_Tanara_dld = read_excel("~/Documents/MRC_postdoc/Cancer_5FU/RNA-seq/summary/complete_stats.xlsx",
                              sheet = "DLD-1")

stats_Tanara_lovo = read_excel("~/Documents/MRC_postdoc/Cancer_5FU/RNA-seq/summary/complete_stats.xlsx",
                               sheet = "LoVo")

stats_Tanara_sw948 = read_excel("~/Documents/MRC_postdoc/Cancer_5FU/RNA-seq/summary/complete_stats.xlsx",
                              sheet = "SW948")



stats.full = results.complete %>% 
  bind_rows(stats_Tanara_hct) %>% 
  bind_rows(stats_Tanara_dld) %>%
  bind_rows(stats_Tanara_lovo) %>%
  bind_rows(stats_Tanara_sw948) 



# FULL HEATMAPS -----------------------------------------------------------

#### Pyrimidine metabolism ------------------------

pyrimidine_genes = kegg.links %>% 
  filter(str_detect(Description,"Pyrimidine")) %>% 
  distinct(entrezid) %>% pull(entrezid)

pyrimidine_genes_sig = stats.full %>% 
  filter(entrezid %in% pyrimidine_genes) %>% 
  filter(padj < 0.05) %>% 
  dplyr::count(entrezid) %>% 
  filter(n >= 1) %>% 
  pull(entrezid)


# plot heatmap
stats.full %>% 
  mutate(p_adj_stars = str_replace(p_adj_stars, '.', '')) %>% 
  filter(entrezid %in% pyrimidine_genes) %>% 
  mutate(Contrast = factor(Contrast, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                       'HT29', 'SK-CO-1', 'SW1417', 
                                      'SW948'))) %>%
  ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
  scale_fill_gradient2(low = "blue",
                       mid = "grey99",
                       high = "red",
                       midpoint = 0) +
  geom_tile() +
  geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
  labs(x = NULL,
       y = NULL,
       fill = 'Fold Change\n(log2)') +
  theme(axis.text.y = NULL)

ggsave(here('summary','FULL_version','heatmaps', 
            'pyrimidine_full_heatmap.pdf'), height = 8, width = 10)


# sig genes heatmap
stats.full %>% 
  mutate(p_adj_stars = str_replace(p_adj_stars, '.', '')) %>% 
  filter(entrezid %in% pyrimidine_genes_sig) %>% 
  mutate(Contrast = factor(Contrast, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                      'HT29', 'SK-CO-1', 'SW1417', 
                                      'SW948'))) %>%
  ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
  scale_fill_gradient2(low = "blue",
                       mid = "grey99",
                       high = "red",
                       midpoint = 0) +
  geom_tile() +
  geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
  labs(x = NULL,
       y = NULL,
       fill = 'Fold Change\n(log2)') +
  theme(axis.text.y = NULL)

ggsave(here('summary','FULL_version','heatmaps', 
            'pyrimidine_sig_heatmap.pdf'), height = 5, width = 7)


# sig genes heatmap
stats.full %>% 
  mutate(p_adj_stars = str_replace(p_adj_stars, '.', '')) %>% 
  filter(entrezid %in% pyrimidine_genes_sig) %>% 
  mutate(Contrast = factor(Contrast, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                      'HT29', 'SK-CO-1', 'SW1417', 
                                      'SW948'))) %>%
  ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
  scale_fill_gradient2(low = "blue",
                       mid = "grey99",
                       high = "red",
                       midpoint = 0) +
  geom_tile() +
  # geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
  labs(x = NULL,
       y = NULL,
       fill = 'Fold Change\n(log2)') +
  theme(axis.text.y = NULL)

ggsave(here('summary','FULL_version','heatmaps', 
            'pyrimidine_sig_nostats_heatmap.pdf'), height = 5, width = 7)


#### purine metabolism ------------------------

purine_genes = kegg.links %>% 
  filter(str_detect(Description,"Purine")) %>% 
  distinct(entrezid) %>% pull(entrezid)

purine_genes_sig = stats.full %>%
  filter(entrezid %in% purine_genes) %>% 
  filter(padj < 0.05) %>% 
  dplyr::count(entrezid) %>% 
  filter(n >= 1) %>% 
  pull(entrezid)


# plot heatmap

stats.full %>% 
  mutate(p_adj_stars = str_replace(p_adj_stars, '.', '')) %>% 
  filter(entrezid %in% purine_genes) %>% 
  filter(log2FoldChange < 3) %>% 
  mutate(Contrast = factor(Contrast, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                      'HT29', 'SK-CO-1', 'SW1417', 
                                      'SW948'))) %>%
  ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
  scale_fill_gradient2(low = "blue",
                       mid = "grey99",
                       high = "red",
                       midpoint = 0) +
  geom_tile() +
  geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
  labs(x = NULL,
       y = NULL,
       fill = 'Fold Change\n(log2)') +
  theme(axis.text.y = NULL)

ggsave(here('summary','FULL_version','heatmaps', 
            'purine_full_heatmap.pdf'), height = 8, width = 10)


# sig genes heatmap

stats.full %>% 
  mutate(p_adj_stars = str_replace(p_adj_stars, '.', '')) %>% 
  filter(entrezid %in% purine_genes_sig) %>% 
  filter(log2FoldChange < 3) %>% 
  mutate(Contrast = factor(Contrast, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                      'HT29', 'SK-CO-1', 'SW1417', 
                                      'SW948'))) %>%
  ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
  scale_fill_gradient2(low = "blue",
                       mid = "grey99",
                       high = "red",
                       midpoint = 0) +
  geom_tile() +
  geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
  labs(x = NULL,
       y = NULL,
       fill = 'Fold Change\n(log2)') +
  theme(axis.text.y = NULL)

ggsave(here('summary','FULL_version','heatmaps', 
            'purine_sig_heatmap.pdf'), height = 8, width = 7)


# sig genes heatmap

stats.full %>% 
  mutate(p_adj_stars = str_replace(p_adj_stars, '.', '')) %>% 
  filter(entrezid %in% purine_genes_sig) %>% 
  filter(log2FoldChange < 3) %>% 
  mutate(Contrast = factor(Contrast, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                      'HT29', 'SK-CO-1', 'SW1417', 
                                      'SW948'))) %>%
  ggplot(aes(y = reorder(gene_name, log2FoldChange), x = Contrast, fill = log2FoldChange)) +
  scale_fill_gradient2(low = "blue",
                       mid = "grey99",
                       high = "red",
                       midpoint = 0) +
  geom_tile() +
  # geom_text(aes(label = p_adj_stars),nudge_y = -0.2) +
  labs(x = NULL,
       y = NULL,
       fill = 'Fold Change\n(log2)') +
  theme(axis.text.y = NULL)

ggsave(here('summary','FULL_version','heatmaps', 
            'purine_sig_nostats_heatmap.pdf'), height = 8, width = 7)




# shared enrichment  ------------------------------------------------------


library(readxl)

ht29.enrich 
sw1417.enrich 
skco1.enrich


hct.enrich = read_excel("String_analysis/HCT116/HCT116_output.xlsx", sheet = "KEGG") %>% 
  dplyr::select(-`...1`) %>% mutate(Cell_line = "HCT116")

dld.enrich = read_excel("String_analysis/DLD/DLD_output.xlsx", sheet = "KEGG") %>% 
  dplyr::select(-`...1`) %>% mutate(Cell_line = "DLD-1")

lovo.enrich = read_excel("String_analysis/LoVo/LoVo_output.xlsx", sheet = "KEGG") %>% 
  dplyr::select(-`...1`) %>% mutate(Cell_line = "LoVo")

# keep paths that show up at least twice
selected_kegg_enrich =
  ht29.enrich %>% 
  bind_rows(sw1417.enrich) %>%
  bind_rows(skco1.enrich) %>%
  bind_rows(hct.enrich) %>%
  bind_rows(dld.enrich) %>%
  bind_rows(lovo.enrich) %>%
  dplyr::count(description) %>% 
  filter(n >= 3) %>% 
  filter(!str_detect(description, "disease")) %>% 
  arrange(description) %>% 
  pull(description)




full.enrich = ht29.enrich %>% 
  bind_rows(sw1417.enrich) %>%
  bind_rows(skco1.enrich) %>% 
  bind_rows(hct.enrich) %>%
  bind_rows(dld.enrich) %>%
  bind_rows(lovo.enrich)
  
# create a dummy line for SW948

sw948.enrich = sw1417.enrich %>% 
  mutate(Cell_line = "SW948") %>% 
  mutate(description = str_replace(description, "SW1417", "SW948")) %>% 
  filter(term == "hsa03050") %>% 
  mutate(fdr = 1, p_value = 1)

full.enrich = full.enrich %>% 
  bind_rows(sw948.enrich) %>% 
  mutate(Cell_line = factor(Cell_line, 
                           levels = c('HCT116', 'DLD-1', 'LoVo', 
                                      'HT29', 'SKCO1', 'SW1417', 
                                      'SW948')))


full.enrich %>%
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
  facet_wrap(~Cell_line, nrow = 1) +
  panel_border()

ggsave("summary/FULL_version/enrichment/KEGG_terms_full.pdf", height = 10, width = 10)



## clean a bit the categories

cat_removals = c("Amyotrophic lateral sclerosis", "Bladder cancer",
                 "Bacterial invasion of epithelial cells", "Cellular senescence",
                 "Epstein-Barr virus infection", "Human papillomavirus infection",
                 "MicroRNAs in cancer",
                 "Pathogenic Escherichia coli infection", "Proteoglycans in cancer",
                 "Salmonella infection", "Tight junction")

full.enrich %>% 
  filter(!description %in% cat_removals) %>%
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
    y = NULL
  ) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~Cell_line, nrow = 1) +
  panel_border()

ggsave("summary/FULL_version/enrichment/KEGG_terms_full_subset.pdf", 
       height = 8, width = 10)


full.enrich %>% 
  filter(!description %in% cat_removals) %>%
  filter((description %in% selected_kegg_enrich)) %>%
  mutate(direction=factor(direction, levels = c('DOWN', 'UP'))) %>%
  mutate(FDR = cut(fdr, 
                   labels = c('0.001', '0.01', '0.05', 'ns'),
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf)),
         .before = term) %>% write_csv("summary/FULL_version/enrichment/KEGG_terms_full_subset.csv")

# Volcano plots FULL ---------------------------------------------------



stats.full %>% 
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

ggsave(here('summary', 'FULL_version', 'volcano_full.pdf'), height = 7, width = 10)
