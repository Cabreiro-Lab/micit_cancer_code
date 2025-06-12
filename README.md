# Chemotherapy Modulation by a Cancer-Associated Microbiota Metabolite

In this repository, we provide the code and data used in key figures from the manuscript "Chemotherapy Modulation by a Cancer-Associated Microbiota Metabolite" by **Daniel Martinez-Martinez, Tanara V. Peres, Kristin Gehling _et al._**

Leading contact: Filipe Cabreiro (f.cabreiro@uni-koeln.de)

## Abstract

Understanding how the microbiota produces regulatory metabolites is of significance for cancer and cancer therapy. Using a host-microbe-drug-nutrient 4-way screening approach, we evaluated the role of nutrition at the molecular level in the context of 5-fluorouracil toxicity. Notably, our screens identified the metabolite 2-methylisocitrate which was found to be produced and enriched in human tumor-associated microbiota. 2-methylisocitrate exhibits anti-proliferative properties across genetically- and tissue-diverse cancer cell lines, 3D spheroids, and an in vivo Drosophila gut tumor model, where it reduced tumor dissemination and increased survival. Chemical landscape interaction screens identified drug-metabolite signatures and highlighted the synergy between 5-fluorouracil and 2-methylisocitrate. Multi-omic analyses revealed that 2-methylisocitrate acts via multiple cellular pathways linking metabolism and DNA damage to regulate chemotherapy. Finally, we converted 2-methylisocitrate into its trimethyl ester, thereby enhancing its potency. This work highlights the great impact of microbiome-derived metabolites on tumor proliferation, and their potential as promising co-adjuvants for cancer treatment.

## Data

The data used in this study is available in the `supplementary_tables` folder. This data is mostly present in the Supplemental Tables but also here for convenience. The Zenodo repository from this article may contain some necessary data to run these scripts. Additional tables are provided here to perform some of the analyses.

The data is organized as follows:

**KEGGEcoCycClass.csv**: List of KEGG and EcoCyc pathways and their corresponding classes for the nutrients in the Biolog plates.

**ecocyc_paths.csv**: List of pathways from the EcoCyc database.

**core_tree_mod.treefile**: Phylogenetic tree from Figure 1G

**core_tree_metadata.csv**: Metadata needed to plot the phylogenetic tree


## Code

The code used in this study is available in the `code` folder. The code is organized as follows:

**4_way_plots.R**: R script to generate the 4-way plots in Figure 2 and S2.

**gene_nutrient_screen.R**: R script to generate the gene-nutrient screen in Figure 3 and S3. Double mutant plot in Figure S4. 

**micit_producers.R**: R script to generate the micit production in the different cohorts in Figure 5 and S5. 

**cancer_cell_viability.R**: R script to generate the cancer cell viability in Figure 5 and cell growth in Figure S5.

**phylo_tree.R**: R script to generate the phylogenetic tree with worm and bacterial phenotypes from Figure 1. 

**Biolog_cancer_cells.R**: R script to generate the original analysis for the Biolog data and the chemical space multivariate analysis. *Legacy code*.

**HCT116_wt_p53.R**: R script to generate the HCT116 wild-type and p53 knockout cell line analysis. *Legacy code*.

**isomer_analysis.R**: R script to generate the isomer analysis in Figure 7. *Legacy code*.

**proteomics_5FU_micit.R**: R script to generate the proteomics analysis in Figure 7. *Legacy code*.

**RNA_seq_CRC.R**: R script to generate the RNA-seq analysis in Figure 6. *Legacy code*.

NOTE: Figures may not look as in the final version as those were further edited in Adobe Illustrator for publication purposes. The code is provided as is, and may require additional dependencies to run.