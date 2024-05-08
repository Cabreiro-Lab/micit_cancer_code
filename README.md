# Chemotherapy Modulation by a Cancer-Associated Microbiota Metabolite

In this repository, we provide the code and data used in the manuscript "Chemotherapy Modulation by a Cancer-Associated Microbiota Metabolite" by **Daniel Martinez-Martinez, Tanara V. Peres _et al._**

## Abstract

Coming soon...

## Data

The data used in this study is available in the `data` folder. This data is mostly present in the Supplemental Tables but also here for convenience. You may want to also download the data uploaded to Zenodo, but it is not necessary for the code to run. Additional tables are provided here to perform some of the analyses.
The data is organized as follows:

**KEGGEcoCycClass.csv**: List of KEGG and EcoCyc pathways and their corresponding classes for the nutrients in the Biolog plates.

**ecocyc_paths.csv**: List of pathways from the EcoCyc database.

## Code

The code used in this study is available in the `code` folder. The code is organized as follows:

**4_way_plots.R**: R script to generate the 4-way plots in Figure 2 and S2.

**gene_nutrient_screen.R**: R script to generate the gene-nutrient screen in Figure 3 and S3. Double mutant plot in Figure S4. 

**micit_producers.R**: R script to generate the micit production in the different cohorts in Figure 5 and S5. 

**cancer_cell_viability.R**: R script to generate the cancer cell viability in Figure 5 and cell growth in Figure S5.

**Biolog_cancer_cells.R**: R script to generate the original analysis for the Biolog data and the chemical space multivariate analysis. *Legacy code*.

**HCT116_wt_p53.R**: R script to generate the HCT116 wild-type and p53 knockout cell line analysis. *Legacy code*.

**isomer_analysis.R**: R script to generate the isomer analysis in Figure 7. *Legacy code*.

**proteomics_5FU_micit.R**: R script to generate the proteomics analysis in Figure 7. *Legacy code*.

**RNA_seq_CRC.R**: R script to generate the RNA-seq analysis in Figure 6. *Legacy code*.

NOTE: Figures may not look as in the final version as those were further edited in Adobe Illustrator for publication purposes. The code is provided as is, and may require additional dependencies to run.