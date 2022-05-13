# BF528-Project-5-Individual
This is the last project for BF528 as an individual final 

# Project Description
Colorectal cancer (CRC) is the third most common type of cancer and the fourth most common cause of death in the world. Pathological staining is the most widely used method to determine the presence of CC, however, it does not accurately predict recurrence of CRC. Researchers have looked into Gene Expression Profiles (GEPs) through the use of microarrays. The Marisa et al study established a classification of the CC subtypes based on their molecular features by exploiting “genome-wide mRNA expression analysis” through the use of microarrays. Initially only three subtypes of CC were identified, but through the Marisa et al study, six subtypes were classified, more accurately reflecting the molecular heterogeneity of CC. 

This analysis focuses only on reproducing the programmer and biologist results from, the comparison of the C3 and C4 tumor subtypes. A short analysis of microarray dataset in Marisa, et. al. https://pubmed.ncbi.nlm.nih.gov/23700391/

# Referenece
Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

# Contributors
Italo Duran: Programmer & Biologist (duran01@bu.edu)

# Programmer
This folder contains the programmers deliverables:
The programmer.R script normalizes the data by using the RMA method to compute standard quality control metrics on the normalized data and visualize the distribution of samples using Principal Component Analysis (PCA).
Along with the RLE, NUSE Histograms and the PCA scatter plot. 

# Biologist
This folder contains the biologist deliverables:
The biologist seeks to understand the biological significance of the different gene expression profiles for each tumor subtype using gene set enrichment analysis. And reproduces an analysis using KEGG gene sets and the differential expression results from 5.6 or the given sample data.
It also includes:
The tables for the top 10 up- and down-regulated probesets with gene symbol, t-statistic, nominal p-value, and adjusted p-value columns.
A description table of the gene set databases used, specifying the number of gene sets considered in each.
A table containing the top 3 enriched gene sets for each geneset type.



