# BF528-Project-5
This is the last project for BF528 as an individual final 

# Project Description
Colorectal cancer (CRC) is the third most common type of cancer and the fourth most common cause of death in the world. Pathological staining is the most widely used method to determine the presence of CC, however, it does not accurately predict recurrence of CRC. About 20% of patients who are diagnosed with stage II or III CRC develop recurrence. Researchers have looked into Gene Expression Profiles (GEPs) through the use of microarrays. The Marisa et al study established a classification of the CC subtypes based on their molecular features by exploiting “genome-wide mRNA expression analysis” through the use of microarrays. Initially only three subtypes of CC were identified, but through the Marisa et al study, six subtypes were classified, more accurately reflecting the molecular heterogeneity of CC. 

This analysis focuses only on reproducing the programmer and biologist results from the comparison of the C3 and C4 tumor subtypes from Marisa et al.

# Referenece:
Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

# Contributors
Italo Duran: Programmer & Biologist (duran01@bu.edu)

# Programmer
The programmer.R script normalizes the data by using the RMA method to compute standard quality control metrics on the normalized data and visualize the distribution of samples using Principal Component Analysis (PCA).

# Biologist
The biologist seeks to understand the biological significance of the different gene expression profiles for each tumor subtype using gene set enrichment analysis. And reproduces an analysis using KEGG gene sets and the differential expression results from 5.6 or the given sample data. 


