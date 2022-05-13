## BF528 - FINAL PROJECT - INDIVIDUAL PROJECT
## Author: Italo Duran
## email: duran01@bu.edu

### installing the required mangers for the packages we need:
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(version = "3.14")
#BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))

## packages needed:
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)

#setting a working directory:
setwd("/projectnb/bf528/students/duran01/individual_project/programmer")

## 3. 
##Read in CEL files:
CEL_DATA_FILES = ReadAffy(celfile.path="/projectnb/bf528/students/duran01/individual_project/programmer/CEL_files")
# Using RMA function to normalize CEL_files:
RMA_DATA = rma(CEL_DATA_FILES, background=FALSE)

## 4. Quality Assessment:
normlz_CELdata_PLM = fitPLM(CEL_DATA_FILES,normalize=TRUE,background=TRUE)
# Relative Log Expression (RLE):
RLE_CELdata = RLE(normlz_CELdata_PLM,type = 'stats')
# Normalized Unscaled Standard Error (NUSE):
NUSE_CELdata = NUSE(normlz_CELdata_PLM,type = 'stats')
# Plot RLE & NUSE histograms:
jpeg(file = 'RLE.jpeg')
hist(RLE_CELdata[1,], main = 'Median Relative Log Expression (RLE)', xlab ='Median RLE Score', c = '#EEB462')
dev.off()
jpeg(file = 'NUSE.jpeg')
hist(NUSE_CELdata[1,], main = 'Median Normalized Unscaled Standard Error (NUSE)', xlab ='Median NUSE Score', c = '#138086')
dev.off()

## 5: Correcting the metadata for batch effects:
# Loading metadata:
proj_metadata = read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
norm_data_batch = as.numeric(factor(proj_metadata$normalizationcombatbatch)) 
norm_mod_data =  as.numeric(factor(proj_metadata$normalizationcombatmod))
# combat normalized batch corrected data:
ComBat(RMA_DATA,batch = norm_data_batch,mod = norm_mod_data)
data_fram = exprs(RMA_DATA)
write.csv(data_fram,"/projectnb/bf528/students/duran01/individual_project/programmer/RMA_normalized.csv", row.names = TRUE)

## 6: Principal component analysis:
## Transpose for scaling the data frame:
combt_trns_scalDT = scale(t(data_fram))
scalDT = t(combt_trns_scalDT)
princomps = prcomp(scalDT,center = FALSE,scale. =FALSE)
pca1 = princomps$rotation[,1]
pca2 = princomps$rotation[,2]
jpeg(file ='pca.jpeg')
plot(pca1,pca2, xlab ='PC1 14.5%', ylab ='PC2 9.54%', main ='Principal component analysis (PCA) Plot',pch=16, col='#138086')
dev.off()
## 7: xamine the percent variability explained by each principal component
princomps$importance
summary(princomps)