## BF528 - FINAL PROJECT - INDIVIDUAL PROJECT
## Author: Italo Duran
## email: duran01@bu.edu

library(BiocManager)
library(hgu133plus2.db)
library(GSEABase)
library(AnnotationDbi)
library(affy)
library(tidyverse)
#library(dplyr)
setwd("/projectnb/bf528/students/duran01/individual_project/biologist")

### 6.1: mapping the probeset IDs to gene symbols by specifying the appropriate key and column arguments:
de_data = read.csv("/project/bf528/project_1/data/differential_expression_results.csv",row.names=1,header=TRUE)
de_data = de_data %>% arrange(desc(t))
view(de_data)
de_matches = AnnotationDbi::select(hgu133plus2.db,keys=as.character(row.names(de_data)),columns = ("SYMBOL"))
dup_de = de_matches[!duplicated(de_matches[1]),] #taking out duplicates from data
de_data = cbind(dup_de, de_data)
de_data = de_data %>%
  group_by(SYMBOL)%>%
  filter(padj==min(padj))%>%
  ungroup(SYMBOL)

### 6.3 & 6.4: selecting the top 1000 up- and down-regulated genes:
#load the data from the 3 sets:
hall_set = getGmt('h.all.v7.5.1.symbols.gmt.txt')
go_set = getGmt('c5.go.v7.5.1.symbols.gmt.txt')
kegg_set = getGmt('c2.cp.kegg.v7.5.1.symbols.gmt.txt')
de_data = de_data[!is.na(de_data $ SYMBOL),] #cleaning the data, deleting the empty values
up_regulated_1K = head(de_data, 1000) #getting the top 1000 for up-regulated
dwn_regulated_1K = tail(de_data, 1000) #getting the top 1000 for down-regulated
up_regulated_top10 = head(up_regulated_1K, 10) #getting the top 10 for up-regulated
down_regulated_top10 = tail(dwn_regulated_1K, 10) #getting the top 10 for down-regulated
write.csv(up_regulated_top10, "Top10_up_reg_genes.csv") #export the top 10 up regulated
write.csv(down_regulated_top10, "Top10_down_reg_genes.csv") #export the top 10 down regulated
#cleaning, removing the genes that were not used in the data set:
notexp_up = subset(de_data, ! de_data $ SYMBOL %in% up_regulated_1K $ SYMBOL)
notexp_dwn = subset(de_data, ! de_data $ SYMBOL %in% dwn_regulated_1K $ SYMBOL)

### 6.5: Using the fisher.test function to compute hypergeometric statistics and p-values:
#defining fishers test on the data set:
fishertest = function(gl, gs, nde){ 
  diffexp_in_gs = length(intersect(gl,gs))    
  diffexp_not_gs = length(gl) - diffexp_in_gs   
  notexp_in_gs = length(intersect(nde,gs))      
  notexp_not_gs = length(nde) - notexp_in_gs      
  return(c(diffexp_in_gs,diffexp_not_gs,notexp_in_gs,notexp_not_gs)) }  

## fisher test for hallmark geneset:
fisher_hallmark_resl = data.frame(setname=character(),pvalue=numeric(),estimate=numeric(),exp=character(),stringsAsFactors=FALSE)
for (i in 1:length(hall_set)){
  set_genid = set_genids(hall_set[i])
  up_reg_fishc_result = fishertest(up_regulated_1K $ SYMBOL,set_genid[[names(set_genid)]], notexp_up$SYMBOL)
  dwn_reg_fishc_result = fishertest(dwn_regulated_1K $ SYMBOL,set_genid[[names(set_genid)]], notexp_dwn$SYMBOL)
  up_Reg_df = fisher.test(matrix(up_reg_fishc_result,nrow=2))
  down_Reg_df = fisher.test(matrix(dwn_reg_fishc_result, nrow=2))
  fisher_hallmark_resl[nrow(fisher_hallmark_resl) +1, ] = c(names(set_genid),up_Reg_df $ p.value,up_Reg_df $ estimate, 'UP')
  fisher_hallmark_resl[nrow(fisher_hallmark_resl) +1, ] = c(names(set_genid),down_Reg_df $ p.value,down_Reg_df $ estimate, 'DOWN') }
fisher_hallmark_resl = fisher_hallmark_resl %>% 
  mutate(pvalue = as.numeric(pvalue),estimate = as.numeric(estimate))
View(fisher_hallmark_resl)

## fisher test for kegg gene-set:
fisch_kegg_resul = data.frame(setname=character(),pvalue=numeric(),estimate=numeric(),exp=character(),stringsAsFactors=FALSE)
for (i in 1:length(kegg_set)){
  set_genid = set_genids(kegg_set[i])
  up_reg_fishc_result = fishertest(up_regulated_1K $ SYMBOL,set_genid[[names(set_genid)]],notexp_up$SYMBOL)
  dwn_reg_fishc_result = fishertest(dwn_regulated_1K $ SYMBOL,set_genid[[names(set_genid)]], notexp_dwn$SYMBOL)
  up_Reg_df = fisher.test(matrix(up_reg_fishc_result,nrow=2))
  down_Reg_df = fisher.test(matrix(dwn_reg_fishc_result, nrow=2))
  fisch_kegg_resul[nrow(fisch_kegg_resul) +1, ] = c(names(set_genid),up_Reg_df $ p.value,up_Reg_df $ estimate, 'UP')
  fisch_kegg_resul[nrow(fisch_kegg_resul) +1, ] = c(names(set_genid),down_Reg_df $ p.value,down_Reg_df $ estimate, 'DOWN') }
fisch_kegg_resul = fisch_kegg_resul %>% mutate(pvalue = as.numeric(pvalue),estimate = as.numeric(estimate))
View(fisch_kegg_resul)

#fisher test for GO gene-set:
fisch_go_resul = data.frame(setname=character(),pvalue=numeric(),estimate=numeric(),exp=character(),stringsAsFactors=FALSE)
for (i in 1:length(go_set)){
  set_genid = set_genids(go_set[i])
  up_reg_fishc_result = fishertest(up_regulated_1K $ SYMBOL,set_genid[[names(set_genid)]], notexp_up$SYMBOL)
  dwn_reg_fishc_result = fishertest(dwn_regulated_1K $ SYMBOL,set_genid[[names(set_genid)]], notexp_dwn$SYMBOL)
  up_Reg_df = fisher.test(matrix(up_reg_fishc_result,nrow=2))
  down_Reg_df = fisher.test(matrix(dwn_reg_fishc_result, nrow=2))
  fisch_go_resul[nrow(fisch_go_resul) +1, ] = c(names(set_genid),up_Reg_df $ p.value,up_Reg_df $ estimate,'UP')
  fisch_go_resul[nrow(fisch_go_resul) +1, ] = c(names(set_genid),down_Reg_df $ p.value,down_Reg_df $ estimate,'DOWN') }
fisch_go_resul = fisch_go_resul %>% mutate(pvalue = as.numeric(pvalue),estimate = as.numeric(estimate))
View(fisch_go_resul)

#getting the total for the p-value adjusted for GO:
fisch_go_resul$BH = p.adjust(fisch_go_resul$pvalue, method = "BH",n = length(fisch_go_resul $ pvalue))
write.csv(fisch_go_resul, "go_total.csv")     
#getting the total for the p-value adjusted for Kegg:
fisch_kegg_resul$BH = p.adjust(fisch_kegg_resul$pvalue, method = "BH",n = length(fisch_kegg_resul $ pvalue))
write.csv(fisch_kegg_resul, "kegg_total.csv")
#getting the total for the p-value adjusted for Hallmark:
fisher_hallmark_resl$BH = p.adjust(fisher_hallmark_resl$pvalue, method = "BH",n = length(fisher_hallmark_resl $ pvalue))
write.csv(fisher_hallmark_resl, "hallmark_total.csv")

# Hallmark, Statistically significant:
hall_sig_enrch = fisher_hallmark_resl[fisher_hallmark_resl$pvalue < 0.05,]
hall_sig_enrchTL = length(hall_sig_enrch $ gene_set)
# Kegg, Statistically significant:
kegg_sig_enrch = fisch_kegg_resul[fisch_kegg_resul$pvalue < 0.05,]
kegg_sig_enrchTL = length(kegg_sig_enrch$gene_set)
# GO, Statistically significant:
go_sig_enrch = fisch_go_resul[fisch_go_resul$pvalue < 0.05,]
go_sig_enrchTL = length(go_sig_enrch$gene_set)
# Top 3 sets for each one:
hall_top3_sets = slice_min(fisher_hallmark_resl, order_by = pvalue, n=3)
kegg_top3_sets = slice_min(fisch_kegg_resul, order_by = pvalue, n=3)
go_top3_sets = slice_min(fisch_go_resul, order_by = pvalue, n=3)
#all top gene sets for all 3: 
All_top3_sets = rbind(kegg_top3_sets, go_top3_sets, hall_top3_sets)
print(All_top3_sets)
write.csv(All_top3_sets, file="top3_enriched_gene_sets.csv")