library(tidyverse)
library(vroom)
library(plotly)

#again go to the directory where the folders "comparisons_with_imputed", "calls" and "QC" are located
setwd("/your/local/folder")

#path to ancestry files
path_anc<-"/your/ancestry/files"

#ancestries
afr<-scan(paste0(path_anc, "ukb.afrIDsPCA.txt"))
amr<-scan(paste0(path_anc, "ukb.amrIDsPCA.txt"))
eas<-scan(paste0(path_anc, "ukb.easIDsPCA.txt"))
eur<-scan(paste0(path_anc, "ukb.eurIDsPCA.txt"))
sas<-scan(paste0(path_anc, "ukb.sasIDsPCA.txt"))

comp<-c()
for(i in 10:60){
  comp<-vroom(paste0("comparisons_with_imputed/comparisons_imputed_called_folder_", i, ".tsv.gz")) %>%
    bind_rows(comp,.)
  if(i == 10){
    num_in_hla_imp<-readRDS(paste0("comparisons_with_imputed/num_alleles_in_imp_ref_", i, ".RDS"))
  } else {
    tmp<-readRDS(paste0("comparisons_with_imputed/num_alleles_in_imp_ref_", i, ".RDS"))
    for(g in names(num_in_hla_imp)){
      for(j in 2:ncol(num_in_hla_imp[[g]])){
        num_in_hla_imp[[g]][,j]<-num_in_hla_imp[[g]][,j]+tmp[[g]][,j]
      }
    }
  }
}

max_concordant_alleles<-data.frame(ancestry=c("all", "afr", "amr", "eas", "eur", "sas"))
for(g in names(num_in_hla_imp)){
  max_concordant_alleles<-max_concordant_alleles %>%
    mutate(!!sym(g):=NA)
  max_concordant_alleles[,g]<-num_in_hla_imp[[g]] %>% reframe(`Two matches`+`One match, other uncalled and unimputed`+
                                       `One match, other unimputed`+`One match, one mismatch imputation`+`One match, other uncalled but imputed`)
}

write_tsv(max_concordant_alleles, "comparisons_with_imputed/max_concordant_alleles.tsv")
