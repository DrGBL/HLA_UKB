# This provides code to check the adjusted concordance

library(tidyverse)
library(vroom)

#again go to the directory where the folders "calls" and "QC" are located
setwd("/path/to/folder/")

path_anc<-"/path/to/ancestry/files/"
path_munged_imputed<-"/path/to/hla_imputation_df_munged.tsv.gz"


#ancestries
afr<-scan(paste0(path_anc, "ukb.afrIDsPCA.txt"))
amr<-scan(paste0(path_anc, "ukb.amrIDsPCA.txt"))
eas<-scan(paste0(path_anc, "ukb.easIDsPCA.txt"))
eur<-scan(paste0(path_anc, "ukb.eurIDsPCA.txt"))
sas<-scan(paste0(path_anc, "ukb.sasIDsPCA.txt"))

#imputed data frame
imputed<-vroom(path_munged_imputed, 
               col_types = cols(.default = "c"))
               
#sequenced data frames
sequenced<-c()              
for(i in 10:60){
  sequenced<-vroom(paste0("calls/hla_df_batch_qced_", i, ".tsv.gz")) %>%
  filter(ID %in% imputed$ID) %>%
  dplyr::select(colnames(imputed)) %>%
  bind_rows(sequenced,.)
}


sequenced<-sequenced %>%
  mutate(across(-c(ID),~str_extract(.x,"[A-Z0-9]*\\*[0-9]*:[0-9]*")))

imputed<-imputed %>%
  filter(ID %in% sequenced$ID)

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
  max_concordant_alleles[1,g]<-sum(num_in_hla_imp[[g]][1,2:ncol(num_in_hla_imp[[g]])])
  max_concordant_alleles[2,g]<-sum(num_in_hla_imp[[g]][2,2:ncol(num_in_hla_imp[[g]])])
  max_concordant_alleles[3,g]<-sum(num_in_hla_imp[[g]][3,2:ncol(num_in_hla_imp[[g]])])
  max_concordant_alleles[4,g]<-sum(num_in_hla_imp[[g]][4,2:ncol(num_in_hla_imp[[g]])])
  max_concordant_alleles[5,g]<-sum(num_in_hla_imp[[g]][5,2:ncol(num_in_hla_imp[[g]])])
  max_concordant_alleles[6,g]<-sum(num_in_hla_imp[[g]][6,2:ncol(num_in_hla_imp[[g]])])
}

write_tsv(max_concordant_alleles, "/path/to/max_concordant_alleles.tsv")
