# this code adjusts the summary statistics to make them ready for meta-analysis

library(tidyverse)

setwd("/path/to/regenie_step_2_outputs/")

pheno<-c("asthma", "auto_immune_thyroid_disorders", "dm1", "lupus", "rheumatoid_arthritis", "ulcerative_colitis",
          "ankylosing_spondylitis", "coeliac", "crohn", "ms_demyelinating", "myasthenia_gravis", "pernicious_anemia",
          "pmr_gcc", "psoriasis", "sjogren")

anc<-c("afr", "amr", "eas", "eur", "sas")

digit<-c("four", "six")

for(p in pheno){
  for(a in anc){
    for(d in digit){
      if(file.exists(paste0("hla_ukb_wes_", d, "_digit_", a, "_", p, ".regenie.gz"))){
        tmp<-read_tsv(paste0("hla_ukb_wes_", d, "_digit_", a, "_", p, ".regenie.gz")) %>%
          filter(!is.na(Effect)) %>%
          separate(Info, into=c("beta", "se", "mac"), sep=";") %>%
          mutate(beta=as.numeric(substr(beta, nchar("REGENIE_BETA=")+1, 100))) %>%
          mutate(se=as.numeric(substr(se, nchar("REGENIE_SE=")+1, 100))) %>%
          dplyr::select(-c(Effect, LCI_Effect, UCI_Effect, mac)) %>%
          mutate(beta=ifelse(Ref=="A", -1*beta, beta)) %>%
          mutate(AAF=ifelse(Ref=="A", 1-AAF, AAF)) %>%
          mutate(Cases_Ref_tmp=ifelse(Ref=="A", Cases_Alt, Cases_Ref)) %>%
          mutate(Cases_Alt_tmp=ifelse(Ref=="A", Cases_Ref, Cases_Alt)) %>%
          mutate(Controls_Ref_tmp=ifelse(Ref=="A", Controls_Alt, Controls_Ref)) %>%
          mutate(Controls_Alt_tmp=ifelse(Ref=="A", Controls_Ref, Controls_Alt)) %>%
          mutate(Cases_Alt=Cases_Alt_tmp) %>%
          mutate(Cases_Ref=Cases_Ref_tmp) %>%
          mutate(Controls_Alt=Controls_Alt_tmp) %>%
          mutate(Controls_Ref=Controls_Ref_tmp) %>%
          mutate(Ref="C") %>%
          mutate(Alt="A") %>%
          dplyr::select(-c(Cases_Ref_tmp, Cases_Alt_tmp, Controls_Ref_tmp, Controls_Alt_tmp))
        
        write_tsv(tmp,paste0("hla_ukb_wes_", d, "_digit_", a, "_", p, "_adj.tsv.gz"))
      }
    }
  }
}
