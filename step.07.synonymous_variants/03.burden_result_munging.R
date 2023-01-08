library(tidyverse)
library(vroom)

setwd("path/to/burden/test/results/")

pheno<-c("asthma", "auto_immune_thyroid_disorders", "dm1", "rheumatoid_arthritis", "ulcerative_colitis",
         "ankylosing_spondylitis", "coeliac", "crohn", "ms_demyelinating", "pmr_gcc", "psoriasis")

anc<-c("afr", "amr", "eas", "eur", "sas")

for(p in pheno){
  for(a in anc){
    if(file.exists(paste0("hla_ukb_wes_burden_", a, "_", p, ".regenie.gz"))){
      tmp<-vroom(paste0("hla_ukb_wes_burden_", a, "_", p, ".regenie.gz"), skip=1) %>%
        filter(!is.na(Effect)) %>%
        separate(Info, into=c("beta", "se", "mac"), sep=";") %>%
        mutate(beta=as.numeric(substr(beta, nchar("REGENIE_BETA=")+1, 100))) %>%
        mutate(se=as.numeric(substr(se, nchar("REGENIE_SE=")+1, 100)))
      
      write_tsv(tmp,paste0("hla_ukb_wes_burden_", a, "_", p, "_adj.tsv.gz"))
    }
  }
}
