# this computes the heritability for each trait. It will be used by LDpred in the next steps.

library(tidyverse)
library(vroom)
library(bigsnpr)

pathOut<-"/path/to/gwas_ld_pred_ready/"
pathHapMap<-"/path/to/map_hm3_ldpred2.rds"

#options(bigstatsr.check.parallel.blas = FALSE)

info <- readRDS(pathHapMap) %>%
  dplyr::select(c(rsid,ld)) %>%
  rename(SNP=rsid) %>%
  rename(L2=ld)

ldsc_res<-data.frame(phenotype=c("asthma", "celiac", "dm1", "ms", "psoriasis", "ra", "uc"),
                      h2=NA)

for(x in c("asthma", "celiac", "dm1", "ms", "psoriasis", "ra", "uc")){
  gwas<-vroom(paste0(pathOut, x, "_with_hla.tsv.gz")) %>%
    filter(!is.na(beta_se) & !is.na(OR) & beta_se>0)

  
  ld_scores<-info %>%
    filter(SNP %in% gwas$rsid)
  
  gwas<-ld_scores %>%
    rename(rsid=SNP) %>%
    merge(gwas,.)
  
  ldsc<-snp_ldsc(ld_score=gwas$L2,
                  ld_size=nrow(gwas),
                  chi2 = (log(gwas$OR) / gwas$beta_se)^2,
                  sample_size = gwas$n_eff,
                  blocks = NULL)
    
  ldsc_res$h2[which(ldsc_res$phenotype==x)]<-ldsc[["h2"]]
}

write_tsv(ldsc_res, paste0(pathOut, "ldsc_h2.tsv"))
