library(tidyverse)
library(vroom)
library(bigsnpr)
library(Matrix)

pathIn<-"/path/to/gwas_ld_pred_ready/"
pathOut<-"/path/to/ld_pred_output/"

h2<-read_tsv(paste0(pathIn, "ldsc_h2.tsv"))


mapping<-readRDS(paste0(pathIn, "ld_pred_ld_eur/map.rds"))

tmp <- tempfile(tmpdir = paste0(pathIn,"tmp-data"))
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

file.remove(paste0(tmp, ".sbk"))

snps<-c()

for(ana in c("_with_hla", "_without_hla")){
  for(x in c("dm1", "psoriasis", "asthma", "celiac", "ms",  "ra", "uc")){   #"dm1", "psoriasis", later
    gwas<-vroom(paste0(pathIn, x, ana, ".tsv.gz")) %>%
      mutate(beta=log(OR)) %>%
      filter(beta_se>0) %>%
      group_by(rsid) %>%
      arrange(desc(MAF)) %>%
      filter(row_number()==1) %>%
      ungroup()
    
    for(chrom in 1:22){
      corr0<-readRDS(paste0(pathIn, "ld_pred_ld_eur/LD_with_blocks_chr", chrom, ".rds"))
      
      map_tmp<-mapping %>% filter(chr==chrom)
      
      included_snps<-which(map_tmp$rsid %in% gwas$rsid)
      
      corr0<-corr0[included_snps,included_snps]
      
      if (chrom == 1) {
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        corr$add_columns(corr0, nrow(corr))
      }
    }
    
    beta_inf <- snp_ldpred2_inf(corr=corr, df_beta=gwas[,c("beta", "beta_se", "n_eff")], h2 = h2$h2[which(h2$phenotype==x)])
    
    ld_res<-data.frame(rsid=gwas$rsid, effect_allele=gwas$a1, beta_inf=beta_inf) %>%
      filter(is.finite(beta_inf))
    
    vroom_write(ld_res, paste0(pathOut,"ld_pred_",x,ana,".tsv.gz"))
    
    snps<-ld_res %>% dplyr::select(rsid) %>% bind_rows(snps)
    
    file.remove(paste0(tmp, ".sbk"))
  }
}

snps<-snps %>% distinct()
vroom_write(snps, paste0(pathOut, "snps_list.tsv"), col_names=FALSE)
