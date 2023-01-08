#this preps all the gwas summary statistics pulled from the gwas catalog for the follow-up LDpred analyses

library(tidyverse)
library(vroom)

pathHapMap<-"/path/to/map_hm3_ldpred2.rds" #pre-compiled hapmap3 reference
pathGWAS<-"/path/to/gwas-catalog/gwas/summary_stats/"
pathOutput<-"/path/to/gwas_LDpred_ready/"

info <- readRDS(pathHapMap) %>%
  mutate(MAF=ifelse(af_UKBB<0.5,af_UKBB,1-af_UKBB))

snp_list<-c()

n_eff_fun<-function(n_case, n_control){
  return(4 / (1 / n_control + 1 / n_case))
}

#grch37 chr6:28,477,797-33,448,354
#grch38 chr6:28,510,120-33,480,577

#asthma 37
asthma<-vroom(paste0(pathGWAS, "asthma_eur_32296059.txt.gz")) %>%
  rename(chr=CHR) %>%
  rename(pos=BP) %>%
  rename(rsid=SNP) %>%
  rename(a1=EA) %>%
  rename(a0=NEA) %>%
  mutate(n_eff=n_eff_fun(n_case=64538, n_control=329321)) %>%
  mutate(beta_se=(log(OR_95U)-log(OR_95L))/3.92) %>%
  rename(p=P) %>%
  rename(OR=OR) %>%
  rename(INFO=INFO) %>%
  mutate(MAF=ifelse(EAF<0.5,EAF,1-EAF)) %>%
  dplyr::select(c(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR, INFO, MAF)) %>%
  filter(rsid %in% info$rsid)

snp_list<-unique(c(snp_list, asthma$rsid))

vroom_write(asthma, paste0(pathOutput, "asthma_with_hla.tsv.gz"))

asthma %>%
  filter(!(chr=6 & pos > 27977797 & pos<33948354)) %>%
  vroom_write(., paste0(pathOutput, "asthma_without_hla.tsv.gz"))

#celiac 37
celiac<-vroom(paste0(pathGWAS, "celiac_eur_20190752.txt.gz")) %>%
  rename(chr=chromosome) %>%
  rename(pos=base_pair_location) %>%
  rename(rsid=variant_id) %>%
  rename(a1=effect_allele) %>%
  rename(a0=other_allele) %>%
  mutate(n_eff=n_eff_fun(n_case=4533, n_control=10750)) %>%
  rename(beta_se=standard_error) %>%
  rename(p=p_value) %>%
  rename(OR=odds_ratio) %>%
  mutate(INFO=rep(1,nrow(.))) %>%
  mutate(MAF=ifelse(effect_allele_frequency<0.5,effect_allele_frequency,1-effect_allele_frequency)) %>%
  dplyr::select(c(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR, INFO, MAF)) %>%
  filter(rsid %in% info$rsid)

snp_list<-unique(c(snp_list, celiac$rsid))

vroom_write(celiac, paste0(pathOutput, "celiac_with_hla.tsv.gz"))

celiac %>%
  filter(!(chr=6 & pos > 27977797 & pos<33948354)) %>%
  vroom_write(., paste0(pathOutput, "celiac_without_hla.tsv.gz"))
  
#dm1 38
dm1<-vroom(paste0(pathGWAS, "dm1_eur_34012112.tsv.gz")) %>%
  rename(chr=chromosome) %>%
  rename(pos=base_pair_location) %>%
  rename(rsid=variant_id) %>%
  rename(a1=effect_allele) %>%
  rename(a0=other_allele) %>%
  mutate(n_eff=n_eff_fun(n_case=18942, n_control=501638)) %>%
  rename(beta_se=standard_error) %>%
  rename(p=p_value) %>%
  mutate(OR=exp(beta)) %>%
  mutate(INFO=rep(1,nrow(.))) %>%
  mutate(MAF=ifelse(effect_allele_frequency<0.5,effect_allele_frequency,1-effect_allele_frequency)) %>%
  dplyr::select(c(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR, INFO, MAF)) %>%
  filter(rsid %in% info$rsid)

snp_list<-unique(c(snp_list, dm1$rsid))

vroom_write(dm1, paste0(pathOutput, "dm1_with_hla.tsv.gz"))

dm1 %>%
  filter(!(chr=6 & pos > 28010120 & pos<33980577)) %>%
  vroom_write(., paste0(pathOutput, "dm1_without_hla.tsv.gz"))
   
#ms 37
ms<-vroom(paste0(pathGWAS, "ms_eur_24076602.tsv.gz")) %>%
  rename(chr=chrom) %>%
  rename(pos=pos) %>%
  rename(rsid=rsid) %>%
  rename(a1=effect_allele) %>%
  rename(a0=other_allele) %>%
  mutate(n_eff=n_eff_fun(n_case=14498, n_control=24091)) %>%
  rename(beta_se=se) %>%
  rename(p=p) %>%
  rename(OR=OR) %>%
  mutate(INFO=rep(1,nrow(.))) %>%
  merge(.,info[,c("rsid","MAF")]) %>%
  dplyr::select(c(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR, INFO, MAF)) %>%
  filter(rsid %in% info$rsid)

snp_list<-unique(c(snp_list, ms$rsid))

vroom_write(ms, paste0(pathOutput, "ms_with_hla.tsv.gz"))

ms %>%
  filter(!(chr=6 & pos > 27977797 & pos<33948354)) %>%
  vroom_write(., paste0(pathOutput, "ms_without_hla.tsv.gz"))
     
#psoriasis 37
psoriasis<-vroom(paste0(pathGWAS, "psoriasis_eur_23143594.tsv.gz")) %>%
  rename(chr=chrom) %>%
  rename(pos=pos) %>%
  rename(rsid=rsid) %>%
  rename(a1=effect_allele) %>%
  rename(a0=other_allele) %>%
  mutate(n_eff=n_eff_fun(n_case=10588, n_control=22806)) %>%
  rename(beta_se=se) %>%
  rename(p=p) %>%
  rename(OR=OR) %>%
  mutate(INFO=rep(1,nrow(.))) %>%
  merge(.,info[,c("rsid","MAF")]) %>%
  dplyr::select(c(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR, INFO, MAF)) %>%
  filter(rsid %in% info$rsid)

snp_list<-unique(c(snp_list, psoriasis$rsid))

vroom_write(psoriasis, paste0(pathOutput, "psoriasis_with_hla.tsv.gz"))

psoriasis %>%
  filter(!(chr=6 & pos > 27977797 & pos<33948354)) %>%
  vroom_write(., paste0(pathOutput, "psoriasis_without_hla.tsv.gz"))
    
  
#rheumatoid arthritis 37
ra<-vroom(paste0(pathGWAS, "RA_eur_eas_24390342.txt.gz")) %>%
  filter(SNPID %in% info$rsid) %>%
  separate(col=`OR_95%CIup-OR_95%CIlow`, into=c("OR_high", "OR_low"), sep="-") %>%
  mutate(OR_high=as.numeric(OR_high)) %>%
  mutate(OR_low=as.numeric(OR_low)) %>%
  rename(chr=Chr) %>%
  rename(pos=`Position(hg19)`) %>%
  rename(rsid=SNPID) %>%
  rename(a1=A1) %>%
  rename(a0=A2) %>%
  mutate(n_eff=n_eff_fun(n_case=29880, n_control=73758)) %>%
  mutate(beta_se=(log(OR_high)-log(OR_low))/3.92) %>%
  rename(p=`P-val`) %>%
  rename(OR=`OR(A1)`) %>%
  mutate(INFO=rep(1,nrow(.))) %>%
  merge(.,info[,c("rsid","MAF")]) %>%
  dplyr::select(c(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR, INFO, MAF))

snp_list<-unique(c(snp_list, ra$rsid))

vroom_write(ra, paste0(pathOutput, "ra_with_hla.tsv.gz"))

ra %>%
  filter(!(chr=6 & pos > 27977797 & pos<33948354)) %>%
  vroom_write(., paste0(pathOutput, "ra_without_hla.tsv.gz"))
    
#ulcerative colitis 37
uc<-vroom(paste0(pathGWAS, "UC_eur_26192919.txt.gz")) %>%
  rename(chr=Chr) %>%
  rename(pos=Pos) %>%
  rename(rsid=SNP) %>%
  rename(a1=A1_effect) %>%
  rename(a0=A2_other) %>%
  mutate(n_eff=n_eff_fun(n_case=6968, n_control=20464)) %>%
  rename(beta_se=se_EUR) %>%
  rename(p=P_EUR) %>%
  mutate(OR=exp(beta_EUR)) %>%
  mutate(INFO=rep(1,nrow(.))) %>%
  merge(.,info[,c("rsid","MAF")]) %>%
  dplyr::select(c(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR, INFO, MAF)) %>%
  filter(rsid %in% info$rsid)

snp_list<-unique(c(snp_list, uc$rsid))

vroom_write(uc, paste0(pathOutput, "uc_with_hla.tsv.gz"))

uc %>%
  filter(!(chr=6 & pos > 27977797 & pos<33948354)) %>%
  vroom_write(., paste0(pathOutput, "uc_without_hla.tsv.gz"))
  
#save list of snps
write_lines(snp_list, paste0(pathOutput, "snp_list.txt.gz"))
