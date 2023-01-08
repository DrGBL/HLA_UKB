library(tidyverse)
library(vroom)
library(xgboost)

pathAnc<-"/path/to/ancestry/"
pathPred<-"/path/to/path_out_xgboost/"
pathPhenoCovar<-"/path/to/regenie_step_1_inputs/"
pathPCs<-"/path/to/pcs/"
pathLDPred<-"/path/to/ld_pred_output/"

#test ids
test_ids<-read_tsv(paste0(pathPred, "test_ids.tsv"))

#extract ancestry IDs
afr<-scan(paste0(pathAnc, "ukb.afrIDsPCA.txt"))
amr<-scan(paste0(pathAnc, "ukb.amrIDsPCA.txt"))
eas<-scan(paste0(pathAnc, "ukb.easIDsPCA.txt"))
eur<-scan(paste0(pathAnc, "ukb.eurIDsPCA.txt"))
sas<-scan(paste0(pathAnc, "ukb.sasIDsPCA.txt"))


#extract PCs
eigenvec <- vroom(paste0(pathPCs,'pca.eigenvec'), 
                  col_names = FALSE,
                  col_types = cols(.default = "c"))
colnames(eigenvec)<-c("sample", "IID", paste0("PC",c(1:20)))
eigenvec<-eigenvec %>%
  filter(sample>0) %>%
  dplyr::select(-c(IID, paste0("PC", c(11:20)))) %>%
  mutate_all(as.numeric)

#extract covar and pheno)
covar<-vroom(paste0(pathPhenoCovar,
                    "covar_file.tsv")) %>%
  dplyr::select(c(FID, age,sex)) %>%
  rename(sample=FID) %>%
  merge(., eigenvec) %>%
  mutate(sample=as.integer(sample)) %>%
  arrange(sample) %>%
  filter(sample %in% test_ids$sample)

rm(eigenvec)


#phenotype names
pheno_names<-data.frame(pheno=c("asthma",
                                "ms_demyelinating",
                                "psoriasis",
                                "dm1",
                                "rheumatoid_arthritis",
                                "ulcerative_colitis",
                                "coeliac"),
                        pheno_proper=c("Asthma",
                                       "MS-Demyelinating",
                                       "Psoriasis",
                                       "DM1",
                                       "Rheumatoid Arthritis",
                                       "Ulcerative Colitis",
                                       "Coeliac"),
                        pheno_name_gwas=c("asthma",
                                          "ms",
                                          "psoriasis",
                                          "dm1",
                                          "ra",
                                          "uc",
                                          "celiac"))

#loop to evaluate prediction
for(p in 1:nrow(pheno_names)){
  
  pheno<-vroom(paste0(pathPhenoCovar,
                      "pheno_file.tsv"))  %>%
    rename(sample=FID) %>%
    dplyr::select(-IID) %>%
    filter(sample %in% test_ids$sample) %>%
    dplyr::select(c(sample, pheno_names$pheno[p])) %>%
    mutate(sample=as.integer(sample)) %>%
    arrange(sample)
  
  #snp only
  print(paste("snp", pheno_names$pheno_proper[p]))
  pgs_snp<-vroom(paste0(pathLDPred, "ld_pred_scored_", pheno_names$pheno_name_gwas[p], "_with_hla.sscore")) %>%
    rename(sample=IID) %>%
    rename(ld_pred=SCORE1_SUM) %>%
    dplyr::select(c(sample, ld_pred)) %>%
    merge(covar,.) %>%
    filter(sample %in% test_ids$sample) %>%
    mutate(sample=as.integer(sample)) %>%
    arrange(sample)
  
  model_snp<-readRDS(paste0(pathPred, "xgboost_training_result_", pheno_names$pheno[p], "_hla_snp.RDS"))
  
  #imputed
  print(paste("imputed", pheno_names$pheno_proper[p]))
  imputed_hla<-readRDS(paste0(pathPred, "imputed_hla.rds"))
  
  model_imp<-readRDS(paste0(pathPred, "xgboost_training_result_", pheno_names$pheno[p], "_hla_imputed.RDS"))
  
  pgs_imp<-vroom(paste0(pathLDPred, "ld_pred_scored_", pheno_names$pheno_name_gwas[p], "_without_hla.sscore")) %>%
    rename(sample=IID) %>%
    rename(ld_pred=SCORE1_SUM) %>%
    dplyr::select(c(sample, ld_pred)) %>%
    merge(covar,.) %>%
    merge(., imputed_hla) %>%
    dplyr::select(c("sample", model_imp$feature_names)) %>%
    mutate(sample=as.integer(sample)) %>%
    arrange(sample) %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(sample %in% test_ids$sample)
  
  rm(imputed_hla)
  
  #sequenced
  print(paste("sequenced", pheno_names$pheno_proper[p]))
  sequenced_hla<-readRDS(paste0(pathPred, "sequenced_hla_four_digit.rds"))
  
  model_seq<-readRDS(paste0(pathPred, "xgboost_training_result_four_digit", pheno_names$pheno[p], "_hla_sequenced.RDS"))
  
  pgs_seq<-vroom(paste0(pathLDPred, "ld_pred_scored_", pheno_names$pheno_name_gwas[p], "_without_hla.sscore")) %>%
    rename(sample=IID) %>%
    rename(ld_pred=SCORE1_SUM) %>%
    dplyr::select(c(sample, ld_pred)) %>%
    merge(covar,.) %>%
    merge(., sequenced_hla) %>%
    dplyr::select(c("sample", model_seq$feature_names)) %>%
    mutate(sample=as.integer(sample)) %>%
    arrange(sample) %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(sample %in% test_ids$sample)
  
  #sequenced using imp predictions
  pgs_seq_imp<-vroom(paste0(pathLDPred, "ld_pred_scored_", pheno_names$pheno_name_gwas[p], "_without_hla.sscore")) %>%
    rename(sample=IID) %>%
    rename(ld_pred=SCORE1_SUM) %>%
    dplyr::select(c(sample, ld_pred)) %>%
    merge(covar,.) %>%
    merge(., sequenced_hla) %>%
    dplyr::select(c("sample", "age", "sex", paste0("PC",c(1:10)), "ld_pred",
                    model_imp$feature_names[which(model_imp$feature_names %in% colnames(sequenced_hla))])) %>% 
    mutate(sample=as.integer(sample)) %>%
    arrange(sample) %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(sample %in% test_ids$sample)
  
  for(missing_allele in model_imp$feature_names[which(!model_imp$feature_names %in% 
                                                  c(colnames(sequenced_hla), 
                                                    "sample", 
                                                    "age", 
                                                    "sex", 
                                                    paste0("PC",c(1:10)), 
                                                    "ld_pred"))]){
    pgs_seq_imp<-pgs_seq_imp %>%
      mutate(!!sym(missing_allele) := 0)
  }
  
  pgs_seq_imp<-pgs_seq_imp %>% relocate(colnames(pgs_imp))
  
  rm(sequenced_hla)
  

  #predictions
  print(paste("predictions", pheno_names$pheno_proper[p]))
  X.bm_snp<-as(as.matrix(pgs_snp[,-1]),"dgCMatrix")
  X.bm_imp<-as(as.matrix(pgs_imp[,-1]),"dgCMatrix")
  X.bm_seq<-as(as.matrix(pgs_seq[,-1]),"dgCMatrix")
  X.bm_seq_imp<-as(as.matrix(pgs_seq_imp[,-1]),"dgCMatrix")
  
  xgb_pred_snp<-predict(model_snp, X.bm_snp)
  xgb_pred_imp<-predict(model_imp, X.bm_imp)
  xgb_pred_seq<-predict(model_seq, X.bm_seq)
  xgb_pred_seq_imp<-predict(model_imp, X.bm_seq_imp)
  
  pred<-data.frame(sample=pgs_snp$sample,
                   xgb_pred_snp=xgb_pred_snp,
                   xgb_pred_imp=xgb_pred_imp,
                   xgb_pred_seq=xgb_pred_seq,
                   xgb_pred_seq_imp=xgb_pred_seq_imp) %>%
    merge(.,pheno) %>%
    mutate(Ancestry=NA) %>%
    mutate(Ancestry=ifelse(sample %in% afr, "afr", Ancestry)) %>%
    mutate(Ancestry=ifelse(sample %in% amr, "amr", Ancestry)) %>%
    mutate(Ancestry=ifelse(sample %in% eas, "eas", Ancestry)) %>%
    mutate(Ancestry=ifelse(sample %in% eur, "eur", Ancestry)) %>%
    mutate(Ancestry=ifelse(sample %in% sas, "sas", Ancestry))
  
  write_tsv(pred, paste0(pathPred, "predictions_", pheno_names$pheno[p], ".tsv"))
  
}
