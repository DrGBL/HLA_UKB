library(tidyverse)
library(vroom)
library(xgboost)
library(rBayesianOptimization)

pathIn<-"/path/to/pgs_inputs/"
pathOut<-"/path/to/path_out_xgboost/"
pathFreqs<-"/path/HLA/allele/frequencies/"
pathPCs<-"/path/to/pcs/"
pathPhenoCovar<-"/path/to/regenie_step_1_inputs/"
pathLDPred<-"/path/to/ld_pred_output/"

#Extract allele frequencies and choose threshold
freq_thresh<-0.0001

#the file here is from plink's freq counts function. Available in step.05.misc section
all_freq<-read.csv(paste0(pathFreqs,"freq_four_digit_all.frq.gz"), sep="")%>%
    mutate(AF=ifelse(A1=="A", MAF, 1-MAF)) %>%
    dplyr::select(-c(CHR, A1, A2, MAF, NCHROBS)) %>%
    filter(AF>freq_thresh) %>%
    pull(SNP)

#extract list of samples ids
#this is a column with UKB identifiers for participants with WES (with header "sample")
sequenced_ids<-read_tsv(paste0(pathIn, "sequenced_sample_ids.tsv.gz")) %>% pull(sample)

#extract PCs
#obtained from projecting on 1000G, like when we assigned genetic ancestries. Please see main paper for details.
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
  filter(sample %in% sequenced_ids)

rm(eigenvec)

#split training and testing samples
set.seed(1234)
test_id<-covar$sample[sample.int(n=nrow(covar), size=nrow(covar)/5)]
train_id<-covar %>% filter(!(sample %in% test_id)) %>% pull(sample)

covar %>% dplyr::select(sample) %>% filter(sample %in% test_id) %>% write_tsv(., paste0(pathOut, "test_ids.tsv"))
covar %>% dplyr::select(sample) %>% filter(sample %in% train_id) %>% write_tsv(., paste0(pathOut, "train_ids.tsv"))

pheno_name_pheno_file<-c("dm1", "psoriasis", "asthma", "coeliac", "ms_demyelinating",  "rheumatoid_arthritis", "ulcerative_colitis")
pheno_name_gwas<-c("dm1", "psoriasis", "asthma", "celiac", "ms",  "ra", "uc")

for(ana in c("hla_sequenced", "hla_imputed", "hla_snp")){
  if(ana=="hla_imputed"){
    #this is a file with first column (header "ID") is all the imputed HLA alleles, and the next columns are 0, 1, or 2 based on an additive model (and the headers are the participants identifiers).
    imputed_hla<-readRDS(paste0(pathOut, "imputed_hla.rds"))
  }
  if(ana=="hla_sequenced"){
    #this is the same file format as hla_imputed, but for sequenced HLA alleles (so way more rows).
    sequenced_hla<-readRDS(paste0(pathOut, "sequenced_hla_four_digit.rds")) %>%
      dplyr::select(c(sample, all_of(all_freq)))
  }
  
  for(p in 1:length(pheno_name_gwas)){
    pheno<-vroom(paste0(pathPhenoCovar,
                        "pheno_file.tsv")) %>%
      rename(sample=FID) %>%
      dplyr::select(-IID) %>%
      mutate(sample=as.integer(sample)) %>%
      filter(sample %in% covar$sample) %>%
      arrange(sample) %>%
      dplyr::select(c(sample, pheno_name_pheno_file[p])) %>%
      filter(sample %in% train_id)
    
    if(ana=="hla_snp"){
      pgs<-vroom(paste0(pathLDPred, "ld_pred_scored_", pheno_name_gwas[p], "_with_hla.sscore")) %>%
        rename(sample=IID) %>%
        rename(ld_pred=SCORE1_SUM) %>%
        dplyr::select(c(sample, ld_pred)) %>%
        merge(covar,.) %>%
        filter(sample %in% train_id)
    }
    if(ana=="hla_imputed"){
      pgs<-vroom(paste0(pathLDPred, "ld_pred_scored_", pheno_name_gwas[p], "_without_hla.sscore")) %>%
        rename(sample=IID) %>%
        rename(ld_pred=SCORE1_SUM) %>%
        dplyr::select(c(sample, ld_pred)) %>%
        merge(covar,.) %>%
        merge(., imputed_hla) %>%
        dplyr::select(where(~ any(. != 0))) %>%
        mutate(sample=as.integer(sample)) %>%
        arrange(sample) %>%
        filter(sample %in% train_id) %>%
        mutate(across(everything(), as.numeric))
    }
    if(ana=="hla_sequenced"){
      pgs<-vroom(paste0(pathLDPred, "ld_pred_scored_", pheno_name_gwas[p], "_without_hla.sscore")) %>%
        rename(sample=IID) %>%
        rename(ld_pred=SCORE1_SUM) %>%
        dplyr::select(c(sample, ld_pred)) %>%
        merge(covar,.) %>%
        merge(., sequenced_hla) %>%
        dplyr::select(where(~ any(. != 0))) %>%
        mutate(sample=as.integer(sample)) %>%
        arrange(sample) %>%
        filter(sample %in% train_id) %>%
        mutate(across(everything(), as.character)) %>%
        mutate(across(everything(), as.numeric))
    }
    
    #XGboost starts here
    #transform to dgCMatrix format
    X.bm_train<-as(as.matrix(pgs[,-1]),"dgCMatrix")
    
    Y_train<-pheno %>% pull(pheno_name_pheno_file[p])
    
    dtrain <- xgb.DMatrix(data = X.bm_train, label = Y_train)
    
    
    #finding best hyperparameters with bayesian optimization (https://stackoverflow.com/questions/35050846/xgboost-in-r-how-does-xgb-cv-pass-the-optimal-parameters-into-xgb-train)
    cv_folds <- KFold(Y_train, nfolds = 5, stratified = FALSE, seed = 1004)
    
    
    xgb_cv_bayes <- function(nround,max.depth, min_child_weight, subsample,eta,gamma,colsample_bytree,max_delta_step) {
      param<-list(booster = "gbtree",
                  max_depth = max.depth,
                  min_child_weight = min_child_weight,
                  eta=eta,gamma=gamma,
                  subsample = subsample, colsample_bytree = colsample_bytree,
                  max_delta_step=max_delta_step,
                  lambda = 1, alpha = 0,
                  objective = "binary:logistic",
                  eval_metric = "logloss")
      cv <- xgb.cv(params = param, data = dtrain, folds = cv_folds,nthread = 10,nrounds = 10000,early_stopping_rounds = 30, verbose = FALSE, maximize=FALSE)
      
      list(Score = -cv$evaluation_log$test_logloss_mean[cv$best_iteration],
           Pred=cv$best_iteration)
      # we don't need cross-validation prediction and we need the number of rounds.
      # a workaround is to pass the number of rounds(best_iteration) to the Pred, which is a default parameter in the rbayesianoptimization library.
    } 
    
    OPT_Res <- BayesianOptimization(xgb_cv_bayes,
                                    bounds = list(max.depth =c(3L, 10L),min_child_weight = c(1L, 40L),
                                                  subsample = c(0.6, 0.9),
                                                  eta=c(0.01,0.3),gamma = c(0.0, 0.2),
                                                  colsample_bytree=c(0.5,0.8),max_delta_step=c(1L,10L)),
                                    init_grid_dt = NULL, init_points = 10, n_iter = 10,
                                    acq = "ucb", kappa = 2.576, eps = 0.0)
    
    best_param <- list(
      booster = "gbtree",
      eval.metric = "logloss",
      objective = "binary:logistic",
      max_depth = OPT_Res$Best_Par["max.depth"],
      eta = OPT_Res$Best_Par["eta"],
      gamma = OPT_Res$Best_Par["gamma"],
      subsample = OPT_Res$Best_Par["subsample"],
      colsample_bytree = OPT_Res$Best_Par["colsample_bytree"],
      min_child_weight = OPT_Res$Best_Par["min_child_weight"],
      max_delta_step = OPT_Res$Best_Par["max_delta_step"])
    
    
    # number of rounds should be tuned using CV
    #https://www.hackerearth.com/practice/machine-learning/machine-learning-algorithms/beginners-tutorial-on-xgboost-parameter-tuning-r/tutorial/
    # However, nrounds can not be directly derivied from the bayesianoptimization function
    # Here, OPT_Res$Pred, which was supposed to be used for cross-validation, is used to record the number of rounds
    nrounds=OPT_Res$Pred[[which.max(OPT_Res$History$Value)]]
    xgb_model <- xgb.train (params = best_param, data = dtrain, nrounds = nrounds)
    
    #save results
    saveRDS(xgb_model, paste0(pathOut, "xgboost_training_result_four_digit", pheno_name_pheno_file[p], "_", ana, ".RDS"))
    rm(pgs)
  }
  
  if(ana=="hla_imputed"){
    rm(imputed_hla)
  }
  if(ana=="hla_sequenced"){
    rm(sequenced_hla)
  }
  
}




