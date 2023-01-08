library(tidyverse)
library(vroom)
library(plotROC)
library(xgboost)
library(ggpubr)
library(pROC)
library(yardstick)

pathPred<-"/path/to/path_out_xgboost/"

#phenotype names
pheno_names<-data.frame(pheno=c("asthma",
                                "ms_demyelinating",
                                "psoriasis",
                                "dm1",
                                "rheumatoid_arthritis",
                                "ulcerative_colitis",
                                "coeliac"),
                        pheno_proper=c("Asthma",
                                       "Multiple Sclerosis and Demyelinating Diseases",
                                       "Psoriasis",
                                       "Type 1 Diabetes Mellitus",
                                       "Rheumatoid Arthritis",
                                       "Ulcerative Colitis",
                                       "Coeliac Disease"),
                        pheno_name_gwas=c("asthma",
                                          "ms",
                                          "psoriasis",
                                          "dm1",
                                          "ra",
                                          "uc",
                                          "celiac"))

pheno_names<-pheno_names[c(1,7,4,2,3,5,6),]

list_plots<-list()

results<-c()
cov_imp_seq<-c()
cov_imp_seq_imp<-c()
for(p in 1:nrow(pheno_names)){
  df<-vroom(paste0(pathPred, "predictions_", pheno_names$pheno[p], ".tsv"))
  
  roc_ldpred<-roc(response=df %>% pull(pheno_names$pheno[p]), predictor=df %>% pull(xgb_pred_snp))
  roc_ldpred_imp<-roc(response=df %>% pull(pheno_names$pheno[p]), predictor=df %>% pull(xgb_pred_imp))
  roc_ldpred_seq<-roc(response=df %>% pull(pheno_names$pheno[p]), predictor=df %>% pull(xgb_pred_seq))
  roc_ldpred_seq_imp<-roc(response=df %>% pull(pheno_names$pheno[p]), predictor=df %>% pull(xgb_pred_seq_imp))
  
  #first print auc and variance
  results<-data.frame(Phenotype=pheno_names$pheno[p],
                      Analyses=c("LDpred", "LDpred+Imputed HLA", "LDpred+Sequenced HLA", "LDpred+Sequenced HLA with Imputed Genes"),
                      AUC=c(roc_ldpred$auc*1, roc_ldpred_imp$auc*1, roc_ldpred_seq$auc*1, roc_ldpred_seq_imp$auc*1),
                      Var=c(var(roc_ldpred), var(roc_ldpred_imp), var(roc_ldpred_seq), var(roc_ldpred_seq_imp))) %>%
    mutate(SE=sqrt(Var)) %>%
    bind_rows(results,.)
    
  cov_imp_seq<-data.frame(Phenotype=pheno_names$pheno[p],
                          cov=pROC::cov(roc_ldpred_imp, roc_ldpred_seq)) %>%
    bind_rows(cov_imp_seq,.)
    
  cov_imp_seq_imp<-data.frame(Phenotype=pheno_names$pheno[p],
                          cov=pROC::cov(roc_ldpred_imp, roc_ldpred_seq_imp)) %>%
    bind_rows(cov_imp_seq_imp,.)
  
  #plot it now
  df_stacked<-data.frame(sample=rep(df$sample,4),
                         xgb_pred=c(df$xgb_pred_snp,
                                    df$xgb_pred_imp,
                                    df$xgb_pred_seq,
                                    df$xgb_pred_seq_imp),
                         truth=rep(df %>% pull(pheno_names$pheno[p]),4),
                         Ancestry=rep(df$Ancestry, 4),
                         Analyses=c(rep("LDpred", nrow(df)),
                                    rep("LDpred+Imputed HLA", nrow(df)),
                                    rep("LDpred+Sequenced HLA", nrow(df)),
                                    rep("LDpred+Sequenced HLA with Imputed Genes", nrow(df))))
  
  colnames(df_stacked)[3]<-pheno_names$pheno[p]
  
  #roc curve
  tmp_roc<-df_stacked %>% 
    ggplot(aes_string(d = pheno_names$pheno[p], m = "xgb_pred", color="Analyses", linetype = "Analyses")) + 
    geom_roc()
  
  df_stacked_final<-df_stacked %>%
    mutate(Analyses=ifelse(Analyses=="LDpred", 
                           paste0("LDpred\nAUC=", format(calc_auc(tmp_roc)$AUC[which(calc_auc(tmp_roc)$Analyses=="LDpred")], digits=3)), 
                           Analyses)) %>%
    mutate(Analyses=ifelse(Analyses=="LDpred+Imputed HLA", 
                           paste0("Imputed HLA\nAUC=", format(calc_auc(tmp_roc)$AUC[which(calc_auc(tmp_roc)$Analyses=="LDpred+Imputed HLA")], digits=3)), 
                           Analyses)) %>%
    mutate(Analyses=ifelse(Analyses=="LDpred+Sequenced HLA", 
                           paste0("Sequenced HLA\nAUC=", format(calc_auc(tmp_roc)$AUC[which(calc_auc(tmp_roc)$Analyses=="LDpred+Sequenced HLA")], digits=3)), 
                           Analyses)) %>%
    mutate(Analyses=ifelse(Analyses=="LDpred+Sequenced HLA with Imputed Genes", 
                           paste0("Sequenced HLA\nwith Imputed Genes\nAUC=", format(calc_auc(tmp_roc)$AUC[which(calc_auc(tmp_roc)$Analyses=="LDpred+Sequenced HLA with Imputed Genes")], digits=3)), 
                           Analyses))
  
  
  list_plots[[p]]<-df_stacked_final %>% 
    ggplot(aes_string(d = pheno_names$pheno[p], m = "xgb_pred", color="Analyses", linetype = "Analyses")) + 
    geom_roc(labels = FALSE, n.cuts=0, size=0.75) +
    scale_color_manual(values=c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"), name=NULL) +
    guides(linetype = "none") +
    theme_bw() +
    theme(legend.key = element_rect(fill = "white", colour = "black"))+
    xlab("") + #xlab("False Positive Rate") +
    ylab("") + #ylab("True Positive Rate") +
    ggtitle(pheno_names$pheno_proper[p])+
    annotate("segment", x=0, y=0, xend=1, yend=1, color="black", linetype="dashed")
  
  ggsave(plot=list_plots[[p]], filename=paste0(pathPred, "final_results/prs_plot_", pheno_names$pheno[p], ".pdf"), width=9.25, height=7)
}

write_tsv(results, paste0(pathPred, "final_results/auc_results_errors.tsv"))
write_tsv(cov_imp_seq, paste0(pathPred, "final_results/cov_imp_seq.tsv"))
write_tsv(cov_imp_seq_imp, paste0(pathPred, "final_results/cov_imp_seq_imp.tsv"))

# write_tsv(AUC_table, file=paste0(pathPred, "final_results/prs_auc_full.tsv"))
ggsave(ggarrange(plotlist=list_plots, nrow=4, ncol=2),filename=paste0(pathPred, "final_results/prs_plot_full.pdf"), width=11, height=11)


#precision recall curves
list_plots_pr<-list()
results_pr<-c()

for(p in 1:nrow(pheno_names)){
  df<-vroom(paste0(pathPred, "predictions_", pheno_names$pheno[p], ".tsv"))
  
  auc_ldpred <- df %>% 
    rename(truth=paste0(pheno_names$pheno[p])) %>%
    arrange(truth) %>% 
    mutate(truth=factor(truth)) %>% 
    pr_auc(data=., truth=truth, xgb_pred_snp, event_level="second") %>%
    pull(.estimate)
  
  auc_ldpred_imp <- df %>% 
    rename(truth=paste0(pheno_names$pheno[p])) %>%
    arrange(truth) %>% 
    mutate(truth=factor(truth)) %>% 
    pr_auc(data=., truth=truth, xgb_pred_imp, event_level="second") %>%
    pull(.estimate)
  
  auc_ldpred_seq <- df %>% 
    rename(truth=paste0(pheno_names$pheno[p])) %>%
    arrange(truth) %>% 
    mutate(truth=factor(truth)) %>% 
    pr_auc(data=., truth=truth, xgb_pred_seq, event_level="second") %>%
    pull(.estimate)
  
  auc_ldpred_seq_imp <- df %>% 
    rename(truth=paste0(pheno_names$pheno[p])) %>%
    arrange(truth) %>% 
    mutate(truth=factor(truth)) %>% 
    pr_auc(data=., truth=truth, xgb_pred_seq_imp, event_level="second") %>%
    pull(.estimate)
  
  results_pr<-data.frame(Phenotype=pheno_names$pheno[p],
                      Analyses=c("LDpred", "LDpred+Imputed HLA", "LDpred+Sequenced HLA", "LDpred+Sequenced HLA with Imputed Genes"),
                      AUC=c(auc_ldpred, auc_ldpred_imp, auc_ldpred_seq, auc_ldpred_seq_imp)) %>%
    bind_rows(results_pr,.)
  
  #now plot results
  df_pr_final<-data.frame(truth=rep(df %>% pull(paste0(pheno_names$pheno[p])),4),
                          method=rep(c(paste0("LDpred\nAUC=", format(auc_ldpred, digits=3)), 
                                       paste0("Imputed HLA\nAUC=", format(auc_ldpred_imp, digits=3)), 
                                       paste0("Sequenced HLA\nAUC=", format(auc_ldpred_seq, digits=3)),
                                       paste0("Sequenced HLA\nwith Imputed Genes\nAUC=", format(auc_ldpred_seq_imp, digits=3))) , 
                                     each=nrow(df)),
                          val=c(df %>% pull(xgb_pred_snp),
                                df %>% pull(xgb_pred_imp),
                                df %>% pull(xgb_pred_seq),
                                df %>% pull(xgb_pred_seq_imp)))
  
  
  list_plots_pr[[p]]<-df_pr_final %>% 
    group_by(method) %>%
    arrange(truth) %>% 
    mutate(truth=factor(truth)) %>% 
    pr_curve(data=., truth=truth, val, event_level="second") %>%
    autoplot() +
    scale_color_manual(values=c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"), name=NULL) +
    guides(linetype = "none") +
    theme_bw() +
    theme(legend.key = element_rect(fill = "white", colour = "black")) +
    ggtitle(pheno_names$pheno_proper[p])+
    xlab("Recall")+
    ylab("Precision")
    
  ggsave(plot=list_plots_pr[[p]], filename=paste0(pathPred, "final_results/pr_prs_plot_", pheno_names$pheno[p], ".pdf"), width=9.25, height=7)
  
}

write_tsv(results_pr, paste0(pathPred, "final_results/pr_auc_results_errors.tsv"))

# write_tsv(AUC_table, file=paste0(pathPred, "final_results/prs_auc_full.tsv"))
ggsave(ggarrange(plotlist=list_plots_pr, nrow=4, ncol=2),filename=paste0(pathPred, "final_results/pr_prs_plot_full.pdf"), width=11, height=11)




#feature importance
for(p in 1:nrow(pheno_names)){
  model_snp<-readRDS(paste0(pathPred, "xgboost_training_result_", pheno_names$pheno[p], "_hla_snp.RDS"))
  model_imp<-readRDS(paste0(pathPred, "xgboost_training_result_", pheno_names$pheno[p], "_hla_imputed.RDS"))
  model_seq<-readRDS(paste0(pathPred, "xgboost_training_result_", pheno_names$pheno[p], "_hla_sequenced.RDS"))
  
  #snp only
  importance_snp<- xgb.importance(model = model_snp) %>%
    mutate(Phenotype=pheno_names$pheno[p]) %>%
    mutate(Model="LDpred")
  
  write_tsv(importance_snp, paste0(pathPred, "final_results/prs_feature_importance_", pheno_names$pheno[p], "_LDpred.tsv.gz"))
  
  pdf(file=paste0(pathPred, "final_results/prs_feature_importance_", pheno_names$pheno[p], "_LDpred.pdf"),
      width=10,
      height=6)
  xgb.plot.importance(importance_snp[1:nrow(importance_snp),], plot=TRUE, main=paste0(pheno_names$pheno_proper[p], ": LDpred\nXGboost Feature Importance"))
  dev.off()
  
  #imputed
  importance_imp<- xgb.importance(model = model_imp) %>%
    mutate(Phenotype=pheno_names$pheno[p]) %>%
    mutate(Model="LDpred_Imputed_HLA")
  
  write_tsv(importance_imp, paste0(pathPred, "final_results/prs_feature_importance_", pheno_names$pheno[p], "_LDpred_hla_imputed.tsv.gz"))
  
  pdf(file=paste0(pathPred, "final_results/prs_feature_importance_", pheno_names$pheno[p], "_LDpred_hla_imputed.pdf"),
      width=10,
      height=6)
  xgb.plot.importance(importance_imp[1:nrow(importance_snp),], plot=TRUE, main=paste0(pheno_names$pheno_proper[p], ": LDpred+Imputed HLA\nXGboost Feature Importance"))
  dev.off()
  
  #sequenced
  importance_seq<- xgb.importance(model = model_seq) %>%
    mutate(Phenotype=pheno_names$pheno[p]) %>%
    mutate(Model="LDpred_Sequenced_HLA")
  
  write_tsv(importance_seq, paste0(pathPred, "final_results/prs_feature_importance_", pheno_names$pheno[p], "_LDpred_hla_sequenced.tsv.gz"))
  
  pdf(file=paste0(pathPred, "final_results/prs_feature_importance_", pheno_names$pheno[p], "_LDpred_hla_sequenced.pdf"),
      width=10,
      height=6)
  xgb.plot.importance(importance_seq[1:nrow(importance_snp),], plot=TRUE, main=paste0(pheno_names$pheno_proper[p], ": LDpred+Sequenced HLA\nXGboost Feature Importance"))
  dev.off()
}
