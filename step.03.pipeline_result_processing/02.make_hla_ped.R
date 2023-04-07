library(tidyverse)
library(vroom)

qc_threshold<-20

setwd("/your/local/folder")

pathCalls<-"calls/"
pathCov<-"QC/"

genes<-c("A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V", "Y", "DMA", "DMB", "DOA", "DOB", "DPA1", "DPA2", "DPB1", "DQA1", "DQB1", "DRA", "DRB1", "DRB3", "DRB4", "DRB5", "DRB6", "DRB7", "DRB2", "DRB8", "DRB9")

for(folder in 10:60){
  print(folder)
  coverage<-vroom(paste0(pathCov, "hla_coverage_batch_", folder, ".tsv.gz"), col_types = cols(.default = "c")) %>% mutate(ID=as.numeric(ID)) %>% arrange(ID) %>% replace(is.na(.),"0") 
  
  hla_calls_all_folder<-vroom(paste0(pathCalls, "hla_df_batch_", folder, ".tsv.gz"),
                              col_types = cols(.default = "c")) %>% mutate(ID=as.numeric(ID)) %>% arrange(ID) %>%
    mutate(across(everything(), gsub, pattern="HLA-", replacement="")) %>%
    replace(is.na(.),"0") 
  
  for(i in 1:nrow(hla_calls_all_folder)){
    if(i %% 500 ==0) { print(i) }
    for(g in genes){
      DP<-as.numeric(pull(coverage[i,g]))
      if(DP<qc_threshold){
        hla_calls_all_folder[i, paste0(g, "_1")]<-"0"
        hla_calls_all_folder[i, paste0(g, "_2")]<-"0"
      } else if (g %in% c("DRB3", "DRB4", "DRB5")) {
        if(coverage_drb[i,paste0(g, "_1")] < coverage_drb[i,paste0(g, "_2")]){
          tmp_a<-hla_calls_all_folder[i,paste0(g, "_1")]
          hla_calls_all_folder[i,paste0(g, "_1")]<-hla_calls_all_folder[i,paste0(g, "_2")]
          hla_calls_all_folder[i,paste0(g, "_2")]<-tmp_a
        }
      }
    }
  }
  
  
  tmp<-hla_calls_all_folder %>% 
    dplyr::select(c(DRB1_1,
                    DRB1_2,
                    DRB3_1,
                    DRB3_2,
                    DRB4_1,
                    DRB4_2,
                    DRB5_1,
                    DRB5_2))
  
  for(i in 1:nrow(tmp)){
    if(i %% 500 ==0) { print(i) }
    drb_df<-data.frame(gene=c("DRB3", "DRB4", "DRB5"),
                       drb1_cov=c(as.numeric(pull(coverage[i,"DRB1"])), as.numeric(pull(coverage[i,"DRB1"])), as.numeric(pull(coverage[i,"DRB1"]))),
                       cov=c(as.numeric(pull(coverage[i,"DRB3"])), as.numeric(pull(coverage[i,"DRB4"])), as.numeric(pull(coverage[i,"DRB5"])))) %>%
      mutate(hom_dr=ifelse(cov>0.6*drb1_cov & cov>qc_threshold, TRUE, FALSE)) %>%
      mutate(het_dr=ifelse(cov>0.3*drb1_cov & cov>qc_threshold, TRUE, FALSE))
    
    check_hom<-drb_df %>% filter(hom_dr)
    if(nrow(check_hom)>0){
      if(nrow(check_hom)==1){
        gene_kept<-check_hom %>% filter(cov==max(cov)) %>% pull(gene)
        gene_removed<-paste0(rep(drb_df %>% filter(gene!=gene_kept) %>% pull(gene), each=2), c("_1", "_2"))
        tmp[i,gene_removed]<-"0"
      } else {
        if(nrow(check_hom==2)){
          gene_kept<-check_hom %>% filter(hom_dr) %>% pull(gene)
          tmp[i,paste0(gene_kept, "_2")]<-"0"
          
          gene_removed<-paste0(rep(drb_df %>% filter(cov==min(cov)) %>% pull(gene), each=2), c("_1", "_2"))
          tmp[i,gene_removed]<-"0"
        } else {
          stop("Three DRB genes above 0.6*coverage DRB1")
        }
      }
    } else {
      check_het<-drb_df %>% filter(het_dr)
      if(nrow(check_het)==1){
        gene_kept<-check_het %>% filter(het_dr) %>% pull(gene)
        tmp[i,paste0(gene_kept, "_2")]<-"0"
        
        gene_removed<-paste0(rep(drb_df %>% filter(!(gene %in% gene_kept)) %>% pull(gene), each=2), c("_1", "_2"))
        tmp[i,gene_removed]<-"0"
      }
      if(nrow(check_het)==2){
        gene_kept<-check_het %>% filter(het_dr) %>% pull(gene)
        tmp[i,paste0(gene_kept, "_2")]<-"0"
        
        gene_removed<-paste0(rep(drb_df %>% filter(cov==min(cov)) %>% pull(gene), each=2), c("_1", "_2"))
        tmp[i,gene_removed]<-"0"
      }
      if(nrow(check_het)==3){
        gene_removed<-check_het %>% filter(cov==min(cov)) %>% pull(gene)
        tmp[i, paste0(gene_removed, "_2")]<-"0"
        
        gene_kept<-check_het %>% filter(gene!=gene_removed) %>% pull(gene)
        tmp[i,paste0(gene_kept, "_2")]<-"0"
      }
    }
  }
  
  hla_calls_all_folder<-hla_calls_all_folder %>% 
    dplyr::select(-c(DRB1_1,
                     DRB1_2,
                     DRB3_1,
                     DRB3_2,
                     DRB4_1,
                     DRB4_2,
                     DRB5_1,
                     DRB5_2)) %>%
    bind_cols(.,tmp)
  
  write_delim(hla_calls_all_folder, paste0(pathCalls, "hla_df_batch_qced_", folder, ".tsv.gz"), delim="\t")
}
