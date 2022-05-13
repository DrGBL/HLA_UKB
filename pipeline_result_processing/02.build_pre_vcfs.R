setwd("/your/local/folder")

library(tidyverse)
library(vroom)

dp_threshold<-10 #read depth threshold used in our analyses

hla_alleles<-as.list(c(1:31))
names(hla_alleles)<-c("A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V", "Y",
                      "DMA", "DMB", "DOA", "DOB", "DPA1", "DPA2", "DPB1", "DQA1", 
                      "DQB1", "DRA", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB6",
                      "DRB7", "DRB8", "DRB9")  

for(folder in 10:60){
  
  hla_df<-vroom(paste0("calls/hla_df_batch_", folder, ".tsv.gz"), col_types = cols(.default = "c"))
  hla_qc<-vroom(paste0("QC/hla_coverage_batch_", folder, ".tsv.gz"), col_types = cols(.default = "c"))
  
  #remove the HLA prefix if needed
  hla_df <- data.frame(lapply(hla_df, function(x){gsub("HLA-", "", x)}))
  
  ######## transform into vcf #########
  # first obtain all alleles
  all_alleles_six<-stack(hla_df[,-1]) %>%
    filter(!is.na(values)) %>%
    dplyr::select(values) %>%
    distinct() %>%
    mutate(count_colon=str_count(values, ":"))
    
  all_alleles_six_or_more<-all_alleles_six %>% filter(count_colon>1)

  hla_vcf_six_tmp<-data.frame(CHROM=rep(6, nrow(all_alleles_six)),
                              POS=NA,
                              ID=all_alleles_six$values,
                              REF=rep("C", nrow(all_alleles_six)),
                              ALT=rep("A", nrow(all_alleles_six)),
                              QUAL=rep(".", nrow(all_alleles_six)),
                              FILTER=rep(".", nrow(all_alleles_six)),
                              INFO=rep(".", nrow(all_alleles_six)),
                              FORMAT=rep("GT", nrow(all_alleles_six)))
  colnames(hla_vcf_six_tmp)[1]<-"#CHROM"
  
  gt_six<-matrix("0/0", nrow=nrow(all_alleles_six), ncol=nrow(hla_df))
  colnames(gt_six)<-hla_df$ID
  rownames(gt_six)<-all_alleles_six$values
  
  for(s in 1:nrow(hla_df)){
    #print(s)
    for(a in 2:ncol(hla_df)){
      #print(s)
      #print(a)
      if(!is.na(hla_df[s,a])){
        gene<-sub("\\*.*","",hla_df[s,a])
        coverage<-hla_qc[s,gene]
        if(coverage>=dp_threshold){
          if(gt_six[hla_df[s,a],s]=="0/0"){
            gt_six[hla_df[s,a],s]<-"0/1"
          } else {
            gt_six[hla_df[s,a],s]<-"1/1"
          }
        }
      } 
    }
  }

  hla_vcf_six <- cbind(hla_vcf_six_tmp, gt_six) %>%
    mutate(tmp_hla=ID) %>%
    separate(tmp_hla, into=c("tmp_gene", "tmp_allele"), sep="[*]") %>%
    mutate(POS=paste0(hla_alleles[tmp_gene], str_remove_all(tmp_allele, "[:]"))) %>%
    dplyr::select(-c(tmp_gene,tmp_allele)) %>%
    mutate(POS=str_remove_all(POS, "[NLSCAQ]")) %>%
    mutate(POS=as.numeric(POS)) %>%
    dplyr::arrange(POS) %>%
    filter(ID %in% all_alleles_six_or_more$values)
  
  #save result
  vroom_write(hla_vcf_six, paste0("vcf/pre_vcf/hla_six_digit_", folder, ".pre_vcf.tsv.gz"))
  
  #now four digits
  hla_df_four <- lapply(hla_df, function(x) ifelse(str_count(x,":")>1, sub(":[0-9A-Z]*$", "", x), x)) %>%
    as.data.frame()

  all_alleles_four<-stack(hla_df_four[,-1]) %>%
    filter(!is.na(values)) %>%
    dplyr::select(values) %>%
    distinct()
  
  hla_vcf_four_tmp<-data.frame(CHROM=rep(6, nrow(all_alleles_four)),
                               POS=NA,
                               ID=all_alleles_four$values,
                               REF=rep("C", nrow(all_alleles_four)),
                               ALT=rep("A", nrow(all_alleles_four)),
                               QUAL=rep(".", nrow(all_alleles_four)),
                               FILTER=rep(".", nrow(all_alleles_four)),
                               INFO=rep(".", nrow(all_alleles_four)),
                               FORMAT=rep("GT", nrow(all_alleles_four)))
  colnames(hla_vcf_four_tmp)[1]<-"#CHROM"
  
  gt_four<-matrix("0/0", nrow=nrow(all_alleles_four), ncol=nrow(hla_df_four))
  colnames(gt_four)<-hla_df_four$ID
  rownames(gt_four)<-all_alleles_four$values
  
  
    
  for(s in 1:nrow(hla_df_four)){
    #print(s)
    for(a in 2:ncol(hla_df_four)){
      #print(s)
      #print(a)
      if(!is.na(hla_df_four[s,a])){
        gene<-sub("\\*.*","",hla_df_four[s,a])
        coverage<-hla_qc[s,gene]
        if(coverage>=dp_threshold){
          if(gt_four[hla_df_four[s,a],s]=="0/0"){
            gt_four[hla_df_four[s,a],s]<-"0/1"
          } else {
            gt_four[hla_df_four[s,a],s]<-"1/1"
          }
        }
      } 
    }
  }

  hla_vcf_four <- cbind(hla_vcf_four_tmp, gt_four) %>%
    mutate(tmp_hla=ID) %>%
    separate(tmp_hla, into=c("tmp_gene", "tmp_allele"), sep="[*]") %>%
    mutate(POS=paste0(hla_alleles[tmp_gene], str_remove_all(tmp_allele, "[:]"))) %>%
    dplyr::select(-c(tmp_gene,tmp_allele)) %>%
    mutate(POS=str_remove_all(POS, "[NLSCAQ]")) %>%
    mutate(POS=as.numeric(POS)) %>%
    dplyr::arrange(POS)
  
  #save result
  vroom_write(hla_vcf_four, paste0("vcf/pre_vcf/hla_four_digit_",folder, ".pre_vcf.tsv.gz"))
  
  
  #now two digits
  hla_df_two <- hla_df %>%
    mutate(ID=as.numeric(ID)) %>%
    arrange(ID) %>%
    mutate_if(is_character, ~str_extract(string=., pattern="^[A-Z0-9]*\\*[0-9]*")) %>%
    select(where(~n_distinct(.) > 1))
  
  
  all_alleles_two<-stack(hla_df_two[,-1]) %>%
    filter(!is.na(values)) %>%
    dplyr::select(values) %>%
    distinct()
  
  hla_vcf_two_tmp<-data.frame(CHROM=rep(6, nrow(all_alleles_two)),
                              POS=NA,
                              ID=all_alleles_two$values,
                              REF=rep("C", nrow(all_alleles_two)),
                              ALT=rep("A", nrow(all_alleles_two)),
                              QUAL=rep(".", nrow(all_alleles_two)),
                              FILTER=rep(".", nrow(all_alleles_two)),
                              INFO=rep(".", nrow(all_alleles_two)),
                              FORMAT=rep("GT", nrow(all_alleles_two)))
  colnames(hla_vcf_two_tmp)[1]<-"#CHROM"
  
  gt_two<-matrix("0/0", nrow=nrow(all_alleles_two), ncol=nrow(hla_df_two))
  colnames(gt_two)<-hla_df_two$ID
  rownames(gt_two)<-all_alleles_two$values
  
  
  
  for(s in 1:nrow(hla_df_two)){
    #print(s)
    for(a in 2:ncol(hla_df_two)){
      #print(s)
      #print(a)
      if(!is.na(hla_df_two[s,a])){
        gene<-sub("\\*.*","",hla_df_two[s,a])
        coverage<-hla_qc[s,gene]
        if(coverage>=dp_threshold){
          if(gt_two[hla_df_two[s,a],s]=="0/0"){
            gt_two[hla_df_two[s,a],s]<-"0/1"
          } else {
            gt_two[hla_df_two[s,a],s]<-"1/1"
          }
        }
      } 
    }
  }
  
  hla_vcf_two <- cbind(hla_vcf_two_tmp, gt_two) %>%
    mutate(tmp_hla=ID) %>%
    separate(tmp_hla, into=c("tmp_gene", "tmp_allele"), sep="[*]") %>%
    mutate(POS=paste0(hla_alleles[tmp_gene], str_remove_all(tmp_allele, "[:]"))) %>%
    dplyr::select(-c(tmp_gene,tmp_allele)) %>%
    mutate(POS=str_remove_all(POS, "[NLSCAQ]")) %>%
    mutate(POS=as.numeric(POS)) %>%
    dplyr::arrange(POS)
  
  #save result
  vroom_write(hla_vcf_two, paste0("vcf/pre_vcf/hla_two_digit_",folder, ".pre_vcf.tsv.gz"))

}
