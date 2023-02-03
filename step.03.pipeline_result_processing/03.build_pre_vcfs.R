setwd("/project/richards/guillaume.butler-laporte/HLA/ukb_wes")

library(tidyverse)
library(vroom)

#options(warn=2)

#dp_threshold<-10


hla_alleles<-as.list(c((29941260+29945884)/2,
                       (31353872+31367067)/2,
                       (31268749+31272130)/2,
                       (30489509+30494194)/2,
                       (29722775+29738528)/2,
                       (29826967+29831125)/2,
                       (29887752+29890482)/2,
                       (30006606+30009539)/2,
                       (29926459+29929232)/2,
                       (30259625+30261703)/2,
                       (29792234+29793136)/2,
                       (29941260+31367067)/2,
                       (32948613+32969094)/2,
                       (32934629+32941028)/2,
                       (33004182+33009591)/2,
                       (32812763+32820466)/2,
                       (33064569+33080775)/2,
                       (33091753+33097295)/2,
                       (33075990+33089696)/2,
                       (32628179+32647062)/2,
                       (32659467+32668383)/2,
                       (32439878+32445046)/2,
                       (32577902+32589848)/2,
                       (32488465+32502952)/2,
                       (32449765+32462852)/2,
                       (32542598+32557561)/2,
                       (32517353+32530287)/2,
                       (32552713+32560022)/2,
                       (32627102+32646281)/2,
                       (32589215+32591752)/2,
                       (32459821+32473500)/2))

names(hla_alleles)<-c("A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V", "Y",
                      "DMA", "DMB", "DOA", "DOB", "DPA1", "DPA2", "DPB1", "DQA1", 
                      "DQB1", "DRA", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB6",
                      "DRB7", "DRB8", "DRB9")  


all_alleles_six<-c()
all_alleles_four<-c()
all_alleles_two<-c()
for(folder in 10:60){
  hla_df<-vroom(paste0("calls/hla_df_batch_qced_", folder, ".tsv.gz"), col_types = cols(.default = "c")) %>% as.data.frame()
  all_alleles_six<-stack(hla_df[,-1]) %>%
    filter(!is.na(values)) %>%
    dplyr::select(values) %>%
    filter(values!="0") %>%
    bind_rows(all_alleles_six,.) %>%
    distinct()
    
  hla_df_four <- lapply(hla_df, function(x) ifelse(str_count(x,":")>1, sub(":[0-9A-Z]*$", "", x), x)) %>%
    as.data.frame()
  all_alleles_four<-stack(hla_df_four[,-1]) %>%
    filter(!is.na(values)) %>%
    dplyr::select(values) %>%
    filter(values!="0") %>%
    bind_rows(all_alleles_four,.) %>%
    distinct()
    
  hla_df_two <- hla_df %>%
    mutate(ID=as.numeric(ID)) %>%
    arrange(ID) %>%
    mutate_if(is_character, ~str_extract(string=., pattern="^[A-Z0-9]*\\*[0-9]*")) %>%
    select(where(~n_distinct(.) > 1)) %>%
    replace(is.na(.),"0")  
  all_alleles_two<-stack(hla_df_two[,-1]) %>%
    filter(!is.na(values)) %>%
    dplyr::select(values) %>%
    filter(values!="0") %>%
    bind_rows(all_alleles_two,.) %>%
    distinct()
}


for(folder in 10:60){
#for(folder in 23){
  print(folder)
  
  hla_df<-vroom(paste0("calls/hla_df_batch_qced_", folder, ".tsv.gz"), col_types = cols(.default = "c")) %>% as.data.frame()
  
  #remove the HLA prefix if needed
  #hla_df <- data.frame(lapply(hla_df, function(x){gsub("HLA-", "", x)}))
  
  ######## transform into vcf #########
  # first obtain all alleles, note that some alleles do not have 3 fields in the database. These are hence reported as (e.g. for HLA-A) A*XX:XX and not as A*XX:XX:01
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
  
  #for DRB345 see here: https://www-sciencedirect-com.proxy3.library.mcgill.ca/science/article/pii/S1357272520301990  
  #a<-26
  #i<-1
  for(s in 1:nrow(hla_df)){
    #print(s)
    for(a in 2:ncol(hla_df)){
      if(hla_df[s,a]!=0){
        if(gt_six[hla_df[s,a],s]=="0/0"){
          gt_six[hla_df[s,a],s]<-"0/1"
        } else {
          gt_six[hla_df[s,a],s]<-"1/1"
        }
      } 
    }
  }

  hla_vcf_six <- cbind(hla_vcf_six_tmp, gt_six) %>%
    mutate(tmp_hla=ID) %>%
    separate(tmp_hla, into=c("tmp_gene", "tmp_allele"), sep="[*]") %>%
    mutate(POS=hla_alleles[tmp_gene]) %>%
    dplyr::select(-c(tmp_gene,tmp_allele)) %>%
    mutate(POS=as.numeric(POS)) %>%
    mutate(POS=POS+c(1:nrow(.)) %>%
    dplyr::arrange(POS)
  
  #save result
  vroom_write(hla_vcf_six, paste0("vcf/pre_vcf/hla_six_digit_", folder, ".pre_vcf.tsv.gz"))
   
  
  #now four digits
  hla_df_four <- lapply(hla_df, function(x) ifelse(str_count(x,":")>1, sub(":[0-9A-Z]*$", "", x), x)) %>%
    as.data.frame()
  
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
      if(hla_df_four[s,a]!=0){
        if(gt_four[hla_df_four[s,a],s]=="0/0"){
          gt_four[hla_df_four[s,a],s]<-"0/1"
        } else {
          gt_four[hla_df_four[s,a],s]<-"1/1"
        }
      } 
    } 
  }

  hla_vcf_four <- cbind(hla_vcf_four_tmp, gt_four) %>%
    mutate(tmp_hla=ID) %>%
    separate(tmp_hla, into=c("tmp_gene", "tmp_allele"), sep="[*]") %>%
    mutate(POS=hla_alleles[tmp_gene]) %>%
    dplyr::select(-c(tmp_gene,tmp_allele)) %>%
    mutate(POS=as.numeric(POS)) %>%
    mutate(POS=POS+c(1:nrow(.)) %>%
    dplyr::arrange(POS)
  
  #save result
  vroom_write(hla_vcf_four, paste0("vcf/pre_vcf/hla_four_digit_",folder, ".pre_vcf.tsv.gz"))
  
  
  #now two digits
  hla_df_two <- hla_df %>%
    mutate(ID=as.numeric(ID)) %>%
    arrange(ID) %>%
    mutate_if(is_character, ~str_extract(string=., pattern="^[A-Z0-9]*\\*[0-9]*")) %>%
    select(where(~n_distinct(.) > 1)) %>%
    replace(is.na(.),"0")
  
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
      if(hla_df_two[s,a]!=0){
        if(gt_two[hla_df_two[s,a],s]=="0/0"){
          gt_two[hla_df_two[s,a],s]<-"0/1"
        } else {
          gt_two[hla_df_two[s,a],s]<-"1/1"
        }
      } 
    }
  }

  
  hla_vcf_two <- cbind(hla_vcf_two_tmp, gt_two) %>%
    mutate(tmp_hla=ID) %>%
    separate(tmp_hla, into=c("tmp_gene", "tmp_allele"), sep="[*]") %>%
    mutate(POS=hla_alleles[tmp_gene]) %>%
    dplyr::select(-c(tmp_gene,tmp_allele)) %>%
    mutate(POS=as.numeric(POS)) %>%
    mutate(POS=POS+c(1:nrow(.)) %>%
    dplyr::arrange(POS)
  
  #save result
  vroom_write(hla_vcf_two, paste0("vcf/pre_vcf/hla_two_digit_",folder, ".pre_vcf.tsv.gz"))
  
}
