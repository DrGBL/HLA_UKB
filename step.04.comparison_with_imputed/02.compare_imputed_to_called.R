library(tidyverse)
library(vroom)

path_anc<-"path_to_anc_folder"

#ancestries
afr<-scan(paste0(path_anc, "ukb.afrIDsPCA.txt"))
amr<-scan(paste0(path_anc, "ukb.amrIDsPCA.txt"))
eas<-scan(paste0(path_anc, "ukb.easIDsPCA.txt"))
eur<-scan(paste0(path_anc, "ukb.eurIDsPCA.txt"))
sas<-scan(paste0(path_anc, "ukb.sasIDsPCA.txt"))

#again go to the directory where the folders "calls" and "QC" are located
setwd("/your/local/folder")

if(!dir.exists("comparisons_with_imputed")){
  dir.create("comparisons_with_imputed")
}

#path to munged imputed HLA data (from step 01)
path_munged_imputed<-"/path/to/munged/imputed"

hla_imputed_munged<-vroom(path_munged_imputed, col_types = cols(.default = "c"))


#option choices, useful for later
options<-c("Both uncalled, one allele imputed", 
  "Both uncalled, both alleles imputed",
  "Both uncalled and unimputed",
  "One match, other uncalled and unimputed",
  "One match, other unimputed",
  "One match, other uncalled but imputed",
  "One match, one mismatch imputation",
  "Two matches",
  "One called, both unimputed",
  "One called, one unimputed, no match",
  "One called, both imputed, no match",
  "Both called, both unimputed",
  "Both called, one unimputed, no match",
  "Both called, both imputed, no match")

#alleles in the imputation panel by ukb https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=2182
imputed_alleles<-data.frame(alleles=scan("imputed_alleles_link_file.txt", what=character())) %>%
  mutate(alleles=ifelse(nchar(str_extract(alleles, "[0-9]*$"))==3,
                        str_replace(alleles,"_","_0"),
                        alleles)) %>%
  mutate(alleles=paste0(str_extract(alleles,
                                    "^[A-Z0-9]*"),
                        "*",
                        str_extract(alleles,
                                    "_[0-9][0-9]"),
                        ":",
                        str_extract(alleles,
                                    "[0-9][0-9]$"))) %>%
  mutate(alleles=str_replace(alleles,"_","")) %>%
  filter(!str_detect(alleles,"99:01"))


#again those are folders (10 to 60) containing the UKB participants
for(folder in 10:60){

  hla_df<-vroom(paste0("calls/hla_df_batch_qced_", folder, ".tsv.gz"), col_types = cols(.default = "c"))

  #now compare to current batch, assuming that all "99:01" alleles are in fact unimputed
  batch_imputed<-hla_imputed_munged %>% 
    filter(ID %in% hla_df$ID)
  batch_imputed[is.na(batch_imputed)]<-"Unimputed"
  batch_imputed<-lapply(batch_imputed, function(x) gsub("[0-9A-Z]*\\*99:01","Unimputed",x)) %>%
    as.data.frame()

  #this trims the allele to two fields, and removes the expression suffix (if present)
  hla_df_comparator_four<-lapply(hla_df, function(x) ifelse(str_count(x,":")>1, sub(":[0-9A-Z]*$", "", x), x)) %>%
    as.data.frame() %>%
    filter(ID %in% batch_imputed$ID)
  hla_df_comparator_four[is.na(hla_df_comparator_four)]<-"Uncalled"

                                 
  #this is the list of imputed genes in the UKB                               
  genes<-c("A",
           "B",
           "C",
           "DRB1",
           "DQA1",
           "DQB1",
           "DPA1",
           "DPB1",
           "DRB3",
           "DRB4",
           "DRB5")
  
  #compare alleles between imputed in UKB and called with HLA-HD
  #"No match" means both errors or uncalled by UKB
  #"One match" means one match
  #"Two match" means two matches
  #"Uncalled" means HLA-HD uncalled
  ID_col<-batch_imputed$ID
  comparisons<- matrix(NA, 
                       nrow=nrow(batch_imputed), 
                       ncol=length(genes), 
                       dimnames=list(list(),genes)) %>%
    as.data.frame() %>%
    cbind(ID_col,.)
  
  for(i in 1:nrow(comparisons)){
    if((i%%1000==0)){
      print(i)
    }
    for(g in genes){
      tmp<-NA
      a1_hd<-hla_df_comparator_four[i,paste0(g,"_1")]
      a2_hd<-hla_df_comparator_four[i,paste0(g,"_2")]
      a1_imputed<-batch_imputed[i,paste0(g,"_1")]
      a2_imputed<-batch_imputed[i,paste0(g,"_2")]
      
      if(a1_hd=="Uncalled"){
        if(a1_imputed != "Unimputed"){
          if(a2_imputed == "Unimputed"){
            tmp<-"Both uncalled, one allele imputed"
          } else {
            tmp<-"Both uncalled, both alleles imputed"
          }
        } else {
          tmp<-"Both uncalled and unimputed"
        }
      } else {
        if(a1_hd == a1_imputed){
          if(a2_hd == a2_imputed){
            tmp<-"Two matches"
          } else {
            if(a2_imputed=="Unimputed"){
              if(a2_hd=="Uncalled"){
                tmp<-"One match, other uncalled and unimputed"
              } else {
                tmp<-"One match, other unimputed"
              }
            } else {
              if(a2_hd=="Uncalled"){
                tmp<-"One match, other uncalled but imputed"
              } else {
                tmp<-"One match, one mismatch imputation"
              }
            }
          }
        } else {
          if(a1_hd == a2_imputed){
            if(a2_hd == a1_imputed){
              tmp<-"Two matches"
            } else {
              if(a2_hd=="Uncalled"){
                tmp<-"One match, other uncalled but imputed"
              } else {
                tmp<-"One match, one mismatch imputation"
              }
            }
          } else {
            if(a2_hd == a1_imputed | a2_hd == a2_imputed){
              if(a2_hd == a1_imputed){
                if(a2_imputed == "Unimputed"){
                  tmp<-"One match, other unimputed"
                } else {
                  tmp<-"One match, one mismatch imputation"
                }
              } else { #i.e. a2_hd is equal to a2_imputed, but a1 are different
                tmp<-"One match, one mismatch imputation"
              }
            } else {
              if(a1_hd != "Uncalled" & a2_hd=="Uncalled"){
                if(a1_imputed=="Unimputed"){
                  tmp<-"One called, both unimputed"
                }
                if(a1_imputed != "Unimputed" & a2_imputed == "Unimputed"){
                  tmp<-"One called, one unimputed, no match"
                }
                if(a1_imputed != "Unimputed" & a2_imputed != "Unimputed"){
                  tmp<-"One called, both imputed, no match"
                }
              } else {
                if(a1_imputed=="Unimputed"){
                  tmp<-"Both called, both unimputed"
                }
                if(a1_imputed != "Unimputed" & a2_imputed == "Unimputed"){
                  tmp<-"Both called, one unimputed, no match"
                }
                if(a1_imputed != "Unimputed" & a2_imputed != "Unimputed"){
                  tmp<-"Both called, both imputed, no match"
                }
              }
            }
          }
        }
      }
      comparisons[i,g]<-tmp
    }
  }
  
  comparisons<-comparisons %>%
    mutate_if(sapply(comparisons, is.character), as.factor)
    #num_in_hla_imp
  num_in_hla_imp<-list()
  
  for(g in genes){
    basic_df<-as.data.frame(matrix(0,nrow=6,ncol=length(options)))
    colnames(basic_df)<-options
    basic_df<-bind_cols(data.frame(ancestry=c("all","afr","amr","eas","eur","sas")),
                        basic_df)
    for(a in c("afr","amr","eas","eur","sas","all")){
      if(a!="all"){
        comp_tmp<-comparisons %>% filter(ID_col %in% !!sym(a))
        hla_df_tmp<-hla_df_comparator_four %>% filter(ID %in% !!sym(a))
        for(o in options){
          ids<-comp_tmp$ID_col[which(comp_tmp[,g]==o)]
          alleles_tmp<-data.frame(alleles=c(hla_df_tmp %>% 
                                              filter(ID %in% ids) %>% 
                                              pull(paste0(g,"_1")),
                                            hla_df_tmp %>% 
                                              filter(ID %in% ids) %>% 
                                              pull(paste0(g,"_2")))) %>%
            filter(alleles %in% imputed_alleles$alleles) %>%
            nrow()
          basic_df[which(basic_df$ancestry==a),o]<-alleles_tmp
        }
      } else {
        comp_tmp<-comparisons %>% filter(ID_col %in% c(afr,amr,eas,eur,sas))
        hla_df_tmp<-hla_df_comparator_four %>% filter(ID_col %in% c(afr,amr,eas,eur,sas))
        for(o in options){
          ids<-comp_tmp$ID_col[which(comp_tmp[,g]==o)]
          alleles_tmp<-data.frame(alleles=c(hla_df_tmp %>% 
                                              filter(ID %in% ids) %>% 
                                              pull(paste0(g,"_1")),
                                            hla_df_tmp %>% 
                                              filter(ID %in% ids) %>% 
                                              pull(paste0(g,"_2")))) %>%
            filter(alleles %in% imputed_alleles$alleles) %>%
            nrow()
          basic_df[which(basic_df$ancestry==a),o]<-alleles_tmp
        }
      }
      num_in_hla_imp[[g]]<-basic_df
    }
  }

  saveRDS(num_in_hla_imp, paste0("comparisons_with_imputed/num_alleles_in_imp_ref_", i, ".RDS"))
                                 
  vroom_write(comparisons, paste0("comparisons_with_imputed/comparisons_imputed_called_folder_", folder, ".tsv.gz"))
  
}
