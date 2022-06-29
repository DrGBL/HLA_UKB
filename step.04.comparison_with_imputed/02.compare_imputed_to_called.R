library(tidyverse)
library(vroom)

#again those are folders containing the UKB participants
for(folder in 10:60){

  hla_df<-vroom(paste0("calls/hla_df_batch_", folder, ".tsv.gz"), col_types = cols(.default = "c"))

  hla_imputed_munged<-vroom("hla_imputation_df_munged.tsv", col_types = cols(.default = "c"))
  
  #now compare to current batch
  batch_imputed<-hla_imputed_munged %>% 
    filter(ID %in% hla_df$ID)
  batch_imputed[is.na(batch_imputed)]<-"Unimputed"
  batch_imputed<-lapply(batch_imputed, function(x) gsub("[0-9A-Z]*\\*99:01","Unimputed",x)) %>%
    as.data.frame()
  
  hla_df_comparator_four<-lapply(hla_df, function(x) ifelse(str_count(x,":")>1, sub(":[0-9A-Z]*$", "", x), x)) %>%
    as.data.frame() %>%
    filter(ID %in% batch_imputed$ID)
  hla_df_comparator_four[is.na(hla_df_comparator_four)]<-"Uncalled"
  
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
  
  vroom_write(comparisons, paste0("comparisons_with_imputed/comparisons_imputed_called_folder_", folder, ".tsv.gz"))
  
}
