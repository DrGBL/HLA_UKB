# The UK Biobank provides HLA imputation from HLA:IMP*2 in data field 22182, with sample row in ressource 1520, and header in Resource 2182.
# Note that alleles that end in "99:01" are to be considered not imputed, refer to the documentation for more details
# Imputed HLA are provided as an allele dosage from 0 to 2, for each allele, and it is recommended that scores below 0.8 be set to 0, if hard calls are needed.
# The following code assumes that the user has already taken these data and made them into a dataframe where the first columns are the participant FID and IID
# and each following columns is the header from ressource 2182.
# Hence, there are as many rows as participants, and each column (other than the first FID and IID columns) is a dosage of the corresponding imputed allele.

# this code munges this data frame to assign hard calls based on dosage, which will be used to compare to sequencing calls in the next steps.

library(tidyverse)
library(vroom)

# path to the imputed file described in the comments above
path_to_imputed_ukb_file<-"/path/to/imputed/ukb/file"

#path to where the munged data should be saved
path_out<-"/path/to/output/folder/"

#now start munging
hla_imputed<-vroom(path_to_imputed_ukb_file) %>%
  rename(ID=IID) %>%
  dplyr::select(-FID)

#modify the colnames
cols_hla<-colnames(hla_imputed)[-1] %>% data.frame(names=.) %>%
  separate(col=names, into=c("gene", "allele"), sep="_") %>%
  mutate(allele=ifelse(nchar(allele)==3,paste0("0", allele),allele)) %>%
  mutate(allele=paste0(substr(allele,1,2),":",substr(allele,3,4)))

colnames(hla_imputed)[-1]<-paste0(cols_hla$gene,"*",cols_hla$allele)

#munge
threshold<-0.8   #if this is above 2/3, then you could get in situation where there will be more than 3 allele calls
hla_imputed_munged<-data.frame(ID=NA)
for(group in c("A", "B", "C", "DRB1", "DRB3", "DRB4",
               "DRB5", "DPA1", "DPB1", "DQA1", "DQB1")){
  cols_to_keep<-colnames(hla_imputed)[which(startsWith(colnames(hla_imputed), group))]
  
  alleles<-c()
  
  for(g in cols_to_keep){
    print(g)
    
    tmp<-hla_imputed %>% 
      dplyr::select(c(ID,all_of(g))) %>%
      filter(!!sym(g)>threshold) %>%
      mutate(n=ifelse(!!sym(g)>2*threshold,2,1)) %>%
      mutate(c=g) %>%
      dplyr::select(c(ID,n,c))
    
    tmp2<-rbind(tmp,
                tmp %>% 
                  filter(n == 2)) %>%
      dplyr::select(-n) %>%
      rename(!!sym(group):=c)
    
    alleles<-alleles %>%
      rbind(.,tmp2)
    
  }
  
  alleles<-alleles %>%
    group_by(ID) %>%
    mutate(n=1:n()) %>%
    ungroup()
  
  alleles_p1<-alleles %>%
    filter(n==1) %>%
    rename(!!paste0(group,"_1"):=!!sym(group)) %>%
    dplyr::select(-n)
  
  alleles_p2<-alleles %>%
    filter(n==2) %>%
    rename(!!paste0(group,"_2"):=!!sym(group)) %>%
    dplyr::select(-n)
  
  hla_imputed_munged<-full_join(hla_imputed_munged,
                   full_join(alleles_p1,
                             alleles_p2))
  
}

hla_imputed_munged<-hla_imputed_munged %>%
  filter(!is.na(ID))
 
write_tsv(hla_imputed_munged, paste0(path_out, "hla_imputation_df_munged.tsv"))
