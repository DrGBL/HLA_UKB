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

#build new data frame
hla_imputed_munged<-data.frame(ID=hla_imputed$ID,
                               A_1=rep(NA,nrow(hla_imputed)),
                               A_2=rep(NA,nrow(hla_imputed)),
                               B_1=rep(NA,nrow(hla_imputed)),
                               B_2=rep(NA,nrow(hla_imputed)),
                               C_1=rep(NA,nrow(hla_imputed)),
                               C_2=rep(NA,nrow(hla_imputed)),
                               DRB1_1=rep(NA,nrow(hla_imputed)),
                               DRB1_2=rep(NA,nrow(hla_imputed)),
                               DQA1_1=rep(NA,nrow(hla_imputed)),
                               DQA1_2=rep(NA,nrow(hla_imputed)),
                               DQB1_1=rep(NA,nrow(hla_imputed)),
                               DQB1_2=rep(NA,nrow(hla_imputed)),
                               DPA1_1=rep(NA,nrow(hla_imputed)),
                               DPA1_2=rep(NA,nrow(hla_imputed)),
                               DPB1_1=rep(NA,nrow(hla_imputed)),
                               DPB1_2=rep(NA,nrow(hla_imputed)),
                               DRB3_1=rep(NA,nrow(hla_imputed)),
                               DRB3_2=rep(NA,nrow(hla_imputed)),
                               DRB4_1=rep(NA,nrow(hla_imputed)),
                               DRB4_2=rep(NA,nrow(hla_imputed)),
                               DRB5_1=rep(NA,nrow(hla_imputed)),
                               DRB5_2=rep(NA,nrow(hla_imputed)))

#this assigns hard calls
for(s in 1:nrow(hla_imputed)){
  if((s%%1000==0)){
    print(paste0(s,"/", nrow(hla_imputed)))
  }
  
  for(a in 2:ncol(hla_imputed)){
    if(hla_imputed[s,a] >= 0.8){
      gene_base <- sub("\\*.*", "", colnames(hla_imputed)[a])
      if(hla_imputed[s,a] >= 1.6){
        hla_imputed_munged[s,paste0(gene_base, "_1")]<-colnames(hla_imputed)[a]
        hla_imputed_munged[s,paste0(gene_base, "_2")]<-colnames(hla_imputed)[a]
      } else {
        if(is.na(hla_imputed_munged[s,paste0(gene_base, "_1")])){
          hla_imputed_munged[s,paste0(gene_base, "_1")]<-colnames(hla_imputed)[a]
        } else {
          hla_imputed_munged[s,paste0(gene_base, "_2")]<-colnames(hla_imputed)[a]
        }
      }
    }
  }
}
  
write_tsv(hla_imputed_munged, paste0(path_out, "hla_imputation_df_munged.tsv"))
