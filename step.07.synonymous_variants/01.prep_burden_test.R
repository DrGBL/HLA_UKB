library(tidyverse)
library(vroom)

setwd("/path/to/HLA/allele/meta-analysis_results/")

n_pheno<-11

pheno<-data.frame(pheno=c("asthma",
                          "ms_demyelinating",
                          "pmr_gcc",
                          "psoriasis",
                          "auto_immune_thyroid_disorders",
                          "dm1",
                          "rheumatoid_arthritis",
                          "ulcerative_colitis",
                          "ankylosing_spondylitis",
                          "coeliac",
                          "crohn"),
                  order=c(1:n_pheno)) %>%
  arrange(pheno)


pheno_proper<-c("Asthma",
                "MS-Demyelinating",
                "PMR-GCA",
                "Psoriasis",
                "Auto-immune Thyroid Disorders",
                "DM1",
                "Rheumatoid Arthritis",
                "Ulcerative Colitis",
                "Ankylosing Spondylitis",
                "Coeliac",
                "Crohn's")
pheno_proper<-pheno_proper[pheno$order]

signif_four<-c()

for(p in pheno$pheno){
  signif_four<-vroom(paste0("four_meta_", 
                            p,
                            "1.txt.gz")) %>%
    filter(Pvalue<5e-8/n_pheno) %>% 
    mutate(phenotype=p) %>%
    rename(Effect_four=Effect,
           StdErr_four=StdErr,
           Pvalue_four=Pvalue,
           HetPVal_four=HetPVal,
           Num_Cases_four=Num_Cases,
           Cases_Ref_four=Cases_Ref,
           Cases_Het_four=Cases_Het,
           Cases_Alt_four=Cases_Alt,
           Num_Controls_four=Num_Controls,
           Controls_Ref_four=Controls_Ref,
           Controls_Het_four=Controls_Het,
           Controls_Alt_four=Controls_Alt) %>%
    dplyr::select(c(MarkerName,
                    Effect_four,
                    StdErr_four,
                    Pvalue_four,
                    HetPVal_four,
                    Num_Cases_four,
                    Cases_Ref_four,
                    Cases_Het_four,
                    Cases_Alt_four,
                    Num_Controls_four,
                    Controls_Ref_four,
                    Controls_Het_four,
                    Controls_Alt_four,
                    phenotype)) %>%
    rbind(signif_four,.)
}



#now check 6 digits
signif_six<-c()

for(p in pheno$pheno){
  signif_six<-vroom(paste0("six_meta_", 
                           p,
                           "1.txt.gz")) %>%
    mutate(phenotype=p) %>%
    mutate(MarkerName_four = str_extract(MarkerName, "[A-Z0-9]*\\*[0-9]*:[0-9]*")) %>%
    rename(MarkerName_six=MarkerName) %>%
    rename(MarkerName=MarkerName_four) %>%
    merge(.,signif_four) %>%
    rbind(signif_six,.)
}

write_tsv(signif_six, "/scratch/richards/guillaume.butler-laporte/hla_ukb_wes/regenie_burden_outputs/four_to_six_digits_alelles.tsv.gz", col_names = TRUE)

#now prep the burden tests
anno_file<-signif_six %>%
  filter(Pvalue>=5e-8/n_pheno) %>%
  dplyr::select(c(MarkerName_six, MarkerName)) %>%
  mutate(type=rep("all", nrow(.))) %>%
  distinct()

#check which ones have only two 6 digit allele in the file above (so only one to compare to the most significant one), and print for future use
anno_file_single_allele<-anno_file %>%
  group_by(MarkerName) %>%
  filter(n()<3) %>%
  ungroup()

marker_unique<-anno_file %>% arrange(MarkerName) %>% pull(MarkerName) %>% unique()
set_list<-data.frame(marker=marker_unique,
                     chr=rep(6,length(marker_unique)),
                     pos=c(1:length(marker_unique)),
                     alleles=NA)
for(m in 1:nrow(set_list)){
  set_list$alleles[m]<-anno_file %>%
    filter(MarkerName==set_list$marker[m]) %>%
    pull(MarkerName_six) %>%
    paste(., collapse=",")
}

write_tsv(anno_file, "anno_file.tsv.gz", col_names = FALSE)
write_tsv(set_list, "set_list.tsv.gz", col_names = FALSE)
write_tsv(anno_file_single_allele, "anno_file_single_allele.tsv.gz", col_names = TRUE)
