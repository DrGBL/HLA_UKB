library(tidyverse)
library(vroom)

n_pheno<-11

anno_file_single_allele<-vroom("anno_file_single_allele.tsv.gz") #from step 01

path_meta_burden<-"/path/to/meta-analysis_burden_test_results/"
path_meta_six_digits<-"/path/to/meta-analysis_6_digit_hla_alleles/"

pheno_df<-data.frame(phenotype=c("asthma",
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
                     Phenotype=c("Asthma",
                                 "MS-Demyelinating",
                                 "PMR-GCA",
                                 "Psoriasis",
                                 "Auto-immune Thyroid Disorders",
                                 "DM1",
                                 "Rheumatoid Arthritis",
                                 "Ulcerative Colitis",
                                 "Ankylosing Spondylitis",
                                 "Coeliac",
                                 "Crohn's"))

welsh_gwas<-function(beta1, beta2, se1, se2, n1, n2){
  t<-(beta1-beta2)/(sqrt(se1^2+se2^2))
  
  v<-(( (se1^2)/n1 +  (se2^2)/n2 )^2) / ( (se1^4)/( (n1^2)*(n1-1)) + (se2^4)/( (n2^2)*(n2-1)) )
  
  res<-2*pt(q=abs(t), df=v, lower.tail=FALSE)
  
  return(res)
}

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

burden_res<-c()

coding_genes<-c("A", "B", "C", "E", "F", "G",  
                "DMA", "DMB", "DOA", "DOB",
                "DPA1", "DPB1", "DQA1", "DQB1",
                "DRA", "DRB1", "DRB3", "DRB4", "DRB5")

for(p in pheno$pheno){
  burden_res<-vroom(paste0(path_meta_burden, "burden_meta_", 
                           p,
                           "1.txt.gz")) %>%
    mutate(Effect=ifelse(Allele1=="ref" | Allele1=="REF", -Effect, Effect)) %>%
    mutate(gene=str_extract(MarkerName, "[A-Z0-9]*")) %>%
    mutate(coding=ifelse(gene %in% coding_genes, "Coding", "Non-coding")) %>%
    mutate(phenotype=p) %>%
    mutate(MarkerName=str_extract(MarkerName, "[A-Z0-9]*\\*[0-9]*:[0-9]*")) %>%
    rename(Effect_burden=Effect,
           StdErr_burden=StdErr,
           Pvalue_burden=Pvalue,
           Num_Cases_burden=Num_Cases,
           Cases_Ref_burden=Cases_Ref,
           Cases_Het_burden=Cases_Het,
           Cases_Alt_burden=Cases_Alt,
           Num_Controls_burden=Num_Controls,
           Controls_Ref_burden=Controls_Ref,
           Controls_Het_burden=Controls_Het,
           Controls_Alt_burden=Controls_Alt) %>%
    filter(Allele1=="mask1.0.999999") %>%
    dplyr::select(c(MarkerName,
                    Effect_burden,
                    StdErr_burden,
                    Pvalue_burden,
                    Num_Cases_burden,
                    Cases_Ref_burden,
                    Cases_Het_burden,
                    Cases_Alt_burden,
                    Num_Controls_burden,
                    Controls_Ref_burden,
                    Controls_Het_burden,
                    Controls_Alt_burden,
                    coding,
                    phenotype)) %>%
    mutate(AF_burden=(Cases_Het_burden+2*Cases_Alt_burden+Controls_Het_burden+2*Controls_Alt_burden)/(2*Num_Cases_burden+2*Num_Controls_burden)) %>%
    rbind(burden_res,.)
}

#first check those with only two significant 6 digit alleles
six_digit_sig<-c()
six_digit_non_sig<-c()
for(p in pheno$pheno){
  six_digit_sig<-vroom(paste0(path_meta_six_digits, "six_meta_", 
                              p,
                              "1.txt.gz")) %>%
    mutate(Effect=ifelse(Allele1=="c" | Allele1=="C", -Effect, Effect)) %>%
    mutate(phenotype=p) %>%
    rename(Effect_six=Effect,
           StdErr_six=StdErr,
           Pvalue_six=Pvalue,
           Num_Cases_six=Num_Cases,
           Cases_Ref_six=Cases_Ref,
           Cases_Het_six=Cases_Het,
           Cases_Alt_six=Cases_Alt,
           Num_Controls_six=Num_Controls,
           Controls_Ref_six=Controls_Ref,
           Controls_Het_six=Controls_Het,
           Controls_Alt_six=Controls_Alt) %>%
    dplyr::select(c(MarkerName,
                    Effect_six,
                    StdErr_six,
                    Pvalue_six,
                    Num_Cases_six,
                    Cases_Ref_six,
                    Cases_Het_six,
                    Cases_Alt_six,
                    Num_Controls_six,
                    Controls_Ref_six,
                    Controls_Het_six,
                    Controls_Alt_six,
                    phenotype)) %>%
    mutate(AF=(Cases_Het_six+2*Cases_Alt_six+Controls_Het_six+2*Controls_Alt_six)/(2*Num_Cases_six+2*Num_Controls_six)) %>%
    filter(Pvalue_six<5e-8/n_pheno) %>%
    mutate(four_digit=str_extract(MarkerName, "^[A-Z0-9]*\\*[0-9]*:[0-9]*")) %>%
    filter(MarkerName != four_digit) %>%
    rbind(six_digit_sig,.)
  
  six_digit_non_sig<-vroom(paste0(path_meta_six_digits, "six_meta_", 
                              p,
                              "1.txt.gz")) %>%
    mutate(Effect=ifelse(Allele1=="c" | Allele1=="C", -Effect, Effect)) %>%
    mutate(phenotype=p) %>%
    rename(Effect_six=Effect,
           StdErr_six=StdErr,
           Pvalue_six=Pvalue,
           Num_Cases_six=Num_Cases,
           Cases_Ref_six=Cases_Ref,
           Cases_Het_six=Cases_Het,
           Cases_Alt_six=Cases_Alt,
           Num_Controls_six=Num_Controls,
           Controls_Ref_six=Controls_Ref,
           Controls_Het_six=Controls_Het,
           Controls_Alt_six=Controls_Alt) %>%
    dplyr::select(c(MarkerName,
                    Effect_six,
                    StdErr_six,
                    Pvalue_six,
                    Num_Cases_six,
                    Cases_Ref_six,
                    Cases_Het_six,
                    Cases_Alt_six,
                    Num_Controls_six,
                    Controls_Ref_six,
                    Controls_Het_six,
                    Controls_Alt_six,
                    phenotype)) %>%
    mutate(AF=(Cases_Het_six+2*Cases_Alt_six+Controls_Het_six+2*Controls_Alt_six)/(2*Num_Cases_six+2*Num_Controls_six)) %>%
    filter(AF>0.01 & AF<0.99) %>%
    mutate(four_digit=str_extract(MarkerName, "^[A-Z0-9]*\\*[0-9]*:[0-9]*")) %>%
    filter(MarkerName != four_digit) %>%
    filter(four_digit %in% (six_digit_sig %>% filter(phenotype==p) %>% pull(four_digit))) %>%
    filter(Pvalue_six>5e-8/n_pheno) %>%
    rbind(six_digit_non_sig,.)
}

#note that using filter(n()>1) gave same result as filter(n()==2)
six_digit_sig_two_or_more<-six_digit_sig %>%
  group_by(phenotype, four_digit) %>%
  filter(n()==2) %>%
  arrange(phenotype,MarkerName) %>%
  ungroup()

six_digit_sig_het_test<-data.frame(allele_four=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   phenotype=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   allele_six_1=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   beta1=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   se1=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   pval1=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   AF1=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   N1=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   allele_six_2=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   beta2=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   se2=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   pval2=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   AF2=rep(NA,nrow(six_digit_sig_two_or_more)/2),
                                   N2=rep(NA,nrow(six_digit_sig_two_or_more)/2))

for(i in 1:(nrow(six_digit_sig_two_or_more)/2)){
  six_digit_sig_het_test$allele_four[i]<-six_digit_sig_two_or_more$four_digit[2*i-1]
  six_digit_sig_het_test$phenotype[i]<-six_digit_sig_two_or_more$phenotype[2*i-1]
  
  six_digit_sig_het_test$allele_six_1[i]<-six_digit_sig_two_or_more$MarkerName[2*i-1]
  six_digit_sig_het_test$allele_six_2[i]<-six_digit_sig_two_or_more$MarkerName[2*i]
  
  six_digit_sig_het_test$beta1[i]<-six_digit_sig_two_or_more$Effect_six[2*i-1]
  six_digit_sig_het_test$beta2[i]<-six_digit_sig_two_or_more$Effect_six[2*i]
  
  six_digit_sig_het_test$se1[i]<-six_digit_sig_two_or_more$StdErr_six[2*i-1]
  six_digit_sig_het_test$se2[i]<-six_digit_sig_two_or_more$StdErr_six[2*i]
  
  six_digit_sig_het_test$pval1[i]<-six_digit_sig_two_or_more$Pvalue_six[2*i-1]
  six_digit_sig_het_test$pval2[i]<-six_digit_sig_two_or_more$Pvalue_six[2*i]
  
  six_digit_sig_het_test$AF1[i]<-six_digit_sig_two_or_more$AF[2*i-1]
  six_digit_sig_het_test$AF2[i]<-six_digit_sig_two_or_more$AF[2*i]
  
  six_digit_sig_het_test$N1[i]<-six_digit_sig_two_or_more$Num_Cases_six[2*i-1]+six_digit_sig_two_or_more$Num_Controls_six[2*i-1]
  six_digit_sig_het_test$N2[i]<-six_digit_sig_two_or_more$Num_Cases_six[2*i]+six_digit_sig_two_or_more$Num_Controls_six[2*i]
}

six_digit_sig_het_test<-six_digit_sig_het_test %>%
  mutate(het_p=welsh_gwas(beta1=beta1, 
                          beta2=beta2, 
                          se1=se1, 
                          se2=se2, 
                          n1=N1, 
                          n2=N2))

#now check significant alleles against non-significant ones
six_digit_non_sig_het_test<-six_digit_sig %>%
  dplyr::select(c(phenotype, 
                  MarkerName,
                  Effect_six,
                  StdErr_six,
                  Pvalue_six, 
                  four_digit,
                  Num_Cases_six,
                  Num_Controls_six,
                  AF)) %>%
  mutate(N1=Num_Cases_six+Num_Controls_six) %>%
  rename(allele_six_1=MarkerName,
         beta1=Effect_six,
         se1=StdErr_six,
         pval1=Pvalue_six,
         AF1=AF) %>%
  dplyr::select(-c(Num_Cases_six,
                   Num_Controls_six)) %>%
  merge(., six_digit_non_sig %>%
          dplyr::select(c(phenotype, 
                          MarkerName,
                          Effect_six,
                          StdErr_six,
                          Pvalue_six, 
                          four_digit,
                          Num_Cases_six,
                          Num_Controls_six,
                          AF)) %>%
          mutate(N2=Num_Cases_six+Num_Controls_six) %>%
          rename(allele_six_2=MarkerName,
                 beta2=Effect_six,
                 se2=StdErr_six,
                 pval2=Pvalue_six,
                 AF2=AF) %>%
          dplyr::select(-c(Num_Cases_six,
                           Num_Controls_six))) %>%
  mutate(het_p=welsh_gwas(beta1=beta1, 
                          beta2=beta2, 
                          se1=se1, 
                          se2=se2, 
                          n1=N1, 
                          n2=N2))
  

#now test the burden test effects against those significant alleles
total_df<-burden_res %>%
  rename(four_digit=MarkerName) %>%
  merge(six_digit_sig %>% rename(AF_six=AF),.) %>%
  mutate(het_p=welsh_gwas(beta1=Effect_six, 
                          beta2=Effect_burden, 
                          se1=StdErr_six, 
                          se2=StdErr_burden, 
                          n1=Num_Cases_six+Num_Controls_six, 
                          n2=Num_Cases_burden+Num_Controls_burden))

total_df_multiple_alleles<-total_df %>%
  filter(!(four_digit %in% anno_file_single_allele$MarkerName))

#now build the final data drame

tmp_all_six<-six_digit_sig_het_test %>%
  bind_rows(.,six_digit_non_sig_het_test) %>%
  dplyr::select(c(phenotype,
                  allele_six_1,
                  beta1,
                  se1,
                  pval1,
                  AF1,
                  allele_six_2,
                  beta2,
                  se2,
                  pval2,
                  AF2,
                  het_p)) %>%
  rename("Six-digit allele"=allele_six_1) %>%
  rename("Six-digit allele beta"=beta1) %>%
  rename("Six-digit allele s.e."=se1) %>%
  rename("Six-digit allele p-value"=pval1) %>%
  rename("Six-digit allele AF"=AF1) %>%
  rename("Comparator"=allele_six_2) %>%
  rename("Comparator beta"=beta2) %>%
  rename("Comparator s.e."=se2) %>%
  rename("Comparator p-value"=pval2) %>%
  rename("Comparator allele AF"=AF2) %>%
  rename("Heterogeneity p-value"=het_p)

tmp_burden<-total_df_multiple_alleles %>%
  mutate(Comparator=rep("Burden", nrow(.))) %>%
  dplyr::select(c(phenotype,
                  MarkerName,
                  Effect_six,
                  StdErr_six,
                  Pvalue_six,
                  AF_six,
                  Comparator,
                  Effect_burden,
                  StdErr_burden,
                  Pvalue_burden,
                  AF_burden,
                  het_p)) %>%
  rename("Six-digit allele"=MarkerName) %>%
  rename("Six-digit allele beta"=Effect_six) %>%
  rename("Six-digit allele s.e."=StdErr_six) %>%
  rename("Six-digit allele p-value"=Pvalue_six) %>%
  rename("Six-digit allele AF"=AF_six) %>%
  rename("Comparator beta"=Effect_burden) %>%
  rename("Comparator s.e."=StdErr_burden) %>%
  rename("Comparator p-value"=Pvalue_burden) %>%
  rename("Comparator allele AF"=AF_burden) %>%
  rename("Heterogeneity p-value"=het_p)

heterogeneity_df<-bind_rows(tmp_all_six,tmp_burden) %>%
  mutate(Significant=`Heterogeneity p-value` < 0.05/nrow(.)) %>%
  merge(pheno_df,.) %>%
  dplyr::select(-phenotype) %>%
  mutate(Direction=ifelse(`Six-digit allele beta`*`Comparator beta`>0, "Same", "Opposite")) %>%
  mutate(Direction=ifelse(`Comparator p-value`>5e-8/11, "Nil effect", Direction))

#plot these results
plot_het<-heterogeneity_df %>%
  mutate(Comparator=ifelse(Comparator=="Burden", "Dummy allele\nburden test", "Other\n3-field allele")) %>%
  mutate(`Effect category`=ifelse(Significant, "Heterogeneous effect", "Non-heterogeneous\neffect")) %>%
  mutate(`Effect category`=ifelse(Significant & Direction=="Same", "Heterogeneus effect\nin same directions", `Effect category`)) %>%
  mutate(`Effect category`=ifelse(Significant & Direction=="Opposite", "Heterogeneus effect\nin opposite directions", `Effect category`)) %>%
  mutate(`Effect category`=ifelse(Significant & Direction=="Nil effect", "Heterogeneus effect,\none allele with null effect", `Effect category`)) %>%
  ggplot(data=., 
         aes(x=Comparator, fill=`Effect category`)) +
  geom_bar() +
  xlab("")+
  ylab("Number of pairs of heterogeneity tests") +
  scale_fill_manual(values = c("#d7191c", "#abd9e9", "#fdae61", "#2c7bb6") )+
  theme_bw()

ggsave("plot_het.pdf", plot_het, height=7, width=7)

write_tsv(file="four_vs_six_digits_burden_all_res.tsv.gz",
            x=heterogeneity_df)

