library(tidyverse)
library(vroom)
library(genpwr)
library(ggcorrplot)
library(ggpubr)
library(cowplot)

setwd("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code")

path_to_hla_spread_csv<-"path/to/HLA_SPREAD.csv"  #this file can be found on the git (HLA_SPREAD.zip)
path_wes_meta_analysis<-"path/to/meta-analysis/from/wes/alleles/"
path_imputed_meta_analysis<-"path/to/meta-analysis/from/imputed/alleles/"

n_pheno<-11

#first need allele frequency files
#these are files with two columns with headers
#first column is the allele ID with header "MarkerName"
#second column is the allele frequency with header "MAF" (to be specific, if the allele's AF is > 50%, then MAF for this allele is 1-AF).
hla_four_freq<-vroom(path_four_digit_allele_freq)

hla_six_freq<-vroom(path_six_digit_allele_freq)

#known alleles (hla spread)
hla_spread_tmp <- read_csv(path_to_hla_spread_csv, col_types = cols(.default = "c")) %>%
  filter(category=="Diseases") %>%
  filter(inner_group!="HAPLOTYPE" & outer_group!="HAPLOTYPE") %>%
  mutate(pos=ifelse(str_detect(third_group,"HLA-[A-Z0-9]+:[0-9]+"),
                    nchar(str_extract(third_group,"HLA-[A-Z0-9]+:")),
                    NA)) %>%
  mutate(third_group=ifelse(!is.na(pos),
                            paste0(substr(third_group,1,pos-1),
                                   "*",
                                   substr(third_group,pos+1,nchar(third_group))),
                            third_group)) %>%
  filter(str_detect(third_group,"HLA-[A-Z0-9]*\\*[0-9]+:[0-9]+")) %>%
  mutate(tmp=str_extract(third_group,"[A-Z0-9]*\\*[0-9]+:[0-9]+")) %>%
  filter(updated_class=="positive") %>%
  mutate(phenotype=rep(NA, nrow(.))) %>%
  mutate(phenotype=ifelse(disease=="ankylosing spondylitis", "ankylosing_spondylitis", phenotype)) %>%
  mutate(phenotype=ifelse(disease=="asthma" |
                            disease=="bronchial asthma", 
                          "asthma", 
                          phenotype)) %>%
  mutate(phenotype=ifelse(disease=="thyroiditis" | 
                            disease=="graves' disease" | 
                            disease=="thyroid diseases" |
                            disease=="hashimoto thyroiditis" |
                            disease=="goiter" |
                            disease=="thyroid associated ophthalmopathy",
                          "auto_immune_thyroid_disorders", 
                          phenotype)) %>%
  mutate(phenotype=ifelse(disease=="celiac disease", "coeliac", phenotype)) %>%
  mutate(phenotype=ifelse(disease=="crohns disease", "crohn", phenotype)) %>%
  mutate(phenotype=ifelse(disease=="type 1 diabetes" |
                            disease=="insulin-dependent diabetes mellitus" |
                            disease=="type 1 diabetes mellitus", 
                          "dm1", 
                          phenotype)) %>%
  mutate(phenotype=ifelse(disease=="multiple sclerosis", "ms_demyelinating", phenotype)) %>%
  mutate(phenotype=ifelse(disease=="giant cell arteritis", "pmr_gcc", phenotype)) %>%
  mutate(phenotype=ifelse(disease=="psoriasis" |
                            disease=="psoriatic arthritis",
                          "psoriasis",
                          phenotype)) %>%
  mutate(phenotype=ifelse(disease=="rheumatoid arthritis" |
                            disease=="rheumatoid nodules", 
                          "rheumatoid_arthritis", 
                          phenotype)) %>%
  mutate(phenotype=ifelse(disease=="ulcerative colitis",
                          "ulcerative_colitis",
                          phenotype)) %>%
  filter(!is.na(phenotype)) %>%
  dplyr::select(c(third_group, phenotype)) %>%
  rename(allele=third_group) %>%
  mutate(allele=gsub("HLA-","",allele)) %>%
  distinct() %>%
  mutate(source="HLA-SPREAD")

#known alleles (imputed allele imputation)
imp_alleles<-c()
for(p in c("ankylosing_spondylitis", "asthma", "auto_immune_thyroid_disorders",
           "coeliac", "coeliac", "dm1", "ms_demyelinating", "pmr_gcc",
           "psoriasis", "rheumatoid_arthritis", "ulcerative_colitis")){
  imp_alleles<-vroom(paste0(path_imputed_meta_analysis, "/four_meta_",
                            p, "1.txt.gz")) %>%
    rename(allele=MarkerName) %>%
    mutate(phenotype=p) %>%
    mutate(source="Imputation") %>%
    bind_rows(imp_alleles,.)
}

hla_spread<-bind_rows(hla_spread_tmp, imp_alleles %>% filter(Pvalue < (5e-8)/11) %>% dplyr::select(c(allele, phenotype, source))) %>% distinct()

#phenotype and gene names
pheno<-data.frame(pheno=c("asthma",
                          "ms_demyelinating",
                          "myasthenia_gravis",
                          "pernicious_anemia",
                          "pmr_gcc",
                          "psoriasis",
                          "sjogren",
                          "auto_immune_thyroid_disorders",
                          "dm1",
                          "lupus",
                          "rheumatoid_arthritis",
                          "ulcerative_colitis",
                          "ankylosing_spondylitis",
                          "coeliac",
                          "crohn"),
                  order=c(1:15)) %>%
  arrange(pheno)


pheno_proper<-c("Asthma",
                "MS-Demyelinating",
                "Myasthenia gravis",
                "Pernicious Anemia",
                "PMR-GCA",
                "Psoriasis",
                "Sjogren",
                "Auto-immune Thyroid Disorders",
                "DM1",
                "Lupus",
                "Rheumatoid Arthritis",
                "Ulcerative Colitis",
                "Ankylosing Spondylitis",
                "Coeliac",
                "Crohn's")
pheno_proper<-pheno_proper[pheno$order]

hla_genes<-c("A", "B", "C", "DMA", "DMB", 
             "DOA", "DOB", "DPA1", "DPA2", "DPB1", 
             "DQA1", "DQB1", "DRA", "DRB1", "DRB2",
             "DRB3", "DRB4", "DRB5", "DRB6", "DRB7", 
             "DRB8", "DRB9", "E", "F", "G",
             "H", "J", "K", "L", "Y")

hla_genes_coding<-c("A", "B", "C", "DMA", "DMB", 
                    "DOA", "DOB", "DPA1", "DPB1", 
                    "DQA1", "DQB1", "DRA", "DRB1",
                    "DRB3", "DRB4", "DRB5", 
                    "E", "F", "G")

classI<-c("A", "B", "C", "E", "F", "G",
          "H", "J", "K", "L", "Y")
classII<-c("DMA", "DMB", 
           "DOA", "DOB", "DPA1", "DPA2", "DPB1",
           "DQA1", "DQB1", "DRA", "DRB1", "DRB2",
           "DRB3", "DRB4", "DRB5", "DRB6", "DRB7", 
           "DRB8", "DRB9")

#load results and prep them
hla_four<-c()
hla_six<-c()
for(p in pheno$pheno){
  hla_four<-vroom(paste0(path_wes_meta_analysis, "four_meta_",
                         p, "1.txt.gz")) %>%
    mutate(original_allele=gsub(".*\\*","",MarkerName)) %>%
    mutate(Phenotype=p) %>%
    mutate(Gene=str_extract(MarkerName,"[A-Z0-9]*")) %>%
    mutate(Class=ifelse(Gene %in% classI, "Class I", "Class II")) %>%
    mutate(Novel=ifelse(MarkerName %in% subset(hla_spread, hla_spread$phenotype==p & hla_spread$source == "HLA-SPREAD")$allele,
                        "Previously reported in\na Pubmed abstract",
                        ifelse(MarkerName %in% subset(hla_spread, hla_spread$phenotype==p & hla_spread$source == "Imputation")$allele,
                               "Also observed using\nthe UKB imputed alleles",
                               "Novel association"))) %>%
    bind_rows(hla_four,.)
  
  hla_six<-vroom(paste0(path_wes_meta_analysis, "six_meta_",
                        p, "1.txt.gz")) %>%
    mutate(original_allele=gsub(".*\\*","",MarkerName)) %>%
    mutate(Phenotype=p) %>%
    mutate(Gene=str_extract(MarkerName,"[A-Z0-9]*")) %>%
    mutate(Class=ifelse(Gene %in% classI, "Class I", "Class II")) %>%
    bind_rows(hla_six,.)
}

hla_four<-hla_four %>% mutate(Novel=factor(Novel, levels=c("Novel association","Previously reported in\na Pubmed abstract","Also observed using\nthe UKB imputed alleles")))

#power curves
sample_sizes<-data.frame(pheno=pheno$pheno,
                         cases=NA_integer_,
                         controls=NA_integer_)
for(i in 1:15){
  sample_sizes[which(sample_sizes$pheno==pheno$pheno[i]),]$cases<-hla_four %>%
    filter(Phenotype==pheno$pheno[i]) %>%
    summarize(max(Num_Cases)) %>% as.numeric()
  sample_sizes[which(sample_sizes$pheno==pheno$pheno[i]),]$controls<-hla_four %>%
    filter(Phenotype==pheno$pheno[i]) %>%
    summarize(max(Num_Controls)) %>% as.numeric()
}
sample_sizes$n_total<-sample_sizes$cases+sample_sizes$controls

mafs<-10^seq(-3.5,-0.40,0.02)


or<-c()
for(i in 1:15){
  print(i)
    or[[paste0(sample_sizes$pheno[i])]] <- genpwr.calc(calc = "es", model = "logistic", ge.interaction = NULL,
                                                       N=sample_sizes$n_total[i], Case.Rate=sample_sizes$cases[i]/sample_sizes$n_total[i],
                                                       MAF=mafs, Power=0.8, Alpha=(5e-8)/n_pheno,
                                                       True.Model="Additive", Test.Model="Additive")
}
 

# if needed
# saveRDS(or,"power_curves_list.RDS")
$ or<-readRDS("power_curves_list.RDS")

or_df<-c()
for(i in 1:length(or)){
  or_df<-or[[i]] %>%
    mutate(Phenotype=names(or)[i]) %>%
    bind_rows(or_df,.)
}

or_df_inv<-or_df %>%
  mutate(MAF=1-MAF)

or_df<-or_df %>%
  #bind_rows(.,or_df_inv) %>%
  mutate(Phenotype=factor(Phenotype, levels=pheno$pheno, labels=pheno_proper)) %>%
  mutate(Phenotype=droplevels(Phenotype))

#plots
all_alleles_six_digit<-hla_six %>%
  mutate(Gene=factor(Gene, levels=hla_genes)) %>%
  mutate(Gene=droplevels(Gene)) %>%
  mutate(Phenotype=factor(Phenotype, levels=pheno$pheno, labels=pheno_proper)) %>%
  mutate(Phenotype=droplevels(Phenotype)) %>%
  filter(Pvalue<5e-8/n_pheno) %>%
  mutate(Gene=droplevels(Gene)) %>%
  merge(.,hla_six_freq) %>%
  mutate(Effect=ifelse(MAF>0.5,-Effect,Effect)) %>%
  mutate(MAF=ifelse(MAF>0.5, 1-MAF, MAF)) %>%
  ggplot(aes(y = exp(Effect),
             x = MAF)) + 
  geom_point(aes(size = -log10(Pvalue),
                 colour=Class))+
  geom_line(data=or_df, aes(y=`OR_at_Alpha_4.54545454545455e-09`, x=MAF), linetype="dotted")+
  geom_line(data=or_df, aes(y=exp(-log(`OR_at_Alpha_4.54545454545455e-09`)), x=MAF), linetype="dotted")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~Phenotype,scales = "free")+
  theme_bw() +
  ylab("Odds Ratios") +
  xlab("Minor Allele Frequencies")+
  geom_hline(yintercept=1, alpha=0.5, linetype = 'dashed') +
  annotation_logticks(sides="b") +
  guides(size=guide_legend(title=expression(-log[10](p-value))))

all_alleles_four_digit<-hla_four %>%
  mutate(Gene=factor(Gene, levels=hla_genes)) %>%
  mutate(Gene=droplevels(Gene)) %>%
  mutate(Phenotype=factor(Phenotype, levels=pheno$pheno, labels=pheno_proper)) %>%
  mutate(Phenotype=droplevels(Phenotype)) %>%
  filter(Pvalue<5e-8/n_pheno) %>%
  mutate(Gene=droplevels(Gene)) %>%
  merge(.,hla_four_freq) %>%
  mutate(Effect=ifelse(MAF>0.5,-Effect,Effect)) %>%
  mutate(MAF=ifelse(MAF>0.5, 1-MAF, MAF)) %>%
  ggplot(aes(y = exp(Effect),
             x = MAF)) + 
  geom_point(aes(size = -log10(Pvalue),
                 colour=Class,
                 shape=Novel))+
  geom_line(data=or_df, aes(y=`OR_at_Alpha_4.54545454545455e-09`, x=MAF), linetype="dotted")+
  geom_line(data=or_df, aes(y=exp(-log(`OR_at_Alpha_4.54545454545455e-09`)), x=MAF), linetype="dotted")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~Phenotype,scales = "free")+
  theme_bw() +
  ylab("Odds Ratios") +
  xlab("Minor Allele Frequencies")+
  geom_hline(yintercept=1, alpha=0.5, linetype = 'dashed') +
  annotation_logticks(sides="b") +
  guides(size=guide_legend(title=expression(-log[10](p-value))))


hla_four %>%
  filter(Pvalue<5e-8/11) %>%
  write_tsv("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code/novel_4_digit.tsv")
  
ggsave("association_studies_results_plots/all_alleles_four_digit.pdf", all_alleles_four_digit, width=13, height=8.5, units="in")
ggsave("association_studies_results_plots/all_alleles_six_digit.pdf", all_alleles_six_digit, width=13, height=8.5, units="in")
