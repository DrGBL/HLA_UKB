library(tidyverse)
library(vroom)
library(asymLD)
library(corrplot)
library(readxl)
library(writexl)

path_anc<-"/path/to/ancestry/"

#ancestries, these are files with one row per ID in each genetic ancestry
afr<-scan(paste0(path_anc, "ukb.afrIDsPCA.txt"))
amr<-scan(paste0(path_anc, "ukb.amrIDsPCA.txt"))
eas<-scan(paste0(path_anc, "ukb.easIDsPCA.txt"))
eur<-scan(paste0(path_anc, "ukb.eurIDsPCA.txt"))
sas<-scan(paste0(path_anc, "ukb.sasIDsPCA.txt"))


#the function that calls asymmetric LD
prep_asymld<-function(gene1, gene2, hla_df_calls){
  
  low_freq_gene1<-as.data.frame(table(c(hla_df_calls[,paste0(gene1, "_1")], 
                                        hla_df_calls[, paste0(gene1, "_2")]))) %>%
    mutate(Freq=Freq/(2*nrow(hla_df_calls))) %>%
    filter(Freq<0.01)
  
  low_freq_gene2<-as.data.frame(table(c(hla_df_calls[,paste0(gene2, "_1")], 
                                        hla_df_calls[, paste0(gene2, "_2")]))) %>%
    mutate(Freq=Freq/(2*nrow(hla_df_calls))) %>%
    filter(Freq<0.01)
  
  #here you assign all alleles with frequency less than 1% to a dummy variable (as per https://doi.org/10.1038/s41588-018-0336-0)
  #similarly, assign uncalled alleles to the same dummy variable (99:99)
  #also remove any letters at the end of each allele
  
  hla_df_calls_fixed<-hla_df_calls %>%
    mutate(!!sym(paste0(gene1, "_1")):=ifelse(!!sym(paste0(gene1, "_1")) %in% low_freq_gene1$Var1,
                                              paste0(gene1, "*99:99"),
                                              !!sym(paste0(gene1, "_1")))) %>%
    mutate(!!sym(paste0(gene1, "_2")):=ifelse(!!sym(paste0(gene1, "_2")) %in% low_freq_gene1$Var1,
                                              paste0(gene1, "*99:99"),
                                              !!sym(paste0(gene1, "_2")))) %>%
    mutate(!!sym(paste0(gene2, "_1")):=ifelse(!!sym(paste0(gene2, "_1")) %in% low_freq_gene2$Var1,
                                              paste0(gene2, "*99:99"),
                                              !!sym(paste0(gene2, "_1")))) %>%
    mutate(!!sym(paste0(gene2, "_2")):=ifelse(!!sym(paste0(gene2, "_2")) %in% low_freq_gene2$Var1,
                                              paste0(gene2, "*99:99"),
                                              !!sym(paste0(gene2, "_2")))) %>%
    mutate(!!sym(paste0(gene1, "_1")):=ifelse(is.na(!!sym(paste0(gene1, "_1"))),
                                              paste0(gene1, "*99:99"),
                                              !!sym(paste0(gene1, "_1")))) %>%
    mutate(!!sym(paste0(gene1, "_2")):=ifelse(is.na(!!sym(paste0(gene1, "_2"))),
                                              paste0(gene1, "*99:99"),
                                              !!sym(paste0(gene1, "_2")))) %>%
    mutate(!!sym(paste0(gene2, "_1")):=ifelse(is.na(!!sym(paste0(gene2, "_1"))),
                                              paste0(gene2, "*99:99"),
                                              !!sym(paste0(gene2, "_1")))) %>%
    mutate(!!sym(paste0(gene2, "_2")):=ifelse(is.na(!!sym(paste0(gene2, "_2"))),
                                              paste0(gene2, "*99:99"),
                                              !!sym(paste0(gene2, "_2")))) %>%
    mutate(!!sym(paste0(gene1, "_1")):=gsub("[A-Z]*$","",!!(sym(paste0(gene1, "_1"))))) %>%
    mutate(!!sym(paste0(gene1, "_2")):=gsub("[A-Z]*$","",!!(sym(paste0(gene1, "_2"))))) %>%
    mutate(!!sym(paste0(gene2, "_1")):=gsub("[A-Z]*$","",!!(sym(paste0(gene2, "_1"))))) %>%
    mutate(!!sym(paste0(gene2, "_2")):=gsub("[A-Z]*$","",!!(sym(paste0(gene2, "_2")))))
  
  
  comb<-c()
  for(chr1 in 1:2){
    for(chr2 in 1:2){
      comb<-c(comb,paste0(hla_df_calls_fixed[,paste0(gene1,"_", chr1)], "-", hla_df_calls_fixed[,paste0(gene2,"_", chr2)]))
    }
  }
  tally_comb<-as.data.frame(table(comb))
  
  asymLD_ready_df<-data.frame(haplo.freq=tally_comb$Freq/(4*nrow(hla_df_calls_fixed)), 
                              locus1=rep(gene1, nrow(tally_comb)), 
                              locus2=rep(gene2, nrow(tally_comb)),
                              allele1=str_extract(tally_comb$comb, "[A-Z,0-9]*\\*[0-9]*:[0-9]*-"),
                              allele2=str_extract(tally_comb$comb, "-[A-Z,0-9]*\\*[0-9]*:[0-9]*")) %>%
    mutate(allele1=gsub("-", "", allele1)) %>%
    mutate(allele2=gsub("-", "", allele2))
    
  return(asymLD_ready_df)
}

#rest of code, this first section only needs to be ran once
setwd("/project/richards/guillaume.butler-laporte/HLA/ukb_wes")

dp_threshold<-10


hla_alleles<-as.list(c(1:31))
names(hla_alleles)<-c("A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V", "Y",
                      "DMA", "DMB", "DOA", "DOB", "DPA1", "DPA2", "DPB1", "DQA1",
                      "DQB1", "DRA", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB6",
                      "DRB7", "DRB8", "DRB9")

hla_four_full_fixed<-c()

for(folder in 10:60){
  hla_df<-vroom(paste0("calls/hla_df_batch_qced", folder, ".tsv.gz"), col_types = cols(.default = "c"))
  hla_qc<-vroom(paste0("QC/hla_coverage_batch_", folder, ".tsv.gz"), col_types = cols(.default = "c"))

  #remove the HLA prefix if needed
  hla_df <- data.frame(lapply(hla_df, function(x){gsub("HLA-", "", x)}))

  hla_df_four <- lapply(hla_df, function(x) ifelse(str_count(x,":")>1, sub(":[0-9A-Z]*$", "", x), x)) %>%
    as.data.frame()

  hla_four_full_fixed<-bind_rows(hla_four_full_fixed, hla_df_four)

}

#if needed
#saveRDS(hla_four_full_fixed, "/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code/hla_four_full_fixed.RDS")

rm(hla_df_four_fixed)
rm(hla_df_four)
rm(hla_df)
rm(hla_qc)

#following also ran only once and saved (it calls the big funtion above for each ancestry)
hla_four_full_fixed<-readRDS("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code/hla_four_full_fixed.RDS") %>%
  dplyr::select(-paste0(rep(c("DRB2", "DRB3", "DRB4", "DRB5", "DRB6", "DRB7", "DRB8", "DRB9"), each=2), c("_1", "_2"))) %>%
  mutate(across(-c(ID), na_if, y="0"))

hla_genes<-substr(colnames(hla_four_full_fixed)[-1], 1, nchar(colnames(hla_four_full_fixed)[-1])-2)[seq(from=2, to=length(colnames(hla_four_full_fixed)[-1]), by=2)]

df_asymld<-as.data.frame(matrix(NA, nrow=length(hla_genes), ncol=length(hla_genes)))
colnames(df_asymld)<-hla_genes
rownames(df_asymld)<-hla_genes

list_df_asymld<-list("all", "afr", "amr", "eas", "eur", "sas")

all<-hla_four_full_fixed %>% pull(ID)

for(anc in c("all", "afr", "amr", "eas", "eur", "sas")){
  for(r in 1:length(hla_genes)){ #rows
    print(paste0("Ancestry: ", anc, ". Row: ", r, "/", length(hla_genes), "."))
    for(c in r:length(hla_genes)){ #columns
      if(r==c){
        df_asymld[r,c]<-1
      } else {
        prepped_df<-prep_asymld(gene1=hla_genes[r],
                                gene2=hla_genes[c],
                                hla_df_calls = hla_four_full_fixed %>% filter(ID %in% !!sym(anc)))
        suppressWarnings(ld_res<-prepped_df %>% filter(haplo.freq>0.0000001) %>% compute.ALD(., tolerance=1))
        #ALD.x.y is asymmetric LD for locus x conditioned on locus y
        #so locus x will be rows, and locus y will be columns
        df_asymld[hla_genes[r], hla_genes[c]]<-ld_res$ALD.1.2
        df_asymld[hla_genes[c], hla_genes[r]]<-ld_res$ALD.2.1
      }
    }
  }
  list_df_asymld[[anc]]<-df_asymld
}

#if needed
#saveRDS(list_df_asymld, "/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code/list_df_asymld.RDS")

#restart here for the plots if saved earlier
#list_df_asymld<-readRDS("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code/list_df_asymld.RDS")


pdf("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code/ld_plots/fig3_ld.pdf", height = 11, width=8.5)
par(mfrow=c(3,2))
for(anc in c("all", "afr", "amr", "eas", "eur", "sas")){
  df_asymld<-list_df_asymld[[anc]]
  df_asymld<-df_asymld[which(!is.na(df_asymld$A)), which(!is.na(df_asymld$A))]
  
  write_xlsx(df_asymld, paste0("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/regenie_code/ld_plots/", anc, "_asym_ld.xlsx"))
  
  df_upper<-as.matrix(df_asymld)
  df_upper[lower.tri(df_upper)]<-t(df_upper)[lower.tri(df_upper)]
  
  df_lower<-as.matrix(df_asymld)
  df_lower[upper.tri(df_lower)]<-t(df_lower)[upper.tri(df_lower)]
  
  df_mean<-(df_upper+df_lower)/2
  
  col_lab<-data.frame(gene=colnames(df_mean)) %>%
    mutate(class=ifelse(gene %in% c("A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "Y"), 
                        "#ef8a62", 
                        "#67a9cf"))
  
  if(anc=="all"){
    num_box<-5
    title_anc<-"All ancestries"
  }
  if(anc=="afr"){
    num_box<-5
    title_anc<-"AFR"
  }
  if(anc=="amr"){
    num_box<-5
    title_anc<-"AMR"
  }  
  if(anc=="eas"){
    num_box<-5
    title_anc<-"EAS"
  }
  if(anc=="eur"){
    num_box<-5
    title_anc<-"EUR"
  }
  if(anc=="sas"){
    num_box<-5
    title_anc<-"SAS"
  }
  
  clust_order<-hclust(as.dist(1-df_mean), method="average")$order
  
  corrplot(df_mean, 
           tl.col = col_lab$class[clust_order],
           order = 'hclust', 
           hclust.method="average",
           addrect = num_box, 
           col.lim=c(0,1), 
           is.corr=FALSE,
           col=COL1(sequential="Blues", n=100),
           title=title_anc,
           mar=c(0,0,1,0))
}
dev.off()


df_asymld<-df_asymld[which(!is.na(df_asymld$A)), which(!is.na(df_asymld$A))]

df_upper<-as.matrix(df_asymld)
df_upper[lower.tri(df_upper)]<-t(df_upper)[lower.tri(df_upper)]

df_lower<-as.matrix(df_asymld)
df_lower[upper.tri(df_lower)]<-t(df_lower)[upper.tri(df_lower)]

df_mean<-(df_upper+df_lower)/2

#df_mean
corrplot(df_mean, 
         tl.col = 'black',
         order = 'hclust', 
         addrect = 5, 
         col.lim=c(0,1), 
         is.corr=FALSE,
         col=COL1(sequential="Blues", n=100))

#df_lower
corrplot(df_lower, 
         tl.col = 'black',
         order = 'hclust', 
         col.lim=c(0,1), 
         is.corr=FALSE,
         col=COL1(sequential="Blues", n=100),
         type="lower",
         tl.srt = 45)

#df_upper
corrplot(df_upper, 
         tl.col = 'black',
         order = 'hclust', 
         col.lim=c(0,1), 
         is.corr=FALSE,
         col=COL1(sequential="Blues", n=100),
         type="upper",
         tl.srt = 45)

