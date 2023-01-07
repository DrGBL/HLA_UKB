library(tidyverse)
library(vroom)
library(logistf)
library(ggpubr)

pathComp<-"/path/to/comparisons_with_imputed/"
pathCalls<-"/path/to/calls/"
path_anc<-"/path/to/folder/with/ancestries/"

#ancestries
#each file is a one row list of ID assigned to each genetic ancestry
list_anc<-c()
list_anc[["afr"]]<-scan(paste0(path_anc, "ukb.afrIDsPCA.txt"))
list_anc[["amr"]]<-scan(paste0(path_anc, "ukb.amrIDsPCA.txt"))
list_anc[["eas"]]<-scan(paste0(path_anc, "ukb.easIDsPCA.txt"))
list_anc[["eur"]]<-scan(paste0(path_anc, "ukb.eurIDsPCA.txt"))
list_anc[["sas"]]<-scan(paste0(path_anc, "ukb.sasIDsPCA.txt"))

comp<-c()
calls<-c()
for(folder in 10:60){
  comp<-vroom(paste0(pathComp, "comparisons_imputed_called_folder_", folder, ".tsv.gz"),
              col_types = cols(.default = "c")) %>%
    rename(ID=ID_col) %>%
    mutate(across(.cols=-c("ID"), 
                  ~ ifelse(.x %in% c("Two matches",
                                     "Both uncalled and unimputed",
                                     "One match, other uncalled and unimputed"),
                           "Concordant",
                           ifelse(startsWith(.x, "Both uncalled"), "Uncalled", "Discordant")))) %>%
    bind_rows(comp,.)
  
  calls<-vroom(paste0(pathCalls, "hla_df_batch_qced_", folder, ".tsv.gz"),
               col_types = cols(.default = "c"),
               col_select = c("ID", 
                              paste0(rep(c("A", 
                                           "B", 
                                           "C", 
                                           "DRB1", 
                                           "DQA1", 
                                           "DQB1", 
                                           "DPA1", 
                                           "DPB1", 
                                           "DRB3", 
                                           "DRB4", 
                                           "DRB5"),
                                         each=2),
                                     c("_1", "_2")))) %>%
    mutate(A_1=ifelse(A_1==0, "A*99:99:99", A_1)) %>%
    mutate(A_2=ifelse(A_2==0, "A*99:99:99", A_2)) %>%
    mutate(B_1=ifelse(B_1==0, "B*99:99:99", B_1)) %>%
    mutate(B_2=ifelse(B_2==0, "B*99:99:99", B_2)) %>%
    mutate(C_1=ifelse(C_1==0, "C*99:99:99", C_1)) %>%
    mutate(C_2=ifelse(C_2==0, "C*99:99:99", C_2)) %>%
    mutate(DRB1_1=ifelse(DRB1_1==0, "DRB1*99:99:99", DRB1_1)) %>%
    mutate(DRB1_2=ifelse(DRB1_2==0, "DRB1*99:99:99", DRB1_2)) %>%
    mutate(DQA1_1=ifelse(DQA1_1==0, "DQA1*99:99:99", DQA1_1)) %>%
    mutate(DQA1_2=ifelse(DQA1_2==0, "DQA1*99:99:99", DQA1_2)) %>%
    mutate(DQB1_1=ifelse(DQB1_1==0, "DQB1*99:99:99", DQB1_1)) %>%
    mutate(DQB1_2=ifelse(DQB1_2==0, "DQB1*99:99:99", DQB1_2)) %>%
    mutate(DPA1_1=ifelse(DPA1_1==0, "DPA1*99:99:99", DPA1_1)) %>%
    mutate(DPA1_2=ifelse(DPA1_2==0, "DPA1*99:99:99", DPA1_2)) %>%
    mutate(DPB1_1=ifelse(DPB1_1==0, "DPB1*99:99:99", DPB1_1)) %>%
    mutate(DPB1_2=ifelse(DPB1_2==0, "DPB1*99:99:99", DPB1_2)) %>%
    mutate(DRB3_1=ifelse(DRB3_1==0, "DRB3*99:99:99", DRB3_1)) %>%
    mutate(DRB3_2=ifelse(DRB3_2==0, "DRB3*99:99:99", DRB3_2)) %>%
    mutate(DRB4_1=ifelse(DRB4_1==0, "DRB4*99:99:99", DRB4_1)) %>%
    mutate(DRB4_2=ifelse(DRB4_2==0, "DRB4*99:99:99", DRB4_2)) %>%
    mutate(DRB5_1=ifelse(DRB5_1==0, "DRB5*99:99:99", DRB5_1)) %>%
    mutate(DRB5_2=ifelse(DRB5_2==0, "DRB5*99:99:99", DRB5_2)) %>%
    mutate(A=ifelse(A_1>A_2, paste0(A_2, "_", A_1), paste0(A_1, "_", A_2))) %>%
    mutate(B=ifelse(B_1>B_2, paste0(B_2, "_", B_1), paste0(B_1, "_", B_2))) %>%
    mutate(C=ifelse(C_1>C_2, paste0(C_2, "_", C_1), paste0(C_1, "_", C_2))) %>%
    mutate(DRB1=ifelse(DRB1_1>DRB1_2, paste0(DRB1_2, "_", DRB1_1), paste0(DRB1_1, "_", DRB1_2))) %>%
    mutate(DRB3=ifelse(DRB3_1>DRB3_2, paste0(DRB3_2, "_", DRB3_1), paste0(DRB3_1, "_", DRB3_2))) %>%
    mutate(DRB4=ifelse(DRB4_1>DRB4_2, paste0(DRB4_2, "_", DRB4_1), paste0(DRB4_1, "_", DRB4_2))) %>%
    mutate(DRB5=ifelse(DRB5_1>DRB5_2, paste0(DRB5_2, "_", DRB5_1), paste0(DRB5_1, "_", DRB5_2))) %>%
    mutate(DQA1=ifelse(DQA1_1>DQA1_2, paste0(DQA1_2, "_", DQA1_1), paste0(DQA1_1, "_", DQA1_2))) %>%
    mutate(DQB1=ifelse(DQB1_1>DQB1_2, paste0(DQB1_2, "_", DQB1_1), paste0(DQB1_1, "_", DQB1_2))) %>%
    mutate(DPA1=ifelse(DPA1_1>DPA1_2, paste0(DPA1_2, "_", DPA1_1), paste0(DPA1_1, "_", DPA1_2))) %>%
    mutate(DPB1=ifelse(DPB1_1>DPB1_2, paste0(DPB1_2, "_", DPB1_1), paste0(DPB1_1, "_", DPB1_2))) %>%
    bind_rows(calls,.)
}

calls<-calls %>%
  filter(ID %in% comp$ID) %>%
  arrange(ID)

comp<-comp %>%
  filter(ID %in% calls$ID) %>%
  arrange(ID)


create_df<-function(gene, calls=calls, comp=comp, anc=NA){
  if(is.na(anc)){
    calls_mod<-calls
    comp_mod<-comp
  } else {
    calls_mod<-calls %>%
      filter(ID %in% list_anc[[anc]])
    comp_mod<-comp %>%
      filter(ID %in% list_anc[[anc]])
  }

  table_gene_comb<-table(calls_mod[,gene]) %>%
    as.data.frame() %>%
    rename(Alleles:=Var1) %>%
    mutate(Concordance=NA)

  table_gene<-c(calls_mod %>% pull(paste0(gene, "_1")), calls_mod %>% pull(paste0(gene, "_2"))) %>%
    table() %>%
    as.data.frame() %>%
    rename(Alleles=".")

  table_allele1<-table_gene %>%
    rename(!!paste0(gene, "_1"):=Alleles) %>%
    rename(Freq1=Freq)

  table_allele2<-table_gene %>%
    rename(!!paste0(gene, "_2"):=Alleles) %>%
    rename(Freq2=Freq)

  table_gene_comb<-table_gene_comb %>%
    mutate(tmp=Alleles) %>%
    separate(col=tmp, into=c(paste0(gene, "_1"), paste0(gene, "_2")), sep="_") %>%
    merge(.,table_allele1) %>%
    merge(.,table_allele2)

  for(i in 1:nrow(table_gene_comb)){
    list_alleles<-calls_mod %>%
      filter(!!sym(gene)==table_gene_comb$Alleles[i]) %>%
      pull(ID)
    table_gene_comb$Concordance[i]<-comp %>%
      filter(ID %in% list_alleles) %>%
      filter(!!sym(gene)=="Concordant") %>%
      nrow()
  }

  if(is.na(anc)){
    saveRDS(table_gene_comb, paste0(pathComp, "table_", gene,".rds"))
  } else {
    saveRDS(table_gene_comb, paste0(pathComp, "table_", gene,"_", anc, ".rds"))
  }
}


for(g in c("A",
           "B",
           "C",
           "DRB1",
           "DQA1",
           "DQB1",
           "DPA1",
           "DPB1",
           "DRB3",
           "DRB4",
           "DRB5")){
  #for(anc in c("afr", "amr", "eas", "eur", "sas")){
  for(anc in c(NA, "afr", "amr", "eas", "eur", "sas")){
    create_df(gene=g, calls=calls, comp=comp, anc=anc)
  }
}

table_df<-c()
for(g in c("A", 
           "B", 
           "C", 
           "DRB1", 
           "DQA1", 
           "DQB1", 
           "DPA1", 
           "DPB1", 
           "DRB3", 
           "DRB4", 
           "DRB5")){
  #for(anc in c("afr", "amr", "eas", "eur", "sas")){
  for(anc in c(NA, "afr", "amr", "eas", "eur", "sas")){
    if(is.na(anc)){
      table_df<-readRDS(paste0(pathComp, "table_", g,".rds")) %>%
        mutate(Gene=paste0("HLA-",g)) %>%
        mutate(Mean_freq= ( ( Freq1/(2*nrow(calls)) )  + ( Freq2/(2*nrow(calls)) )  )/2   ) %>%
        mutate(Concordance_rate=Concordance/Freq) %>%
        mutate(Ancestry="all") %>%
        dplyr::select(c(Alleles, Freq, Concordance, Freq1, Freq2, Gene, Mean_freq, Concordance_rate, Ancestry)) %>%
        bind_rows(table_df,.)
    } else {
      table_df<-readRDS(paste0(pathComp, "table_", g,"_", anc, ".rds")) %>%
        mutate(Gene=paste0("HLA-",g)) %>%
        mutate(Mean_freq= ( ( Freq1/(2*nrow(calls %>% filter(ID %in% list_anc[[anc]]))) )  + ( Freq2/(2*nrow(calls %>% filter(ID %in% list_anc[[anc]]))) )  )/2   ) %>%
        mutate(Concordance_rate=Concordance/Freq) %>%
        mutate(Ancestry=anc) %>%
        dplyr::select(c(Alleles, Freq, Concordance, Freq1, Freq2, Gene, Mean_freq, Concordance_rate, Ancestry)) %>%
        bind_rows(table_df,.)
    }
  }
}




anc_gene_long_function<-function(table_df=table_df, calls=calls, comp=comp, list_anc=list_anc, gene, anc){
  if(anc=="all"){
    comp_tmp<-comp %>% dplyr::select(c(ID, all_of(gene))) %>% rename(comp_gene=gene)
    calls_tmp<-calls
    table_df_tmp<-table_df
  } else {
    comp_tmp<-comp %>% dplyr::select(c(ID, all_of(gene))) %>% filter(ID %in% list_anc[[anc]]) %>% rename(comp_gene=all_of(gene))
    calls_tmp<-calls %>% filter(ID %in% list_anc[[anc]])
    table_df_tmp <-table_df %>%
      filter(Ancestry==anc)
  }
  
  print(dim(calls_tmp))
  long_df<-table_df_tmp %>%
    filter(Gene==paste0("HLA-", gene)) %>%
    rename(!!sym(gene):=Alleles) %>%
    rename(Mean_freq_gene=Mean_freq) %>%
    dplyr::select(c(all_of(gene), Mean_freq_gene)) %>%
    merge(calls_tmp,.) %>%
    merge(comp_tmp,.)
  
  return(long_df %>% dplyr::select(ID, comp_gene, Mean_freq_gene) %>% mutate(Gene=gene) %>% mutate(Ancestry=anc))
}


list_plot_smooth<-c()
for(g in c("A", "B", "C", "DRB1",  "DQA1", "DQB1", "DPA1", "DPB1", "DRB3", "DRB4", "DRB5")){
  df_data<-c()
  for(a in c("all", "afr", "amr", "eas", "eur", "sas")){
    df_data<-anc_gene_long_function(table_df=table_df, gene=g, anc=a, comp=comp, calls=calls, list_anc=list_anc) %>%
      bind_rows(df_data,.)
  }
  
  df_data<-df_data %>%
    mutate(Ancestry=toupper(Ancestry))
  
  list_plot_smooth[[g]]<-df_data %>%
    mutate(comp_gene=ifelse(comp_gene=="Concordant",
                            1,
                            ifelse(comp_gene=="Uncalled",NA,0))) %>%
    ggplot(aes(x=Mean_freq_gene, y=comp_gene)) +
    geom_smooth(aes(colour=Ancestry), method = "glm", method.args = list(family = "binomial"))+
    ggtitle(paste0("HLA-",g)) +
    theme_bw() +
    xlab("Mean Allele Frequency")+
    ylab("Predicted Concordance")
}

glm_plots<-ggarrange(plotlist=list_plot_smooth, ncol=4, nrow=3,common.legend = TRUE)

ggsave(paste0(pathComp, "concord_vs_af.pdf"), glm_plots, height = 8.5, width=11)
