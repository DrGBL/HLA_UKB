#go to the directory where the folders "hla_calling_hla_hd_batch_X" are located
setwd("/your/local/folder")
#now also create a "QC" and "calls" folder
if(!dir.exists("QC")){
  dir.create("QC")
}
if(!dir.exists("calls")){
  dir.create("calls")
}


library(tidyverse)
library(vroom)
library(ff)

hla_alleles<-as.list(c(1:31))
names(hla_alleles)<-c("A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V", "Y",
                      "DMA", "DMB", "DOA", "DOB", "DPA1", "DPA2", "DPB1", "DQA1", 
                      "DQB1", "DRA", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB6",
                      "DRB7", "DRB8", "DRB9")

for(folder in 10:60){
 
  samples<-list.dirs(paste0("hla_calling_hla_hd_batch_", folder, "/"),
                     recursive=FALSE) %>%
    data.frame(samples=.) %>%
    mutate(samples=substr(samples,
                          start=nchar(paste0("hla_calling_hla_hd_batch_", folder, "//"))+1,
                          stop=1000000))
  
  ########### prep the hla calling results #############
  
  #first, pass through all files, to get the list of all genes called, to be 100% sure we're not making mistakes later
  #this also unzips and renames the files to what will be provided to the UK Biobank
  genes_hla<-data.frame(genes=NA_character_)
  for(sample in samples$samples){
    #print(sample)
    
    setwd(paste0("hla_calling_hla_hd_batch_", folder, "/", sample, "/"))
    untar("hla_res.tar.gz")
    file.move(from = "test_call/result", to = paste0(getwd(),"/", "result"))
    unlink("test_call", recursive = TRUE)
    setwd("result")
    old_files <- list.files(full.names = TRUE)
    new_files<- gsub("test_call", sample, old_files)
    file.copy(from = old_files, to = new_files)
    file.remove(old_files)
    
    genes_hla<-read.table(paste0(sample,
                                 "_final.result.txt"),
                          header = FALSE, 
                          sep = "", 
                          col.names = paste0("X",seq_len(1000)), 
                          fill = TRUE,
                          allowEscapes = FALSE) %>%
      dplyr::select(X1) %>%
      dplyr::rename(genes=X1) %>%
      rbind(genes_hla) %>%
      distinct()
    
    setwd(paste0("/your/local/folder/hla_calling_hla_hd_batch_", folder, "/", sample, "/"))
    
    tar(paste0(sample,".tar.gz"), "result", compression = 'gzip', tar="tar")
    file.remove("hla_res.tar.gz")
    unlink("result", recursive = TRUE)
    
    setwd("/your/local/folder/ukb_wes")

    if(nrow(genes_hla) > 34){
      break
    }
    
  }
  
  genes_hla <- genes_hla %>% filter(!is.na(genes) & genes!="Couldn't")
  vroom_write(genes_hla, paste0("genes_hla_", folder,".tsv"))
  #if need to re-extract
  #genes_hla<-read_tsv(paste0("genes_hla_", folder,".tsv"))
  
  
  #now set up data frame
  hla_matrix<-matrix(NA, nrow=nrow(samples), ncol=2*nrow(genes_hla)+1)
  numbering_vector<-rep(c(1,2), nrow(genes_hla))
  gene_names<-paste0(rep(genes_hla$genes,
                         each=2),
                     "_",
                     numbering_vector)
  colnames(hla_matrix)<-c("ID", gene_names)
  hla_df<-as.data.frame(hla_matrix)
  
  hla_matrix_coverage<-matrix(NA, nrow=nrow(samples), ncol=nrow(genes_hla)+1)
  colnames(hla_matrix_coverage)<-c("ID", genes_hla$genes)
  hla_coverage<-as.data.frame(hla_matrix_coverage) %>%
    select(-c("Couldn't"))
  
  hla_drb<-matrix(NA, nrow=nrow(samples), ncol=7)
  numbering_vector_drb<-rep(c(1,2), 3)
  drb_names<-paste0(rep(c("DRB3", "DRB4", "DRB5"),
                         each=2),
                     "_",
                     numbering_vector_drb)
  colnames(hla_drb)<-c("ID", drb_names)
  hla_drb_coverage<-as.data.frame(hla_drb)
  
  #note that this is using the best called alleles at 6 digits. It ignores the other candidates.
  list_error<-list()
  for(s in 1:nrow(samples)){
    if((s%%100==0)){
      print(paste0(s, ", folder ", folder))
    }
    
    setwd(paste0("hla_calling_hla_hd_batch_", folder, "/",
                 samples$samples[s],
                 "/"))
    
    untar(paste0(samples$samples[s],
                 ".tar.gz"))
    
    tmp<-read.table(paste0("result/",
                           samples$samples[s],
                           "_final.result.txt"),
                    header = FALSE,
                    sep = "",
                    col.names = paste0("X",seq_len(1000)),
                    fill = TRUE,
                    allowEscapes = FALSE) %>%
      select(c(X1,X2,X3)) %>%
      filter(X1!="Couldn't")
    
    hla_df$ID[s]<-samples$samples[s]
    
    for(g in genes_hla$genes){
      if(length(tmp[which(tmp$X1 == g), 1])>0){
        if(!(tmp[which(tmp$X1==g),2] %in% c("Not", "typed", "read", "result"))){
          hla_df[s,paste0(g, "_1")] <- tmp[which(tmp$X1==g),2]
        }
        if(!(tmp[which(tmp$X1==g),3] %in% c("Not", "typed", "read", "result"))){
          if(tmp[which(tmp$X1==g),3]=="-"){
            hla_df[s,paste0(g, "_2")] <- tmp[which(tmp$X1==g),2]
          } else {
            hla_df[s,paste0(g, "_2")] <- tmp[which(tmp$X1==g),3]
          }
        }
      }
    }
    
    #now do the coverage
    #for drb2 and drb8, exon3 was used. exon2 was used for all others
    for(gene in colnames(hla_coverage)[-1]){
      #print(gene)
      if(gene %in% c("DRB2", "DRB8")){
        exon<-"exon3"
       } else {
        exon<-"exon2"
       }
      hla_coverage$ID[s]<-samples$samples[s]
      if(file.exists(paste0("result/",
                            samples$samples[s],
                            "_",
                            gene,
                            ".est.txt"))){
        hla_coverage[s,gene]<-read_file(paste0("result/",
                                               samples$samples[s],
                                               "_",
                                               gene,
                                               ".est.txt")) %>%
          data.frame(tmp_cov=.) %>%
          separate_rows(tmp_cov, sep=" ") %>%
          separate_rows(tmp_cov, sep="\t") %>%
          separate_rows(tmp_cov, sep="\n") %>%
          separate_rows(tmp_cov, sep=",") %>%
          filter(str_detect(tmp_cov, exon)) %>%
          separate(tmp_cov, into=c("exon", "coverage", "comp"), sep=":") %>%
          mutate(coverage=as.numeric(coverage)) %>%
          summarize(sum(coverage))
      } else {
        list_error[[samples$samples[s]]]<-c(list_error[[samples$samples[s]]],gene)
      }
    }
   
    #for DRB3-4-5, extract coverage for both allele, not just mean
    exon<-"exon2"
    for(gene in c("DRB3", "DRB4", "DRB5")){
      #print(gene)
      hla_drb_coverage$ID[s]<-samples$samples[s]
      if(file.exists(paste0("result/",
                            samples$samples[s],
                            "_",
                            gene,
                            ".est.txt"))){
        tmp_cov_drb<-read_file(paste0("result/",
                                               samples$samples[s],
                                               "_",
                                               gene,
                                               ".est.txt")) %>%
          data.frame(tmp_cov=.) %>%
          separate_rows(tmp_cov, sep=" ") %>%
          separate_rows(tmp_cov, sep="\t") %>%
          separate_rows(tmp_cov, sep="\n") %>%
          separate_rows(tmp_cov, sep=",") %>%
          filter(str_detect(tmp_cov, exon)) %>%
          separate(tmp_cov, into=c("exon", "coverage", "comp"), sep=":") %>%
          mutate(coverage=as.numeric(coverage))
          
        hla_drb_coverage[s,paste0(gene, "_1")] <- tmp_cov_drb$coverage[1]
        hla_drb_coverage[s,paste0(gene, "_2")] <- tmp_cov_drb$coverage[2]
        
        
      } else {
        list_error[[samples$samples[s]]]<-c(list_error[[samples$samples[s]]],gene)
      }
    }
    
    #now remove the recently untarred tar ball
    unlink("result", recursive = TRUE)
    
    setwd("/your/local/folder/ukb_wes")
    
  }

  vroom_write(hla_df, paste0("calls/hla_df_batch_", folder, ".tsv.gz"))
  vroom_write(hla_coverage, paste0("QC/hla_coverage_batch_", folder, ".tsv.gz"))
  vroom_write(hla_drb_coverage, paste0("QC/hla_drb_coverage_batch_", folder, ".tsv.gz"))
  
  genes_missing<-data.frame(samples=names(list_error)) %>%
    distinct()
  vroom_write(genes_missing, paste0("genes_missing/genes_missing_folder_", folder, ".tsv.gz"))
  
  saveRDS(list_error, paste0("genes_missing/list_genes_missing_folder_", folder, ".rds"))
  
}
