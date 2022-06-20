library(tidyverse)
library(vroom)

#still at the clone git directory in step 04
setwd("/path/to/IMGTHLA/")

#will create a folder where the properly aligned amino acids sequences will be output, for further analyses
if(!dir.exists("alignments_aa_fixed")){
  dir.create("alignments_aa_fixed")
}

for(gene in c("A", "B", "C", "DMA", "DMB", "DOA", "DOB", "DPA1", "DPA2", "DPB1", "DQA1", "DQB1", "DRA", "DRB", "E", "F", "G", "H", "J", "K", "L", "Y")){
#for(gene in "C"){
  if(file.exists(paste0("tmp_aa_align/",gene,1,".txt"))){
    #gene<-"A"
    
    pos_1<-read_file(paste0("tmp_aa_align/pos_string_", gene, ".txt")) %>%
      substr(1,nchar(.)-1) %>%
      substr(18,nchar(.)) %>%
      gsub(" ", "",.) %>%
      nchar(.)
    
    pos_1_no_indel<-read_file(paste0("tmp_aa_align/pos_string_", gene, ".txt")) %>%
      substr(1,nchar(.)-1) %>%
      substr(18,nchar(.)) %>%
      gsub(" ", "",.) %>%
      gsub("\\.","",.) %>%
      nchar(.)
    
    n_indel_space<-pos_1-pos_1_no_indel
    
    pos_first_aa<-ifelse(pos_1>1,-pos_1_no_indel+1,pos_1)
    
    align<-vroom_fwf(paste0("tmp_aa_align/",gene,1,".txt"), col_types = cols(.default = "c")) %>%
      rename(allele=X1) %>%
      mutate(allele=gsub(" ", "", allele)) %>%
      filter(!is.na(allele)) %>%
      unite(col=!!paste0("align",1),-allele,sep="") %>%
      distinct() %>%
      mutate(tmp_order=c(1:nrow(.)))
    
    for(file in 2:10){
      if(file.exists(paste0("tmp_aa_align/",gene,file,".txt"))){
        tmp<-vroom_fwf(paste0("tmp_aa_align/",gene,file,".txt"), col_types = cols(.default = "c")) %>%
          rename(allele=X1) %>%
          mutate(allele=gsub(" ", "", allele)) %>%
          filter(!is.na(allele))
        if(ncol(tmp)>1){
          align<-tmp %>%
            unite(col=!!paste0("align",file),-allele,sep="") %>%
            distinct() %>%
            merge(align,., all.x=TRUE) %>%
            arrange(tmp_order)
        }
      }
    }
    
    align<-align %>% 
      dplyr::select(-tmp_order) %>%
      unite(col=alignment, -allele, sep="") %>%
      separate(col=alignment, sep=c(1:(max(nchar(.$alignment)-1))), into=paste0("X",c(1:max(nchar(.$alignment))))) %>%
      mutate(across(-allele,~ifelse(.x=="-",.x[1],.x))) 
    
    aa_pos<-rep(NA,ncol(align)-1)
    add_pos<-0
    indel_pos<-0
    for(i in 1:length(aa_pos)){
      if(align[1,i+1]!="."){
        if(pos_first_aa+add_pos>0){
          aa_pos[i]<-paste0(gene,"_",pos_first_aa+add_pos)
          add_pos<-add_pos+1
        } else {
          aa_pos[i]<-paste0(gene,pos_first_aa+add_pos)
          add_pos<-add_pos+1
          if(pos_first_aa+add_pos==0){
            add_pos<-add_pos+1
          }
        }
      } else {
        aa_pos[i]<-paste0(gene,"_indel",pos_first_aa+add_pos+indel_pos)
        if(align[1,i+2]=="."){
          indel_pos<-indel_pos+1
        } else {
          indel_pos<-0
        }
      }
    }
    
    colnames(align)[-1]<-aa_pos
    
    align<-Filter(function(x) length(unique(x))>1, align)
    
    vroom_write(x=align, file=paste0("alignments_aa_fixed/", gene, "_fixed.tsv.gz"))
    
    #now do 4 digits only
    align_four<-align %>%
      mutate(tmp_order=c(1:nrow(.))) %>%
      mutate(allele=str_extract(allele,"[A-Z0-9]*\\*[0-9]*:[0-9]*")) %>% 
      group_by(allele) %>% 
      filter(row_number()==1) %>%
      ungroup() %>%
      arrange(tmp_order) %>%
      dplyr::select(-tmp_order)
    
    align_four<-Filter(function(x) length(unique(x))>1, align_four)
    
    vroom_write(x=align_four, file=paste0("alignments_aa_fixed/", gene, "_four_fixed.tsv.gz"))
  }
}
