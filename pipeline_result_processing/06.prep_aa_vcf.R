library(tidyverse)
library(vroom)
library(stringi)

#note: I assume that unknown amino acids ("*") are from reference (i.e. set to 0).

#here you can adjust the different paths
#this is the path where you can find the HLA allele calls
pathCalls<-"/path/to/folder/calls/"
#this is the path to the fixed amino acid alignments from step 05
pathSeq<-"/path/to/alignments_aa_fixed/"
#this is the output folder
pathAA<-"path/to/amino_acids/output/"

dat_HLA<-c()
for(f in 10:60){
  dat_HLA<-vroom(paste0(pathCalls, "hla_df_batch_", f, ".tsv.gz"), col_types = cols(.default = "c")) %>%
    mutate(across(everything(), gsub, pattern="HLA-", replacement="")) %>%
    mutate(across(!ID, str_extract, pattern=".*\\*[0-9]*:[0-9]*")) %>%
    rbind(dat_HLA,.)
}

#If you want to save that intermediate file and reload it after
#vroom_write(dat_HLA, paste0(pathAA,"hla_four_digit_aa_ready.tsv.gz"))
#dat_HLA<-vroom(paste0(pathAA,"hla_four_digit_aa_ready.tsv.gz"), col_types = cols(.default = "c"))

pos_update<-1

#start gene loop here
for(gene in c("A", "B", "C", "E", "F", "G",  
              "DMA", "DMB", "DOA", "DOB",
              "DPA1", "DPB1", "DQA1", "DQB1",
              "DRA", "DRB1", "DRB3", "DRB4", "DRB5")){

#for(gene in c("DPB1", "DQA1", "DQB1",
#              "DRA", "DRB1", "DRB3", "DRB4", "DRB5")){
  #print(gene)
  
  ukb_alleles<-c(unlist(dat_HLA[,paste0(gene,"_1")]), unlist(dat_HLA[,paste0(gene,"_2")])) %>% unique()
  
  #DRB1, DRB3, DRB4, DRB5
  if(gene %in% c("DRB1", "DRB3", "DRB4", "DRB5")){
    seq<-vroom(paste0(pathSeq,"DRB_four_fixed.tsv.gz"), col_types = cols(.default = "c")) %>%
      filter(allele %in% ukb_alleles)
    colnames(seq)<-gsub("DRB",gene,colnames(seq))
  } else {
    seq<-vroom(paste0(pathSeq,gene,"_four_fixed.tsv.gz"), col_types = cols(.default = "c")) %>%
      filter(allele %in% ukb_alleles)
  }
  
  aa1_null<-seq[base::match(unlist(dat_HLA[,paste0(gene,"_1")]),unlist(seq$allele)),-1] %>%
    rbind(seq[1,-1],.) %>%
    mutate(across(everything(), gsub, pattern="[\\*\\.]", replacement=NA))
    #mutate(across(everything(),~ifelse(.x=="*",.x[1],.x)))
  aa1_null<-aa1_null[-1,]
  aa1<-aa1_null
  colnames(aa1)<-paste0(colnames(aa1), "_chr1")
  aa2_null<-seq[base::match(unlist(dat_HLA[,paste0(gene,"_2")]),unlist(seq$allele)),-1] %>%
    rbind(seq[1,-1],.) %>%
    mutate(across(everything(), gsub, pattern="[\\*\\.]", replacement=NA))
    #mutate(across(everything(),~ifelse(.x=="*",.x[1],.x))) 
  aa2_null<-aa2_null[-1,]
  aa2<-aa2_null
  colnames(aa2)<-paste0(colnames(aa2), "_chr2")
  
  aa1_null<-aa1_null %>% distinct()
  aa2_null<-aa2_null %>% distinct()
  
  all_amino_acids<-rbind(aa1_null,aa2_null)
  
  aa_per_pos<-lapply(all_amino_acids, unique)
  
  for(i in 1:length(aa_per_pos)){
    aa_per_pos[[i]]<-aa_per_pos[[i]][!is.na(aa_per_pos[[i]])]
  }
  
  aa_per_pos_var_names<-names(aa_per_pos)[which(lengths(aa_per_pos)>1)]
  aa_per_pos_var<-aa_per_pos[aa_per_pos_var_names]
  
  #now remove NA and build the vcf per amino acid
  
  aa<-cbind(aa1, aa2)[order(c(seq_along(aa1), seq_along(aa2)))]
  aa<-aa %>%
    mutate(ID=dat_HLA$ID) %>%
    relocate(ID)
  
  vroom_write(aa, paste0("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/amino_acids/aa_calls/", gene, "_aa_calls.tsv.gz"))
  
  #now make the vector of variant names (i.e. the different amino acid possibilities)
  for(i in 1:length(aa_per_pos_var)){
    aa_per_pos_var[[i]]<-paste0(names(aa_per_pos_var)[i],"_",aa_per_pos_var[[i]])
  }
  
  var_names<-unlist(aa_per_pos_var) %>% as.vector()
  
  vcf_call<-matrix(NA, nrow=length(var_names), ncol=nrow(aa)) %>%
    as.data.frame(.)
  colnames(vcf_call)<-aa$ID
  rownames(vcf_call)<-var_names
  
  
  for(i in 1:nrow(vcf_call)){
    if(i %% 100 ==1){
      print(paste0("Gene ", gene, ": ",i,"/",nrow(vcf_call)))
    }
    tmp_row<-rownames(vcf_call)[i]
    aa_let<-substr(tmp_row, nchar(tmp_row), nchar(tmp_row))
    pos<-paste0(substr(tmp_row,1,nchar(tmp_row)-1), c("chr1", "chr2"))
    aa_dose<-aa[,pos] %>%
      mutate(across(everything(),
                    stri_replace_all_regex,
                    pattern=aa_let,
                    replacement="1",
                    vectorize=FALSE)) %>%
      mutate(across(everything(),~ifelse(.x!=1 | is.na(.x),0,.x))) %>%
      mutate(!!tmp_row:=paste(.[,1],.[,2],sep="/")) %>%
      dplyr::select(all_of(tmp_row))
    
    vcf_call[colnames(aa_dose)[1],]<-aa_dose[,1]
  }
  
  vcf<-data.frame(CHROM=rep(6, nrow(vcf_call)),
                  POS=c(pos_update:(pos_update+nrow(vcf_call)-1)),
                  ID=rownames(vcf_call),
                  REF=rep("C", nrow(vcf_call)),
                  ALT=rep("A", nrow(vcf_call)),
                  QUAL=rep(".", nrow(vcf_call)),
                  FILTER=rep(".", nrow(vcf_call)),
                  INFO=rep(".", nrow(vcf_call)),
                  FORMAT=rep("GT", nrow(vcf_call))) %>%
    cbind(.,vcf_call) %>%
    rename("#CHROM"=CHROM)
  
  vroom_write(vcf, paste0(pathAA,
                          gene,
                          "_amino_acids.tsv.gz"))
  
  pos_update<-pos_update+nrow(vcf_call)
}

