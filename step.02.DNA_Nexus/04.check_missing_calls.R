#go to the directory where the folders "hla_calling_hla_hd_batch_X" are located
setwd("/your/local/folder")

library(tidyverse)

folder<-42 #the folder you're testing

files<-list.dirs(paste0("hla_calling_hla_hd_batch_", folder)) %>%
  data.frame(files=.) %>%
  mutate(files=paste0(files, "/hla_res.tar.gz")) %>%
  filter(files!=paste0("hla_calling_hla_hd_batch_", folder, "/hla_res.tar.gz")) %>%
  filter(!file.exists(files)) %>%
  mutate(files=substr(files, 
                      nchar(paste0("hla_calling_hla_hd_batch_",
                                   folder,
                                   "/"))+1, 1000000)) %>%
  mutate(files=substr(files, 1, nchar("2100046"))) %>%
  rename(`batch ID`=files)

missing<-read_tsv(paste0("batches/batch_", folder, ".tsv")) %>%
  mutate(`batch ID`=substr(`batch ID`, 1, nchar("2100046"))) %>%
  merge(files,.,all = FALSE)

write_tsv(missing, paste0("batches/batch_", folder, "_missing_calls.tsv"))
