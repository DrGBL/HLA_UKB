setwd("/your/local/folder")

library(tidyverse)

folder<-42 #the batch you're checking

files<-list.dirs(paste0("hla_calling_hla_hd_batch_", folder, "/")) %>%
  data.frame(sample_called=.) %>%
  filter(sample_called!=paste0("hla_calling_hla_hd_batch_", folder, "/")) %>%
  mutate(sample_called=substr(sample_called,
                              start=nchar(paste0("hla_calling_hla_hd_batch_", folder, "//"))+1,
                              stop=1000))

expected<-read_tsv(paste0("batches/batch_", folder, ".tsv"))

missing<-expected %>%
  mutate(`batch ID`=substr(`batch ID`,1,7)) %>%
  filter(!(`batch ID` %in% substr(files$sample_called,1,7)))

write_tsv(missing, paste0("batches/batch_", folder, "_missing.tsv"))

