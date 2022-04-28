setwd("/project/richards/guillaume.butler-laporte/HLA/ukb_wes")

library(tidyverse)

user<-c("","_tianyuan", "_yiheng")
user_index<-2
folder<-42

files<-list.dirs(paste0("hla_calling_hla_hd_batch_", folder, "/")) %>%
  data.frame(sample_called=.) %>%
  filter(sample_called!=paste0("hla_calling_hla_hd_batch_", folder, "/")) %>%
  mutate(sample_called=substr(sample_called,
                              start=nchar(paste0("hla_calling_hla_hd_batch_", folder, "//"))+1,
                              stop=1000))

# files<-read_tsv("03.1.remove.sh", col_names=FALSE) %>%
#   rename(sample_called=X1) %>%
#   mutate(sample_called=substr(sample_called,
#                               start=nchar("dx rm /hla_calling_hla_hd/")+1,
#                               stop=nchar("dx rm /hla_calling_hla_hd/1054188_23153_0_0")))

expected<-read_tsv(paste0("batches", user[user_index], "/batch_", folder, ".tsv"))

missing<-expected %>%
  mutate(`batch ID`=substr(`batch ID`,1,7)) %>%
  filter(!(`batch ID` %in% substr(files$sample_called,1,7)))

write_tsv(missing, paste0("batches", user[user_index], "/batch_", folder, "_missing.tsv"))
         

nrow(missing)

