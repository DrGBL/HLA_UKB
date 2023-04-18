library(tidyverse)
library(vroom)

#where calls and QC folders are
setwd("/path/to/folder/")

coverage<-c()
for(folder in 10:60){
  coverage<-vroom(paste0("QC/hla_coverage_batch_", folder, ".tsv.gz"), col_types = cols(.default = "c")) %>% 
    bind_rows(coverage,.)
}

coverage_long<-coverage %>%
  pivot_longer(cols=-c(ID)) %>%
  mutate(value=as.numeric(value)) %>%
  rename(Genes=name,
         Coverage=value) %>%
  mutate(Genes=paste0("HLA-", Genes)) %>%
  filter(!is.na(Coverage))

cov_hist<-coverage_long %>%
  ggplot(aes(x=Coverage)) +
  geom_histogram(colour="black", fill="white", bins = 60) +
  geom_vline(xintercept=20,
             linetype="dashed",
             colour="red") +
  facet_wrap(~Genes, scales="free") +
  theme_bw()+
  xlab("")+
  ylab("")

alpha<-2
ggsave("supp_fig1_coverage_histogram.pdf", cov_hist, height = alpha*8.5, width = alpha*11)


