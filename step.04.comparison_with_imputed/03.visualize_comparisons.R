library(tidyverse)
library(vroom)
library(plotly)

#again go to the directory where the folders "comparisons_with_imputed", "calls" and "QC" are located
setwd("/your/local/folder")

comp<-c()
for(i in 10:60){
  comp<-vroom(paste0("comparisons_with_imputed/comparisons_imputed_called_folder_", i, ".tsv.gz")) %>%
    bind_rows(comp,.)
}
comp<-comp %>% mutate(across(where(is_character),as_factor))

sankey_df<-as.data.frame(matrix(data=0, nrow=22, ncol=10))
colnames(sankey_df)<-c("WES", 
                       "Gene", 
                       "Two allele matches", 
                       "One allele match, one allele mismatched",
                       "One allele match, one allele unimputed",
                       "Both alleles mismatched",
                       "One allele mismatched, one allele unimputed",
                       "Neither alleles imputed",
                       "One allele imputed",
                       "Both alleles imputed")
sankey_df$Gene<-rep(colnames(comp)[-1],each=2)
sankey_df$WES<-rep(c("Both alleles called", "Neither alleles called"), by=2)

link<-data.frame(table=c("Both called, both unimputed",
                         "One match, other unimputed",
                         "Both called, one unimputed, no match",
                         "Both uncalled and unimputed",
                         "One match, one mismatch imputation",
                         "Two matches",
                         "Both uncalled, one allele imputed",
                         "Both called, both imputed, no match",
                         "Both uncalled, both alleles imputed"),
                 call_res=c(1,
                            1,
                            1,
                            2,
                            1,
                            1,
                            2,
                            1,
                            2),
                 impute_res=c("Neither alleles imputed",
                              "One allele match, one allele unimputed",
                              "One allele mismatched, one allele unimputed",
                              "Neither alleles imputed",
                              "One allele match, one allele mismatched",
                              "Two allele matches",
                              "One allele imputed",
                              "Both alleles mismatched",
                              "Both alleles imputed"))

for(g in colnames(comp)[-1]){
  tmp<-comp %>% dplyr::select(paste0(g)) %>% table(.)
  for(j in 1:length(tmp)){
    sankey_df[which(sankey_df$Gene==g),
              link$impute_res[which(link$table==names(tmp)[j])]][link$call_res[which(link$table==names(tmp)[j])]] <- as.numeric(tmp[j])
  }
}

write_tsv(sankey_df, "comparisons_with_imputed/comparison_results.tsv")

nodes<- data.frame(name = c("Both Alleles Called", "Neither Alleles Called",
                            "Two allele matches", 
                            "One allele match, one allele mismatched",
                            "One allele match, one allele unimputed",
                            "Both alleles mismatched",
                            "One allele mismatched, one allele unimputed",
                            "Neither alleles imputed",
                            "One allele imputed",
                            "Both alleles imputed"))


for(g in colnames(comp)[-1]){
  fig <- plot_ly(
    type = "sankey",
    orientation = "h",
    
    node = list(
      label = c("Both Alleles Called", 
                "Neither Alleles Called",
                "Two allele matches", 
                "One allele match, one allele mismatched",
                "One allele match, one allele unimputed",
                "Both alleles mismatched",
                "One allele mismatched, one allele unimputed",
                "Neither alleles imputed",
                "One allele imputed",
                "Both alleles imputed"),
      #color = c("blue", "blue", "blue", "blue", "blue", "blue"),
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = rep(c(0,1),each=8),
      target = rep(c(2:9), 2),
      value =  c(as.numeric(sankey_df[which(sankey_df$Gene==g),][1,3:10]),
                 as.numeric(sankey_df[which(sankey_df$Gene==g),][2,3:10]))
    )
  )
  
  fig <- fig %>% layout(
    title = paste0("HLA-",g," calling vs imputation"),
    font = list(
      size = 10
    )
  )
  
  fig<-fig %>% plotly_build()
  
  save_image(fig, paste0("comparisons_with_imputed/", g,".pdf"))
  
}

