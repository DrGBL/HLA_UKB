#this code extracts the snps used by LDpred and builds a plink bim/bed/fam file for use in the next step
#it is ran on each chromosome

pathUKB=/path/to/ukb/imputed.v3/bgen/ 
pathSample=/path/to/bgen/sample_file.sample
pathOut=/path/to/gwas_ld_pred_ready/
pathSnps=/path/to/snps_list.tsv

plink2 \
  --bgen ${pathUKB}ukb_imp_chr${x}_v3.bgen ref-first \
  --sample ${pathSample} \
  --extract ${pathSnps} \
  --threads 5 \
  --memory 60000 \
  --make-bed \
  --out ${pathOut}ukb_ld_pred_var_chr${x}
