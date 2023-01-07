# the following is the code used for regenie's step 1
# it needs to be ran once per genetic ancestry
# please see the regenie documentation for the format of each files

regenie \
  --step 1 \
  --bed ukb_plink_file \ #the pruned plink bed/bim/fam file prefix, see manuscript for how this was pruned.
  --covarFile covar_file.tsv \
  --phenoFile pheno_file.tsv \
  --keep ancestry_list.txt \ #to keep the correct participants in the given ancestry
  --bt \
  --lowmem \
  --lowmem-prefix tmp_rg_${anc} \
  --minCaseCount 50 \
  --bsize 1000 \
  --catCovarList sex \
  --force-step1 \
  --threads 10 \
  --gz \
  --out ${pathOutput}hla_ukb_wes_${anc}
