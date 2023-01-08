# this is sample code for the regenie burden test.
# it uses the regenie step 1 result from the previous section
# it should be ran separately for each ancestry
# please refer to the regenie documentation for file format

regenie \
  --step 2 \
  --bed full_hla_six_digit \ #plink bim/bed/fam file with 6 digit HLA alleles
  --covarFile covar_file.tsv \
  --phenoFile pheno_file.tsv \
  --keep ukb.${anc}IDsPCA.plink_final.txt \ # the ancestry file
  --firth 0.01 --approx \
  --anno-file anno_file.tsv.gz \ #from step 01
  --set-list set_list.tsv.gz \ #from step 01
  --mask-def mask_def.txt \ #available on this git
  --aaf-bins 0.999999 \
  --build-mask comphet \
  --pred ${pathOutputStep1}hla_ukb_wes_${anc}_pred.list \ #regenie step 1 file
  --bsize 400 \
  --minCaseCount 50 \
  --threads 10 \
  --bt \
  --minMAC 1 \
  --htp burden \
  --firth-se \
  --gz \
  --out ${pathOutputBurden}hla_ukb_wes_burden_${anc}
