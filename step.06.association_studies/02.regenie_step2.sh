# this is the code for regenie's step 2, shown for 6 digit HLA
# it uses the result in regenie's step 1
# this should be done separately for each ancestry
# please refer to regenie's documentation for file formats

pathGen=/path/to/vcf/

regenie \
  --step 2 \
  --bed ${pathGen}full_hla_six_digit \
  --covarFile covar_file.tsv \
  --phenoFile pheno_file.tsv \
  --keep ancestry_list.txt \
  --firth 0.01 --approx \
  --pred ${pathOutputStep1}hla_ukb_wes_${anc}_pred.list \
  --bsize 400 \
  --minCaseCount 50 \
  --threads 10 \
  --bt \
  --minMAC 1 \
  --htp six_digit \
  --firth-se \
  --gz \
  --out ${pathOutputStep2}hla_ukb_wes_four_digit_${anc}
