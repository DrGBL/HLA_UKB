pathIn=/path/to/gwas_ld_pred_ready/
pathLDPred=/path/to/ld_pred_output/

ls ${pathIn} | grep .bim | sed 's|.bim||g' | sed "s|^|${pathIn}|" > ${pathIn}list_plink_files.txt

plink \
  --merge-list ${pathIn}list_plink_files.txt \
  --make-bed \
  --threads 5 \
  --memory 30000 \
  --out ${pathIn}ukb_ld_pred_var

for pheno in "dm1" "psoriasis" "asthma" "celiac" "ms" "ra" "uc"
do
  for ana in "_with_hla" "_without_hla"
  do
    plink2 \
      --bfile ${pathIn}ukb_ld_pred_var \
      --score ${pathLDPred}ld_pred_${pheno}${ana}.tsv.gz 1 2 3 header cols=-scoreavgs,+scoresums \
      --threads 5 \
      --memory 20000 \
      --out ${pathLDPred}ld_pred_scored_${pheno}${ana}
  done
done
