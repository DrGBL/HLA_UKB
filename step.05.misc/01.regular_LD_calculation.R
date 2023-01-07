pathPlink=/path/to/vcf/
pathOut=/path/to/ld_results/
pathAnc=/path/to/ancestry/files/    #these are files containing FID and IID of each participants in each genetic ancestry

declare -a anc=("afr" "amr" "eas" "eur" "sas")

for a in "${anc[@]}"
do
  plink \
    --r2 square gz \
    --bfile ${pathPlink}full_hla_four_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}ld_four_digit_${a}
  
  plink \
    --r2 square gz \
    --bfile ${pathPlink}full_hla_six_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}ld_six_digit_${a}
    
  plink \
    --r2 square gz \
    --bfile ${pathPlink}full_hla_two_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}ld_two_digit_${a}
done


