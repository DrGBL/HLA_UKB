
pathPlink=/path/to/vcf/
pathPlinkAA=/path/to/amino_acids/
pathOut=/path/out/
pathAnc=/path/to/ancestry/files/    #these are files containing FID and IID of each participants in each genetic ancestry

declare -a anc=("afr" "amr" "eas" "eur" "sas")

for a in "${anc[@]}"
do
  plink \
    --freq gz \
    --bfile ${pathPlink}full_hla_four_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_four_digit_${a}
  
  plink \
    --freq gz \
    --bfile ${pathPlink}full_hla_six_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_six_digit_${a}
    
    plink \
    --freq gz \
    --bfile ${pathPlinkAA}full_hla_aa \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_aa_${a}
    
  plink \
    --freq gz \
    --bfile ${pathPlink}full_hla_two_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_two_digit_${a}
    
  plink \
    --freq counts gz \
    --bfile ${pathPlink}full_hla_four_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_four_digit_${a}
  
  plink \
    --freq counts gz \
    --bfile ${pathPlink}full_hla_six_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_six_digit_${a}
    
    plink \
    --freq counts gz \
    --bfile ${pathPlinkAA}full_hla_aa \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_aa_${a}
    
  plink \
    --freq counts gz \
    --bfile ${pathPlink}full_hla_two_digit \
    --keep ${pathAnc}ukb.${a}IDsPCA.plink_final.txt \
    --out ${pathOut}freq_two_digit_${a}
done


plink \
  --freq gz \
  --bfile ${pathPlink}full_hla_four_digit \
  --out ${pathOut}freq_four_digit_all

plink \
  --freq gz \
  --bfile ${pathPlink}full_hla_six_digit \
  --out ${pathOut}freq_six_digit_all

plink \
  --freq gz \
  --bfile ${pathPlinkAA}full_hla_aa \
  --out ${pathOut}freq_aa_all

plink \
  --freq gz \
  --bfile ${pathPlink}full_hla_two_digit \
  --out ${pathOut}freq_two_digit_all
  
plink \
  --freq counts gz \
  --bfile ${pathPlink}full_hla_four_digit \
  --out ${pathOut}freq_four_digit_all

plink \
  --freq counts gz \
  --bfile ${pathPlink}full_hla_six_digit \
  --out ${pathOut}freq_six_digit_all

plink \
  --freq counts gz \
  --bfile ${pathPlinkAA}full_hla_aa \
  --out ${pathOut}freq_aa_all

plink \
  --freq counts gz \
  --bfile ${pathPlink}full_hla_two_digit \
  --out ${pathOut}freq_two_digit_all  
  
 
