cd /your/local/folder/vcf/
  
for f in {10..60}; do
  zcat vcf_header.txt.gz pre_vcf/hla_four_digit_${f}.pre_vcf.tsv.gz | bgzip > hla_four_digit_${f}.vcf.gz
  tabix --csi -p vcf hla_four_digit_${f}.vcf.gz
  zcat vcf_header.txt.gz pre_vcf/hla_six_digit_${f}.pre_vcf.tsv.gz | bgzip > hla_six_digit_${f}.vcf.gz
  tabix --csi -p vcf hla_six_digit_${f}.vcf.gz;
  zcat vcf_header.txt.gz pre_vcf/hla_two_digit_${f}.pre_vcf.tsv.gz | bgzip > hla_two_digit_${f}.vcf.gz
  tabix --csi -p vcf hla_two_digit_${f}.vcf.gz;
done

bcftools merge --file-list list_four_vcf.txt --missing-to-ref -Ou | \
  bcftools +fill-tags -Oz > full_hla_four_digit.vcf.gz
tabix -p vcf full_hla_four_digit.vcf.gz

bcftools merge --file-list list_six_vcf.txt --missing-to-ref -Ou | \
  bcftools +fill-tags -Oz > full_hla_six_digit.vcf.gz
tabix --csi -p vcf full_hla_six_digit.vcf.gz

bcftools merge --file-list list_two_vcf.txt --missing-to-ref -Ou | \
  bcftools +fill-tags -Oz > full_hla_two_digit.vcf.gz
tabix --csi -p vcf full_hla_two_digit.vcf.gz

rm hla_*

plink --vcf full_hla_four_digit.vcf.gz --double-id --make-bed --out full_hla_four_digit
plink --vcf full_hla_six_digit.vcf.gz --double-id --make-bed --out full_hla_six_digit
plink --vcf full_hla_two_digit.vcf.gz --double-id --make-bed --out full_hla_two_digit

zcat full_hla_four_digit.vcf.gz | tail -n +19 | awk '{print $3}' > list_vars_four_digit.txt
zcat full_hla_six_digit.vcf.gz | tail -n +19 | awk '{print $3}' > list_vars_six_digit.txt
zcat full_hla_two_digit.vcf.gz | tail -n +19 | awk '{print $3}' > list_vars_two_digit.txt
