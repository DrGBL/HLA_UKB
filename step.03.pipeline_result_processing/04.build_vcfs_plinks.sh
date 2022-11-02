#again go to the directory where the folders "vcf" is located
cd /your/local/folder
 
#vcf header to insert
printf "##fileformat=VCFv4.3\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" | gzip > vcf_header.txt.gz 

for f in {10..60}; do
  zcat vcf_header.txt.gz vcf/pre_vcf/hla_four_digit_${f}.pre_vcf.tsv.gz | bgzip > vcf/hla_four_digit_${f}.vcf.gz
  tabix --csi -p vcf vcf/hla_four_digit_${f}.vcf.gz
  zcat vcf_header.txt.gz vcf/pre_vcf/hla_six_digit_${f}.pre_vcf.tsv.gz | bgzip > vcf/hla_six_digit_${f}.vcf.gz
  tabix --csi -p vcf vcf/hla_six_digit_${f}.vcf.gz;
  zcat vcf_header.txt.gz vcf/pre_vcf/hla_two_digit_${f}.pre_vcf.tsv.gz | bgzip > vcf/hla_two_digit_${f}.vcf.gz
  tabix --csi -p vcf vcf/hla_two_digit_${f}.vcf.gz;
done

ls vcf/pre_vcf | grep four | sed 's|pre_vcf|vcf|g' > vcf/list_four_vcf.txt
ls vcf/pre_vcf | grep two | sed 's|pre_vcf|vcf|g' > vcf/list_two_vcf.txt
ls vcf/pre_vcf | grep six | sed 's|pre_vcf|vcf|g' > vcf/list_six_vcf.txt

bcftools merge --file-list vcf/list_four_vcf.txt --missing-to-ref -Ou | \
  bcftools +fill-tags -Oz > vcf/full_hla_four_digit.vcf.gz
tabix -p vcf vcf/full_hla_four_digit.vcf.gz

bcftools merge --file-list vcf/list_six_vcf.txt --missing-to-ref -Ou | \
  bcftools +fill-tags -Oz > vcf/full_hla_six_digit.vcf.gz
tabix --csi -p vcf vcf/full_hla_six_digit.vcf.gz

bcftools merge --file-list vcf/list_two_vcf.txt --missing-to-ref -Ou | \
  bcftools +fill-tags -Oz > vcf/full_hla_two_digit.vcf.gz
tabix --csi -p vcf vcf/full_hla_two_digit.vcf.gz

rm vcf/hla_*

plink --vcf vcf/full_hla_four_digit.vcf.gz --double-id --make-bed --out vcf/full_hla_four_digit
plink --vcf vcf/full_hla_six_digit.vcf.gz --double-id --make-bed --out vcf/full_hla_six_digit
plink --vcf vcf/full_hla_two_digit.vcf.gz --double-id --make-bed --out vcf/full_hla_two_digit

zcat vcf/full_hla_four_digit.vcf.gz | tail -n +19 | awk '{print $3}' > vcf/list_vars_four_digit.txt
zcat vcf/full_hla_six_digit.vcf.gz | tail -n +19 | awk '{print $3}' > vcf/list_vars_six_digit.txt
zcat vcf/full_hla_two_digit.vcf.gz | tail -n +19 | awk '{print $3}' > vcf/list_vars_two_digit.txt
