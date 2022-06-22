#move to the "amino_acids" directory. 
cd /path/to/amino_acids/

#vcf header to insert
printf "##fileformat=VCFv4.3\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" | gzip > vcf_header.txt.gz

declare -a genes_hla=("A" "B" "C" "E" "F" "G" "DMA" "DMB" "DOA" "DOB" "DPA1" "DPB1" "DQA1" "DQB1" "DRA" "DRB1" "DRB3" "DRB4" "DRB5")
  
for g in "${genes_hla[@]}"; do
  zcat vcf_header.txt.gz ${g}_amino_acids.tsv.gz | bgzip > ${g}_amino_acids.vcf.gz
  tabix --csi -p vcf ${g}_amino_acids.vcf.gz
done

ls | grep _amino_acids.vcf.gz | | grep -v csi > list_aa_vcf.txt

bcftools concat --file-list list_aa_vcf.txt -Ou | \
  bcftools sort -Ou | \
  bcftools +fill-tags -Oz > full_hla_aa.vcf.gz
tabix -p vcf full_hla_aa.vcf.gz

plink --vcf full_hla_aa.vcf.gz --double-id --make-bed --out full_hla_aa
