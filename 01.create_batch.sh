conda activate dnanexus-env

# dx login

cd /project/richards/guillaume.butler-laporte/HLA/ukb_wes/

dx cd '/Bulk/Exome sequences/Exome OQFE CRAM files/'

dx ls > folders_full_exomes.txt
sed -i 's/\///' folders_full_exomes.txt

# File cram_file
# File ref_genome
# File ref_genome_index
# File ref_genome_dict
# File ref_genome_gzi
# Int cores_hla
# Int min_read_length

#rm batches_prelim/*

while read folder; do
  dx generate_batch_inputs -icram_file='(.*).cram$' --path 'exome_full:Bulk/Exome sequences/Exome OQFE CRAM files/'${folder}'/' -o 'batches_prelim_tianyuan/folder_'${folder}
done < folders_full_exomes.txt

folder=30
head -n 1 batches_prelim/folder_${folder}.0000.tsv > batches/batch_${folder}.tsv
awk 'FNR>1' batches_prelim/folder_${folder}* >> batches/batch_${folder}.tsv
sed -i 's/cram_file/stage-common.cram_file/g' batches/batch_${folder}.tsv



dx cd /

#rm batches_prelim/*.tsv



#for tianyuan
folder=60
head -n 1 batches_prelim_tianyuan/folder_${folder}.0000.tsv > batches_tianyuan/batch_${folder}.tsv
awk 'FNR>1' batches_prelim_tianyuan/folder_${folder}* >> batches_tianyuan/batch_${folder}.tsv
sed -i 's/cram_file/stage-common.cram_file/g' batches_tianyuan/batch_${folder}.tsv


#for yiheng
head -n 1 batches_prelim_yiheng/folder_43.0000.tsv > batches_yiheng/batch_43.tsv
awk 'FNR>1' batches_prelim_yiheng/folder_43* >> batches_yiheng/batch_43.tsv
sed -i 's/cram_file/stage-common.cram_file/g' batches_yiheng/batch_43.tsv
