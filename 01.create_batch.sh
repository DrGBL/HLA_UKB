#the goal of this code was to create batches to submit to dna nexus (the "batches/batch_${folder}.tsv" files at the end). 
#Column names are important and very tricky, see more here: https://documentation.dnanexus.com/user/running-apps-and-workflows/running-batch-jobs

#I created a conda environment with the dna-nexus python package: https://documentation.dnanexus.com/downloads
conda activate dnanexus-env

dx login #enter your username and password

#below I'm showing what I did for the HLA work, adapt as needed
#it assumes that there's a folder called "batches" and another one named "batches_prelim" in the working directory

cd /project/richards/guillaume.butler-laporte/HLA/ukb_wes/

dx cd '/Bulk/Exome sequences/Exome OQFE CRAM files/'

dx ls > folders_full_exomes.txt
sed -i 's/\///' folders_full_exomes.txt

while read folder; do
  dx generate_batch_inputs -icram_file='(.*).cram$' --path 'exome_full:Bulk/Exome sequences/Exome OQFE CRAM files/'${folder}'/' -o 'batches_prelim/folder_'${folder}
done < folders_full_exomes.txt

read folder; do
  head -n 1 batches_prelim/folder_${folder}.0000.tsv > batches/batch_${folder}.tsv
  awk 'FNR>1' batches_prelim/folder_${folder}* >> batches/batch_${folder}.tsv
  sed -i 's/cram_file/stage-common.cram_file/g' batches/batch_${folder}.tsv
done < folders_full_exomes.txt

dx cd /
