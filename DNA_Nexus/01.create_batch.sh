#the goal of this code was to create batches to submit to dna nexus (the "batches/batch_${folder}.tsv" files at the end). 
#Column names are important and very tricky, see more here: https://documentation.dnanexus.com/user/running-apps-and-workflows/running-batch-jobs

#where you want batches to be recorded (on your local machine)
path_local_directory=/your/local/directory

#where the "Exome OQFE CRAM files" files are located (on your DNAnexus folder)
path_to_cram=/path/to/Exome OQFE CRAM files/

#I created a conda environment with the dna-nexus python package: https://documentation.dnanexus.com/downloads
conda activate dnanexus-env

dx login #enter your username and password, and then specify you dna nexus project

#go to working directory and make the right folders
cd ${path_local_directory}

mkdir -p batches
mkdir -p batches_prelim

#now go to the right location in dna nexus
dx cd ${path_to_cram}  '/Bulk/Exome sequences/Exome OQFE CRAM files/'

#list the folders (10 to 60)
dx ls > folders_full_exomes.txt
sed -i 's/\///' folders_full_exomes.txt

#loop through the folders to extract participant IDs
while read folder; do
  dx generate_batch_inputs -icram_file='(.*).cram$' --path ${path_to_cram}${folder}'/' -o 'batches_prelim/folder_'${folder}
done < folders_full_exomes.txt

while read folder; do
  head -n 1 batches_prelim/folder_${folder}.0000.tsv > batches/batch_${folder}.tsv
  awk 'FNR>1' batches_prelim/folder_${folder}* >> batches/batch_${folder}.tsv
  sed -i 's/cram_file/stage-common.cram_file/g' batches/batch_${folder}.tsv
done < folders_full_exomes.txt

