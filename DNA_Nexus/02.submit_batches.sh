conda activate dnanexus-env

#only need to run the following once, but needs to be re-ran everytime you change the workflow wdl
#it uses the dna-nexus wdl compiling java application, available here: https://github.com/dnanexus/dxCompiler/releases
#note the name of my project here was "exome_full", change this to whatever is needed for you

java -jardxCompiler-2.5.0.jar compile workflow.wdl -project exome_full -folder /my_workflows/

#now move to the directory where you'll work (in my case, where the "batches" folder is located)
cd /project/richards/guillaume.butler-laporte/HLA/ukb_wes/

dx cd '/'

#example where I run batch 10
dx run my_workflows/hla_calling_wf \           #your workflow
  --batch-tsv batches/batch_10.tsv \    #the batch, this is the result of 01.create_batch.sh
  -istage-common.ref_genome="Homo_sapiens_assembly38.fasta.gz" \      #the next four files need to be on your project on dna-nexus, i.e. NOT on your local cluster
  -istage-common.ref_genome_index="Homo_sapiens_assembly38.fasta.gz.fai" \
  -istage-common.ref_genome_dict="Homo_sapiens_assembly38.dict" \
  -istage-common.ref_genome_gzi="Homo_sapiens_assembly38.fasta.gz.gzi" \
  --priority low \            #priority for the workers on dna-nexus, see here: https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/managing-job-priority
  --batch-folders \           #will create a different folder for every result in your batch
  --destination=/hla_calling_hla_hd_batch_10      #the destination of the results (you do not need to create it yourself, dna-nexus will create this folder automatically)

#now remove all fastq files (were created as byproduct of the workflow), if you want to remove them
dx find data --name sample.1.fastq.gz | awk '{print $6}' | sed 's/.*/dx rm &/' > 03.1.remove.sh
dx find data --name sample.2.fastq.gz | awk '{print $6}' | sed 's/.*/dx rm &/' > 03.2.remove.sh

sh 03.1.remove.sh
sh 03.2.remove.sh
