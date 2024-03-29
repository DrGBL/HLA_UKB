# DNA_Nexus

The following folder contains the code used to build the cram file batches, and run the HLA calling pipeline on DNA Nexus. This is done from the command line, on your local cluster.

Note that specific paths may need to be changed for your purposes. 

We assume that the user already has approved UK Biobank and DNA Nexus accounts.

## Details:

***hla_dnanexus.wdl***: this is the wdl script used in DNA Nexus. Note that due to copyrights, the HLA-HD docker image needs to be filled by the user before they can use this script (see "docker" folder in this git for more information). Also note that the user may need to change the DNA Nexus instance type used. Links to required files are provided in comments directly in the script.

***01.create_batch.sh***: this is used to create cram file batches to allow for batch calls (instead of manually starting the pipeline once for each of the UK Biobank participants). Note that the paths shown here may differ based on DNA Nexus data release, and we cannot guarantee that the location of whole-exome sequences cram files will be the same for all projects. Hence, paths may need to be modified based on the user's needs.

***02.submit_batches.sh***: this calls the pipeline for any given batch. Needs to be called for each of the 51 batches (10, 11, 12, ..., 60). Again, paths need to be adjusted. This is also where the wdl script above is compiled and fed to DNA Nexus.

The following functions are called once you have transferred the HLA calling results on your local cluster, and want to see which ones were missed/cancelled/errors:

***03.check_missing.R***: this is an R function provided to find participants for which the pipeline failed, and needs to be recalled. Using our setting (i.e. low priority job), an average of 1% of jobs were cancelled and needed to be restarted. This assumes that you have downloaded the batches obtained from step 02 above in your local cluster.

***04.check_missing_calls.R***: this is a second function to check for situation where the DNA Nexus job was killed before HLA calls were finished, but after a folder was created for a participant's results. This was also rare (included in the 1% figure above). These also need to be restarted.
