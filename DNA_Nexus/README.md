# DNA_Nexus

The following folder contains the code used to build the cram file batches, and run the HLA calling pipeline on DNA Nexus.

Note that specific paths may need to be changed for your purposes. 

We assume that the user already has approved UK Biobank and DNA Nexus accounts.

## Details:

***hla_dnanexus.wdl***: this is the wdl script used in DNA Nexus. Note that due to copyrights, the HLA-HD docker image needs to be filled by the used before they can use this script (see "docker" folder in this git for more information). Also note that the user may need to change the DNA Nexus instance type used. Links to required files are provided in comments directly in the script.
