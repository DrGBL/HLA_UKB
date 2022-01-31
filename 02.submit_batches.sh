conda activate dnanexus-env

#only need to run the following once, but needs to be re-ran everytime you change the workflow wdl
#it uses the dna-nexus wdl compiling java application, available here: https://github.com/dnanexus/dxCompiler/releases
java -jardxCompiler-2.5.0.jar compile workflow.wdl -project exome_full -folder /my_workflows/
