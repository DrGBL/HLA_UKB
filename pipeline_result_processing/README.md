# pipeline_result_processing

This section provides code to retrieve HLA calls and coverage details for each participants. It also cleans the pipeline results in the format that will be provided to the UK Biobank.

Finally, it creates vcf files for two-digit, four-digit, six-digit, and amino acids. These are eventually made into plink files and used in the Regenie analysis later. The amino acids file was created from the four-digit calls, using the IMGT-HLA database v3450.

## Details

***01.calls_coverage_extraction.R***: this function extracts HLA-HD calls and coverage. For all HLA genes, exon 2 was used for coverage, except for DRB2 and DRB7, for which exon 2 is absent. Exon 3 was used in these two cases.

***02.build_pre_vcfs.R***: this reads calls for each gene, discards them if below a coverage threshold (10 is used by default), and writes a tsv file that will become a vcf in the next step. This creates one file for 2-digit (1-field), 4-digit (2-field), and 6-digit(3-field) precision. These were obtained by trimming the extra digits/fields and collapse them appropriately in one allele. For example, individuals with HLA-A*01:01:01 and HLA-A*01:01:02 at the 6-digit precision were collapsed to HLA-A*01:01 at the 4-digit precision.

***03.build_vcfs_plinks.sh***: this takes the pre-vcf files above, adds a header, and makes them into full vcf files. Note that the resulting POS, REF, and ALT columns are dummy variables. Also note that in this vcf, a homozygous individual for the HLA allele in the ID column will be listed as "1/1", and therefore 0/0 individuals do not have the listed allele (or had no calls for that allele with sufficient coverage). This also outputs the full list of alleles for each HLA resolutions.
