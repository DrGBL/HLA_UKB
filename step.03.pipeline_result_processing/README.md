# pipeline_result_processing

This section provides code to retrieve HLA calls and coverage details for each participants. It also cleans the pipeline results in the format that will be provided to the UK Biobank.

Finally, it creates vcf files for two-digit, four-digit, six-digit, and amino acids. These are eventually made into plink files and used in the Regenie analysis later. The amino acids file was created from the four-digit calls, using the IMGT-HLA database v3450.

These are done on your local cluster, once you have downloaded the files resulting from the HLA calling.

## Details

***01.calls_coverage_extraction.R***: this function extracts HLA-HD calls and coverage. For all HLA genes, exon 2 was used for coverage, except for DRB2 and DRB7, for which exon 2 is absent. Exon 3 was used in these two cases.

***02.build_pre_vcfs.R***: this reads calls for each gene, discards them if below a coverage threshold (10 is used by default), and writes a tsv file that will become a vcf in the next step. This creates one file for 2-digit (1-field), 4-digit (2-field), and 6-digit(3-field) precision. These were obtained by trimming the extra digits/fields and collapse them appropriately in one allele. For example, individuals with HLA-A*01:01:01 and HLA-A*01:01:02 at the 6-digit precision were collapsed to HLA-A*01:01 at the 4-digit precision.

***03.build_vcfs_plinks.sh***: this takes the pre-vcf files above, adds a header, and makes them into full vcf files. It also produces plink files from those vcf files. Note that the resulting POS, REF, and ALT columns are dummy variables. Also note that in this vcf, a homozygous individual for the HLA allele in the ID column will be listed as "1/1", and therefore 0/0 individuals do not have the listed allele (or had no calls for that allele with sufficient coverage). This also outputs the full list of alleles for each HLA resolutions.

***04.munge_IMGTHLAv3450_prot.sh***: this processes the IMGT-HLA v3450 protein alignment files.

***05.munge_IMGTHLAv3450_prot_part2.R***: part 2 of the IMGT-HLA processing of protein alignment files. The end result for each gene is a data frame where each column is a different amino acid position, and each row is a different HLA allele for that gene.

***06.prep_aa_vcf.R***: this parses through the HLA allele calls for each participant in the UK Biobank, and assigns them their corresponding amino acid sequence, then builds a data frame that will be made into a vcf file in the next step.

***07.build_vcf_plink_aa.sh***: this builds the vcf and the plink files for the amino acid association studies.

## Note for the amino acid
Amino acids position are numbered according to the reference IMGT-HLA sequence. Therefore, there can be protein "indels", which cannot be directly numbered without disrupting the rest of the protein alignment numbering. Hence, we used the following convention:
- Amino acid variant IDs are named gene_position e.g. for gene A, position 123: A_123.
- For indels, we add "indel" prior to the position, but still increase position by 1 e.g. for an insertion after position 123 in gene A: A_indel124.
- We continue to do this, until we reach the next non-indel position. For example, for a sequence of five amino acid insertion after position 123 in gene A, we would write A_123, A_indel124, A_indel125, A_indel126, A_indel127, A_indel128, A_124.

This way the alignment numbering is preserved when there are no indels (allowing for proper association analyses), but indels are also clearly numbered.
