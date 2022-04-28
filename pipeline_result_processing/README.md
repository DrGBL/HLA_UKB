# pipeline_result_processing

This section provides code to retrieve HLA calls and coverage details for each participants. It also cleans the pipeline results in the format that will be provided to the UK Biobank.

Finally, it creates vcf files for two-digit, four-digit, six-digit, and amino acids. These are eventually made into plink files and used in the Regenie analysis later. The amino acids file was created from the four-digit calls, using the IMGT-HLA database v3450.

## Details

***01.calls_coverage_extraction.R***: this function extracts HLA-HD calls and coverage.
