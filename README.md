# SARS-CoV-2_Pipeline
SARS-CoV-2 raw genomic sequence data processing and analysis pipeline

This contains the folders for the SARS-CoV-2 pipeline work. The folder 'scripts' has the following object files.

clean_mol_seq.py --- quality control step and filtering of ambiguous nucleotides

seq_to_covariate.py --- computing the variant constellation, synonymous nucleotide mutations, non-synonymous nucleotide mutations, and locations of ambiguous regions by protein

variant_dynamics.py --- runs seq_to_covariate.py and collect all variant constellations, single amino acid mutations, and counts by region and date and store in JSON

variant_analysis.py --- computes prevalence and growth rates by date and region for each variant constellation and stores as a panel dataframe

gisaid_metadata.py --- process GISAID metadata file and aggregate region-date counts for covariate sequences or single mutations

variant_scoring.py --- score coviarates or single matations based on the spatial-temporal epidemiological dynamics as well as based on experimental evidence for functional impact

Additional files:

test_pipeline.py --- test entire pipeline, running through all four files to generate region-time prevalence and growth rates data for each variant constellation and single mutation
test_gisaid.py --- used for testing gisaid_metadata.py object
test_scoring.py --- user for testing variant_scoring object
NC_045512.2.fasta --- Wuhan-Hu-1 reference sequences necessary for alignments and computing variants
genomes_test_run.fasta --- 20 genomes sampled from ViPR used for testing the pipeline
SARS2_April_500.fasta -- 500 genomes sampled from ViPR in the month of April 2022 used for testing
