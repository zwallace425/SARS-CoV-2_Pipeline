# SARS-CoV-2_Pipeline
SARS-CoV-2 raw genomic sequence data processing and analysis pipeline

This contains the folders for the SARS-CoV-2 pipeline work. The folder 'scripts' has the following object files.

clean_mol_seq.py --- quality control step and filtering of ambiguous nucleotides
seq_to_covariate.py --- computing the variant constellation, synonymous nucleotide mutations, non-synonymous nucleotide mutations, and locations of ambiguous regions by protein
variant_dynamics.py --- runs seq_to_covariate.py and collect all variant constellations, single amino acid mutations, and counts by region and date and store in JSON
variant_analysis.py --- computes prevalence and growth rates by date and region for each variant constellation and stores as a panel dataframe

Additional files:

test_pipeline.py --- test entire pipeline, running through all four files to generate region-time prevalence and growth rates data for each variant constellation and single mutation
NC_045512.2.fasta --- Wuhan-Hu-1 reference sequences necessary for alignments and computing variants
genomes_test_run.fasta --- 20 genomes sample from ViPR used for testing the pipeline
