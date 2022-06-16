Early Detection of SARS-CoV-2 Variants of Interest

This repository contains the code we make public for scoring emerging SARS-CoV-2 variants based on the algorithms described in 
Wallace ZS et al., "Early Detection of Emerging SARS-CoV-2 Variants of Interest for Experimental Evaluation", Frontiers in
Bioinformatics.

Dependencies: python >=3.x, pandas, edlib, Biopython, collections, typing, compress_json, datetime
Please use pip install to install these dependencies prior to running the software.

INTRODUCTION: 

Users can score and rank SARS-CoV-2 variants as the emerge for priotizing experimental evaluation based on newly updated data 
in the GISAID metadatafile or GenBank/BV-BRC FASTA files.  Rankings are based on a heuristic that combines epidemiological dynamics 
and predicted functional impacts.  Users can rank what we call variant constellations, or covariates for short, which are a constellation 
of amino acid mutations relative to the Wuhan-Hu-1 reference strain concatenated and represented as one string.  Users
can score and rank full proteome variant constellations, Spike protein variant constellations, or non-Spike variant constellations,
such as ranking variant constellations found on NSP3 or NSP5.  In addition, this software allows users to rank individual
amino acid mutations based on the same algorithm.

To run this pipeline on GISAID data, you must have a GISAID account set up and request downloads of the metadata file.  Donwnload
this file and input to the commandline.  Please keep in mind the GISAID metadata file is about 7.5 GB and growing.

To run this pipeline on GenBank/BV-BRC SARS-CoV-2 sequence data in FASTA format, first download batches from
https://www.viprbrc.org/brcDocs/datafiles/public_share/Corona/. As of now, the pipeline will only accept FASTA sequences from
those in this repository due to data parsing requirements. WARNING: The FASTA files are large and compute time is long due to
alignments and processing.  It is highly advised to have access to HPC resources for running on FASTA data.

RUNNING PIPELINE: 

cd path/to/local/dir

git clone https://github.com/zwallace425/SARS-CoV-2_Pipeline.git

cd scripts

To score and rank full proteome variant constellations:

	(1) GISAID: python main.py --gisaid [GISAID Metadata File] --analyze covariates
	(2) GenBank/BV-BRC: python main.py --ref NC_045512.2.fasta --query [FASTA File] --analyze covariates

To score and rank protein-specific variant constellations:

	(3) [GISAID or GenBank/BV-BRC Args] --protein Spike

To score and rank single amino acid mutations:

	(4) GISAID: python main.py --gisaid [GISAID Metadata File] --analyze mutations
	(5) GenBank/BV-BRC: python main.py --ref NC_045512.2.fasta --query [FASTA File] --analyze mutations

To score and rank protein-specific amino acid mutations:

	(6) [GISAID or GenBank/BV-BRC args] --protein Spike

The ranking results will be printed to the commandline and automatically deposited in the "results" directory with the current date
concatenated to the filename.  Please see the file "sample_composite_score.txt" in the "results" directory for an example output.

There are data files that get built and saved during execuation of the pipeline and stored in the "data" directory.  This includes
compressed JSON files for covariate/mutaton region-date counts, sequence isolate counts per region-date, variant constellations
mapping to accessions, and pickled dataframe files for storing variant/matution region-date counts and region-date isolate counts
in dataframe format.  Also initially included in the "data" directory are files necessary for scoring functional impact --- Spike_SFoCs.txt
and non-Spike_SFoCs.txt.  DO NOT remove those files from the "data" directory!

When running the pipeline on GenBank/BV-BRC FASTA data, the stored files in the "data" directory are more extensive.  On top of covariate/mutation
region-date counts and region-date isolate counts, the pipeline will store a file titled "bvbrc_variant_mappings.json.gz", a compressed JSON
file that maps each each nucleotide sequence to the full variant constellation as well as synonymous, non-synonymous mutations, and exact protein
regions where ambiguous nucleotides were identified.  There is also a file titled "bvbrc_SARS-CoV-2_Sequence_QC.json.gz", which is a compressed
JSON file mapping each GenBank accession to the nucleotide sequence to quality control results, such as the N-content, frameshift status, or if 
the record failed any metadata criteria like no region or date.

NOTES:

When running the pipeline to score covariates with the GISAID metadata file, compute time is no more than 15 minutes.
When running the pipeline on FASTA files that includes a batch download for a single month, compute time will be a few hours.

The following are key files in the "scripts" directory: 

clean_mol_seq.py --- quality control step and filtering of ambiguous nucleotides (courtesy of Christian Zmasek)

seq_to_covariate.py --- computing the variant constellation, synonymous nucleotide mutations, non-synonymous nucleotide mutations, and 
	locations of ambiguous regions by protein

variant_dynamics.py --- runs seq_to_covariate.py or is called from gisaid_metadata.py and collects all variant constellations, single 
	amino acid mutations, and counts by region and date and stores in JSON and pickled dataframes

variant_analysis.py --- computes prevalence and growth rates by date and region for each variant constellation or single mutation and stores
	as a panel dataframe

gisaid_metadata.py --- process GISAID metadata file and aggregate region-date counts for covariate sequences or single mutations

variant_scoring.py --- score coviarates or single matations based on the spatial-temporal epidemiological dynamics as well as based on experimental 
	evidence for functional impact

NC_045512.2.fasta --- Wuhan-Hu-1 reference sequences necessary for alignments and computing variants

reference_aa.py --- SARS-CoV-2 Wuhan-Hu-1 reference translated by gene

Additional files in the "testing" directory:

test_pipeline.py --- test entire pipeline, running through all four files to generate region-time prevalence and growth rates data for each variant 
	constellation and single mutation
test_gisaid.py --- used for testing gisaid_metadata.py object
test_scoring.py --- used  for testing variant_scoring object
test_variant_analysis.py --- used or testing variant_analysis.py
genomes_test_run.fasta --- 20 genomes sampled from ViPR used for testing the pipeline
SARS2_April_500.fasta -- 500 genomes sampled from ViPR in the month of April 2022 used for testing

FINAL REMARKS:

We are continually updating this repository to enhance code efficiency, algorithms, analysis options, and user experience.  Currently the plotting analysis
described in the manuscript is not available as that was originally designed for BV-BRC Emerging Variant Report spreadsheets generated from GISAID
data and used for consulting with the NIAID SAVE Consortium, but we cannot publically share those spreadsheets and code due to GISAID licensing agreements.

For any questions or bug catches please contact Zach Wallace at zwallace@jcvi.org or ask questions on GitHub for public discussion.


