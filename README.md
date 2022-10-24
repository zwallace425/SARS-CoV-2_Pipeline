# Early Detection of SARS-CoV-2 Variants of Interest

This repository contains the code we make public for scoring and visualizing emerging SARS-CoV-2 variants based on the algorithms described in 
Wallace ZS, Davis J, Niewiadomska AM, Olson RD, Shukla M, Stevens R, Zhang Y, Zmasek CM and Scheuermann RH (2022) Early detection of emerging 
SARS-CoV-2 variants of interest for experimental evaluation. Front. Bioinform. 2:1020189. doi: 10.3389/fbinf.2022.1020189,
https://www.frontiersin.org/articles/10.3389/fbinf.2022.1020189/full

DEPENDENCIES: 

	python >=3.x, pandas, edlib, Biopython, collections, typing, compress_json, datetime

Please use pip install to install these dependencies prior to running the software.

## Introduction: 

Users can score and rank SARS-CoV-2 variants as the emerge for priotizing experimental evaluation based on newly updated data 
in the GISAID metadatafile or GenBank/BV-BRC FASTA files.  Rankings are based on a heuristic that combines epidemiological dynamics 
and predicted functional impacts.  Users can rank what we call variant constellations, or covariants for short, which are a constellation 
of amino acid mutations relative to the Wuhan-Hu-1 reference strain concatenated and represented as one string.  Users
can score and rank full proteome variant constellations, Spike protein variant constellations, or non-Spike variant constellations,
such as ranking variant constellations found on NSP3 or NSP5.  In addition, this software allows users to rank individual
amino acid mutations based on the same algorithm.  Lastly, this software allows for visualization of covariant, single mutation, or 
PANGO lineage growth dynamics.

To run this pipeline on GISAID data (recommended), you must have a GISAID account set up and request downloads of the metadata file.  Donwnload
this file and input to the commandline.  Please keep in mind the GISAID metadata file is around 10 GB and growing, so minimize downloads
and store appropriately.

To run this pipeline on GenBank/BV-BRC SARS-CoV-2 sequence data in FASTA format, first download batches from
https://www.viprbrc.org/brcDocs/datafiles/public_share/Corona/. As of now, the pipeline will only accept FASTA sequences from
those in this repository due to data parsing requirements. WARNING: The FASTA files are large and compute time is long due to
alignments and processing.  It is highly advised to have access to HPC resources for running on FASTA data.

## Running the Pipeline:

Run the following commands to get the pipeline working in your local dir.  Note, every analysis begins with "python main.py"

cd path/to/local/dir

git clone https://github.com/zwallace425/SARS-CoV-2_Pipeline.git

cd SARS-CoV-2_Pipeline/scripts

### (1) Scoring covariants/mutations/lineages:

To score and rank full proteome variant constellations with the Composite Score (minimum commandline requirments):

	(1) GISAID: python main.py --gisaid [GISAID Metadata File] --analyze covariants
	(2) GenBank/BV-BRC: python main.py --ref NC_045512.2.fasta --query [FASTA File] --analyze covariants

To score and rank protein-specific variant constellations with the Composite Score:

	(3) [GISAID or GenBank/BV-BRC Args] --protein [SARS-CoV-2 Protein]

To score and rank single amino acid mutations with the Mutation Prevalence Score:

	(4) GISAID: python main.py --gisaid [GISAID Metadata File] --analyze mutations
	(5) GenBank/BV-BRC: python main.py --ref NC_045512.2.fasta --query [FASTA File] --analyze mutations

To score and rank protein-specific amino acid mutations with the Mutation Prevalence Score:

	(6) [GISAID or GenBank/BV-BRC args] --protein [SARS-CoV-2 Protein]

To score and rank variant constellations or single amino acid mutations within a specific PANGO Lineage or WHO Name:

	(7) python main.py --gisaid [GISAID Metadata File] --analyze covariants --WHO [WHO Name]
	(8) python main.py --gisaid [GISAID Metadata File] --analyze covariants --PANGO [PANGO Lineage]
	(9) python main.py --gisaid [GISAID Metadata File] --analyze mutations --WHO [WHO Name]
	(10) python main.py --gisaid [GISAID Metadata File] --analyze mutations --PANGO [PANGO Lineage]
	(11) [ANY OF THE ABOVE] --protein [SARS-CoV-2 Protein]

	NOTE: Scoring of covariants/mutations specific to a PANGO Lineage/WHO Name not permitted with FASTA data

To score and rank PANGO Lineages with the Emerging Lineage Score:
	
	(12) python main.py --gisaid [GISAID Metadata File] --analyze lineages

	NOTE: Scoring PANGO Lineages not permitted with FASTA data. Additionally, any '--protein' input will be ignored

### (2) Visualizing growth of covariants/mutations/lineages:

NOTE: The plots will only dispaly the 15 covariants/mutations/lineages with the highest prevlance within a six month window.
NOTE: Plotting only allowed with GISAID Metatdata file

To plot growth of full proteome covariants or mutations:

	(13) Covariants: python main.py --gisaid [GISAID Metadata File] --plot covariants
	(14) Single Mutations: pyton main.py --gisaid [GISAID Metadata File] --plot mutations

To plot growth of protein specific covariants or mutations:

	(15) Covariants: python main.py --gisaid [GISAID Metadata File] --plot covariants --protein [SARS-CoV-2 Protein]
	(16) Single Mutations: python main.py --gisaid [GISAID Metadata File] --plot mutations --protein [SARS-CoV-2 Protein]

To plot growth of PANGO Lineages:

	(17) python main.py --gisaid [GISAID Metadata File] --plot lineages

To plot growth of covariants, mutations, or lineages specific to a region:
	
	(18) [Commands 13 - 16] --region [REGION]


The following are optional arguments that can be inputted with any of the analysis options:
	
	[REQUIRED ARGS] --period [D/W/2W/M] --interval [>=2 and <=6] --n_content [>0 and <1] --seq_length [NT Sequence Length] --min_date [>2019-11-01] --max_date [YYYY-MM-DD]

If a --max_date argument is supplied, then the pipeline will compute a Composite Score or plot up to that date instead of up to the current date,
for the purpose of studying results from the past and comparing to the present.

The ranking results will be printed to the commandline and automatically deposited in the "results" directory with the current date
concatenated to the filename.  Please see the file "sample_composite_score.txt" in the "results" directory for an example output.

The plotting results will be will be displayed to the user's screen and can be manually saved.

There are data files that get built and saved during execuation of the pipeline and stored in the "data" directory.  This includes
compressed JSON files for covariant/mutaton region-date counts, sequence isolate counts per region-date, variant constellations
mapping to accessions, and pickled dataframe files for storing variant/matution region-date counts and region-date isolate counts
in dataframe format.  Also initially included in the "data" directory are files necessary for scoring functional impact --- Spike_SFoCs.txt
and non-Spike_SFoCs.txt.  DO NOT remove those files from the "data" directory!

When running the pipeline on GenBank/BV-BRC FASTA data, the stored files in the "data" directory are more extensive.  On top of covariant/mutation
region-date counts and region-date isolate counts, the pipeline will store a file titled "bvbrc_variant_mappings.json.gz", a compressed JSON
file that maps each each nucleotide sequence to the full variant constellation as well as synonymous, non-synonymous mutations, and exact protein
regions where ambiguous nucleotides were identified.  There is also a file titled "bvbrc_SARS-CoV-2_Sequence_QC.json.gz", which is a compressed
JSON file mapping each GenBank accession to the nucleotide sequence to quality control results, such as the N-content, frameshift status, or if 
the record failed any metadata criteria like no region or date.

## Files:

When running the pipeline to score covariants with the GISAID metadata file, compute time is no more than 15 minutes.
When running the pipeline on FASTA files that includes a batch download for a single month, compute time will be a few hours.

The following are key files in the "scripts" directory: 

clean_mol_seq.py --- quality control step and filtering of ambiguous nucleotides (courtesy of Christian Zmasek)

seq_to_covariate.py --- computing the variant constellation, synonymous nucleotide mutations, non-synonymous nucleotide mutations, and 
	locations of ambiguous regions by protein, (Note, covariate is another term for covariant)

variant_dynamics.py --- runs seq_to_covariant.py or is called from gisaid_metadata.py and collects all variant constellations, single 
	amino acid mutations, and counts by region and date and stores in JSON and pickled dataframes

variant_analysis.py --- computes prevalence and growth rates by date and region for each variant constellation or single mutation and stores
	as a panel dataframe

gisaid_metadata.py --- process GISAID metadata file and aggregate region-date counts for covariant sequences, single mutations, or PANGO lineages

variant_scoring.py --- score coviarates, single matations, or PANGO lineages based on the spatial-temporal epidemiological dynamics as well as 
	based on experimental evidence for functional impact

variant_plots.py --- plot growth over time of covariants, single mutations, or PANGO lineages based on sptail-temporal epidemiological dynamics

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

NOTE: Throughout the python scripts, often times we refer to "covariant" as "covariate". They mean the same thing, variant constellation.

## Supplementary Material:

In the "materials" directory we offer a detailed version of the algorithms utilized in this pipeline in the form of a pseudocode pdf file.  This pseudocode 
pretty much covers the entire pipeline all the way from data processing, alignments, capturing dynamics, to scoring heuristics and is labeled as
"Early_Detection_Pseudocode.pdf".

In addition to pseudocode, the materials directory offers the Spike Sequence Features of Concern and Non-Spike Sequence Features of Concern data as an
additional location other than the "scripts/data" directory, the latter being meant for the pipeline.  These data are labeled as Spike_SFoCs.txt and 
non-Spike_SFoCs.txt

Finally, we offer the Supplementary Materials file from the manuscript labeled as "Manuscript_Supplementary_Materials.pdf".

## Final Remarks:

We are continually updating this repository to enhance code efficiency, algorithms, analysis options, and user experience.  In particular, if new knowledge
for updating the Sequence Features of Concern list is available, such as sites notable for resistance against any therapuetics, we will update that table
for Functional Impact Scoring.

For any questions or bug catches please contact Zach Wallace at zwallace@jcvi.org or ask questions on GitHub for public discussion.


