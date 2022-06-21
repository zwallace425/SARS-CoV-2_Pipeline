# Main program for running the entire SARS-CoV-2 variant tracking and analysis pipeline
# See the README.md for details

import sys
import argparse
import warnings
import pandas as pd
from datetime import date
from clean_mol_seq import CleanMolSeq as cms
from gisaid_metadata import GisaidMetadata as gm
from variant_dynamics import VariantDynamics as vd
from variant_analysis import VariantAnalysis as va
from variant_scoring import VariantScoring as vs

warnings.filterwarnings("ignore")

BVBRC_PROTEINS = ['nsp1', 'nsp2', 'nsp3', 'nsp4', 'nsp5', 'nsp6', 'nsp7', 'nsp8', 'nsp9', 'nsp10', 'nsp11', 'nsp12', 'nsp13', 'nsp14', 'nsp15', 'nsp16',
'E', 'M', 'N', 'Spike', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'ORF10', 'ORF14']

GISAID_PROTEINS = ['NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'NSP11', 'NSP12', 'NSP13', 'NSP14', 'NSP15', 'NSP16',
'E', 'M', 'N', 'Spike', 'NS3', 'NS6', 'NS7a', 'NS7b', 'NS8']

# Requirements for running the pipeline
def ProgramUsage():
	print('\n')
	print("USAGE: Commandline must have the minimum arguments:")
	print('\n')
	print("python main.py --gisaid [GISAID Metadata File] --analyze [covariates/mutations]")
	print("OR")
	print("python main.py --ref NC_045512.2.fasta --query [FASTA file] --analyze [covariates/mutations]")
	print('\n')
	print("OPTIONAL: The following are optional arguments following the required:")
	print('\n')
	print("[REQUIRED ARGS] --period [D/W/2W/M] --interval [>=2 and <=6] --n_content [>=0 and <=1] --seq_length [Sequence Length] --min_date [>2019-11-01] --protein [SARS-CoV-2 protein]")


# Open GISAID metadata file supplied to command line.
def open_gisaid_metadata(file):
	try:
		print("Opening GISAID metadata file ...")
		metadata = pd.read_csv(file, sep = '\t')
		print("Done opening file")
		print('\n')
		return(metadata)
	except:
		raise Exception("Could not open the GISAID metadata file")

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	
	# Required arguments:
	parser.add_argument('--gisaid', dest = 'gisaid', type = str)
	parser.add_argument('--analyze', dest = 'analyze', type = str)
	# Or (following two must be together)
	parser.add_argument('--ref', dest = 'ref', type = str)
	parser.add_argument('--query', dest = 'query', type = str)

	# Optional arguments
	parser.add_argument('--period', dest = 'period', type = str)
	parser.add_argument('--interval', dest = 'interval', type = str)
	parser.add_argument('--n_content', dest = 'n_content', type = str)
	parser.add_argument('--seq_length', dest = 'seq_length', type = str)
	parser.add_argument('--min_date', dest = 'min_date', type = str)
	parser.add_argument('--protein', dest = 'protein', type = str)
	args = parser.parse_args()

	spike_sfoc = pd.read_csv("data/Spike_SFoCs.txt", sep = '\t')
	non_spike_sfoc = pd.read_csv("data/Non-Spike_SFoCs.txt", sep = '\t')

	today = date.today()
	today = today.strftime("%Y-%m-%d")

	# Exit program if required commandline arguments are incorrect
	if (args.gisaid and (args.ref or args.query)):
		sys.exit(ProgramUsage())
	if ((args.ref and (not args.query)) or (args.query and (not args.ref))):
		sys.exit(ProgramUsage())
	if (args.analyze != "covariates" and args.analyze != "mutations"):
		sys.exit(ProgramUsage())
	
	# Initialize optional arguments
	if args.period:
		period = args.period
		if args.period not in ['D', 'W', '2W', 'M']:
			sys.exit(ProgramUsage())
	else:
		period = "M"

	if args.interval:
		interval = int(args.interval)
		if (interval > 6 or interval < 2):
			sys.exit(ProgramUsage())
	else:
		interval = 3
	
	if args.n_content:
		n_content = float(args.n_content)
		if (n_content <= 0 or n_content >= 1):
			sys.exit(ProgramUsage())
	else:
		n_content = 0.01
	
	if args.seq_length:
		seq_length = int(args.seq_length)
	else:
		seq_length = 29400
	
	if args.min_date:
		min_date = args.min_date
	else:
		min_date = "2022-01-01"

	if args.protein:
		if (args.gisaid and (args.protein not in GISAID_PROTEINS)):
			print('\n')
			print("Invalid GISAID SARS-CoV-2 protein entry.")
			print('\n')
			print("Valid GISAID protein entries:", GISAID_PROTEINS)
			print('\n')
			sys.exit()
		elif (args.ref and (args.protein not in BVBRC_PROTEINS)):
			print('\n')
			print("Invalid BV-BRC SARS-CoV-2 protein entry.")
			print('\n')
			print("Valid BV-BRC protein entries:", BVBRC_PROTEINS)
			print('\n')
			sys.exit()



	# Running the pipline on GISAID metadata or GenBank/BV-BRC fasta data
	if args.gisaid:
		gisaid = open_gisaid_metadata(args.gisaid)
		gm_object = gm(gisaid, seq_length, n_content, min_date)
		if (args.analyze == "covariates"):
			if args.protein:
				prot_counts_df, prot_region_dates_df = gm_object.protein_counts(args.protein)
				prot_prevalence_df = va.analyze_dynamics(prot_counts_df, prot_region_dates_df, period)
				if args.protein == "Spike":
					scores = vs(prot_prevalence_df, interval).composite_score(spike_sfoc = spike_sfoc)
				else:
					scores = vs(prot_prevalence_df, interval).composite_score(non_spike_sfoc = non_spike_sfoc)
			else:
				cov_counts_df, cov_region_dates_df = gm_object.covariate_counts()
				cov_prevalence_df = va.analyze_dynamics(cov_counts_df, cov_region_dates_df, period)
				scores = vs(cov_prevalence_df, interval).composite_score(spike_sfoc = spike_sfoc, non_spike_sfoc = non_spike_sfoc)
			scores.to_csv("results/gisaid_composite_scores_"+period+"_"+today+".txt", sep = '\t', index = False)
			print(scores)
		elif (args.analyze == "mutations"):
			aa_counts_df, aa_region_dates_df = gm_object.mutation_counts()
			aa_prevalence_df = va.analyze_dynamics(aa_counts_df, aa_region_dates_df, period)
			if (args.protein):
				if args.protein == "Spike":
					spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(Spike = True)
					spike_scores.to_csv("results/gisaid_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
					print(spike_scores)
				else:
					non_spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(nonSpike = True)
					non_spike_scores = non_spike_scores[non_spike_scores['Variant'].str.contains(args.protein)]
					non_spike_scores.to_csv("results/gisaid_"+args.protein+"_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
					print(non_spike_scores)
			else:
				spike_scores.to_csv("results/gisaid_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
				non_spike_scores.to_csv("results/gisaid_non_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
				print(spike_scores)
				print('\n')
				print(non_spike_scores)


	elif args.reference:
		
		clean_sequences = "data/clean_sequences.fasta"
		kept_acc = cms.clean_mol_seqs(args.query, clean_sequences, seq_length, 1 - n_content, True, 'data/SARS-CoV-2_Sequence_QC.json.gz')
		cov_counts_df, aa_counts_df, prot_counts_dict, region_dates_df = vd.spatial_temporal_prevalence(args.reference, clean_sequences, qc_storage_file = 'data/SARS-CoV-2_Sequence_QC.json.gz')

		if (args.analyze == "covariates"):
			if args.protein:
				cov_prevalence_df = va.analyze_protein_covariates(prot_counts_dict, region_dates_df, args.protein, period)
				if args.protein == "Spike":
					scores = vs(cov_prevalence_df, interval).composite_score(spike_sfoc = spike_sfoc)
				else:
					scores = vs(cov_prevalence_df, interval).composite_score(non_spike_sfoc = non_spike_sfoc)
			else:
				cov_prevalence_df = va.analyze_dynamics(cov_counts_df, region_dates_df, period)
				scores = vs(cov_prevalence_df, interval).composite_score(spike_sfoc = spike_sfoc, non_spike_sfoc = non_spike_sfoc)
			scores.to_csv("results/bvbrc_composite_scores_"+period+"_"+today+".txt", sep = '\t', index = False)
			print(scores)
		elif (args.analyze == "mutations"):
			aa_prevalence_df = va.analyze_dynamics(aa_counts_df, region_dates_df, period)
			if (args.protein):
				if args.protein == "Spike":
					spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(Spike = True)
					spike_scores.to_csv("results/bvbrc_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
					print(spike_scores)
				else:
					non_spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(nonSpike = True)
					non_spike_scores = non_spike_scores[non_spike_scores['Variant'].str.contains(args.protein)]
					non_spike_scores.to_csv("results/bvbrc_"+args.protein+"_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
					print(non_spike_scores)
			else:
				spike_scores.to_csv("results/bvbrc_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
				non_spike_scores.to_csv("results/bvbrc_non_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
				print(spike_scores)
				print('\n')
				print(non_spike_scores)
			







