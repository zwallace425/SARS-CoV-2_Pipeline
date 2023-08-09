# Main program for running the entire SARS-CoV-2 variant tracking and analysis pipeline
# See the README.md for details

import sys
import time
import argparse
import warnings
import pandas as pd
import pickle
from datetime import date
from clean_mol_seq import CleanMolSeq as cms
from gisaid_metadata import GisaidMetadata as gm
from variant_dynamics import VariantDynamics as vd
from variant_analysis import VariantAnalysis as va
from variant_scoring import VariantScoring as vs
from variant_plots import VariantPlots as vp

warnings.filterwarnings("ignore")

BVBRC_PROTEINS = ['nsp1', 'nsp2', 'nsp3', 'nsp4', 'nsp5', 'nsp6', 'nsp7', 'nsp8', 'nsp9', 'nsp10', 'nsp11', 'nsp12', 'nsp13', 'nsp14', 'nsp15', 'nsp16',
'E', 'M', 'N', 'Spike', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'ORF10', 'ORF14']

GISAID_PROTEINS = ['NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'NSP11', 'NSP12', 'NSP13', 'NSP14', 'NSP15', 'NSP16',
'E', 'M', 'N', 'Spike', 'NS3', 'NS6', 'NS7a', 'NS7b', 'NS8']

WHO_NAMES = ['Omicron', 'Delta', 'Alpha', 'Beta', 'Mu', 'Iota', 'Epsilon', 'Gamma', 'Eta', 'Kappa', 'Zeta', 'Lambda', 'Theta']

# Requirements for running the pipeline
def ProgramUsage():
	print('\n')
	print("USAGE: The pipeline must either score or plot variants with the minimum commandline arguments:")
	print('\n')
	print("Scoring GISAID Data:")
	print("python main.py --gisaid [GISAID METADATA FILE] --analyze [covariants/mutations/lineages]")
	print("OR")
	print("Scoring GenBank/BV-BRC Data:")
	print("python main.py --ref NC_045512.2.fasta --query [FASTA FILE] --analyze [covariants/mutations]")
	print('\n')
	print("Plotting GISAID Data:")
	print("python main.py --gisaid [GISAID METADATA FILE] --plot [covariants/mutations/lineages")
	print('\n')
	print("NOTE: NO PLOTS WITH GENBANK/BV-BRC FASTA DATA. INCONVIENT EXPENSIVE COMPUTING.")
	print('\n')
	print("---------------------------------------------------------------------------------------------------")
	print('\n')
	print("OPTIONAL: The following are optional arguments following the required:")
	print('\n')
	print("[GISAID args] --WHO [WHO NAME]")
	print("[GISAID args] --PANGO [PANGO LINEAGE]")
	print("[GISAID or GenBank/BV-BRC args] --period [D/W/2W/M] --interval [>=2 AND <=6] --n_content [>0 AND <1] --seq_length [SEQUENCE LENGTH] --min_date [>2019-11-01] --max_date [YYYY-MM-DD] --protein [SARS-CoV-2 PROTEIN] --region [REGION] --covariants [VARIANTS FILE]")
	print('\n')
	print("NOTE: WHO or PANGO input only allowed with GISAID metadata.")
	print("NOTE: WHO and PANGO input not allowed together.")
	print("NOTE: WHO or PANGO input not allowed with '--analyze lineages' or '--plot lineages'. Unecessary analysis.")
	print("NOTE: --region not allowed with --analyze.  Scoring algorithms require global data.")
	print('\n')


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

# Used for opening the covariants file.  This file must have at least one column named 'Variant'
def open_covariants_file(file, protein): 
	try:
		covariants_file = pd.read_csv(file, sep = '\t', encoding = 'latin-1')
		if ('Variant' in covariants_file.columns):
			if ('Name' in covariants_file.columns):
				covariants_file['Name'] = covariants_file['Name'].str.encode('ascii', 'ignore').str.decode('ascii')
				covariants_file['Name'] = covariants_file['Name'].str.replace('Ê', '')
			covariants_file['Variant'] = covariants_file['Variant'].str.replace(' ', '')
			covariants_file['Variant'] = covariants_file['Variant'].str.replace('Ê', '')
			variants = list(covariants_file['Variant'])
			variants_prot = []
			for variant in variants:
				muts = variant.split(",")
				muts = [protein+"_"+mut for mut in muts]
				muts = ",".join(muts)
				variants_prot.append(muts)
			covariants_file['Protein+Variant'] = variants_prot
			return(covariants_file)
		else:
			sys.exit("Invalid covariants file format. Inputted covariants file must be tabular with one column labeled 'Variant'.  Optional 'Name' or 'Identifier' columns permitted. See ReadMe.md for more details.")
	except:
		raise Exception("Could not open covariants file")

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	
	# Required arguments:
	parser.add_argument('--gisaid', dest = 'gisaid', type = str)
	# Or (following two must be together)
	parser.add_argument('--ref', dest = 'ref', type = str)
	parser.add_argument('--query', dest = 'query', type = str)
	# And one of the following only
	parser.add_argument('--analyze', dest = 'analyze', type = str)
	parser.add_argument('--plot', dest = 'plot', type = str)

	# Optional arguments
	parser.add_argument('--period', dest = 'period', type = str)
	parser.add_argument('--interval', dest = 'interval', type = str)
	parser.add_argument('--n_content', dest = 'n_content', type = str)
	parser.add_argument('--seq_length', dest = 'seq_length', type = str)
	parser.add_argument('--min_date', dest = 'min_date', type = str)
	parser.add_argument('--max_date', dest = 'max_date', type = str)
	parser.add_argument('--protein', dest = 'protein', type = str)
	parser.add_argument('--WHO', dest = 'who', type = str)
	parser.add_argument('--PANGO', dest = 'pango', type = str)
	parser.add_argument('--region', dest = 'region', type = str)
	parser.add_argument('--covariants', dest = 'covariants', type = str)
	args = parser.parse_args()

	spike_sfoc = pd.read_csv("data/Spike_SFoCs.txt", sep = '\t')
	non_spike_sfoc = pd.read_csv("data/Non-Spike_SFoCs.txt", sep = '\t')

	# Default min date first of year, default max date first of current month
	today = date.today().strftime("%Y-%m-%d")
	min_date = date(date.today().year, 1, 1).strftime("%Y-%m-%d")
	max_date = date(date.today().year, date.today().month, 1).strftime("%Y-%m-%d")

	# Exit program if required commandline arguments are incorrect
	if (args.gisaid and (args.ref or args.query)):
		sys.exit(ProgramUsage())
	if (args.analyze and args.plot):
		sys.exit(ProgramUsage())
	if ((args.ref and (not args.query)) or (args.query and (not args.ref))):
		sys.exit(ProgramUsage())
	if (args.ref and (args.analyze == "lineages")):
		sys.exit(ProgramUsage())
	if (args.ref and (args.who or args.pango)):
		sys.exit(ProgramUsage())
	if args.analyze:
		if (args.analyze != "covariants" and args.analyze != "mutations" and args.analyze != "lineages"):
			sys.exit(ProgramUsage())
	if args.plot:
		if (args.plot != "covariants" and args.plot != "mutations" and args.plot != "lineages"):
			sys.exit(ProgramUsage())	
	if (args.who and args.pango):
		sys.exit(ProgramUsage())
	if (args.analyze == "lineages" and (args.who or args.pango)):
		sys.exit(ProgramUsage())
	if (args.plot == "lineages" and (args.who or args.pango)):
		sys.exit(ProgramUsage())	
	if (args.analyze and args.region):
		sys.exit(ProgramUsage())
	if (args.plot and args.ref):
		sys.exit(ProgramUsage())
	
	# Initialize all optional arguments
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
		# Different intialized interval for scoring vs plotting
		if args.analyze:
			interval = 4
		elif args.plot:
			interval = 6
	
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

	if args.max_date:
		max_date = args.max_date

	if args.who:
		who = args.who
		if who not in WHO_NAMES:
			print('\n')
			print("Invalid WHO name.")
			print('\n')
			print("Valid WHO names:", WHO_NAMES)
			print('\n')
			sys.exit()
	else:
		who = ""

	if args.pango:
		pango = args.pango
	else:
		pango = ""

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
		else:
			protein = args.protein
	else:
		protein = ""

	# Specific for plotting, specifies growth dynamics for whole world or individual regions
	if args.plot and args.region:
		region = args.region
		world_growth = False
	elif args.plot and (not args.region):
		region = 'World'
		world_growth = True
	else:
		world_growth = False

	# Inputted covariants files to score
	if args.covariants:
		if (not args.protein):
			print('\n')
			print("Must enter protein argument with covariants file.")
			print('\n')
			sys.exit()
		else:
			covariants_file = open_covariants_file(args.covariants, args.protein)
			covariants = list(covariants_file["Protein+Variant"])
	else:
		covariants = []




	# Running the pipline on GISAID metadata
	if args.gisaid:
		gisaid = open_gisaid_metadata(args.gisaid)
		gm_object = gm(gisaid, seq_length, min_date, max_date, who, pango)
		if (args.analyze == "covariants" or args.plot == "covariants"):
			if protein:
				prot_counts_df, prot_region_dates_df = gm_object.protein_counts(protein)
				prot_prevalence_df = va.analyze_dynamics(prot_counts_df, prot_region_dates_df, period, world_growth)
				if args.analyze:
					if protein == "Spike":
						scores = vs(prot_prevalence_df, interval).composite_score(spike_sfoc = spike_sfoc, covariates = covariants)
					else:
						scores = vs(prot_prevalence_df, interval).composite_score(non_spike_sfoc = non_spike_sfoc, covariates = covariants)
				elif args.plot:
					vp(prot_prevalence_df, interval, region, period).plot_covariates(protein, pango, who)
			else:
				cov_counts_df, cov_region_dates_df, cov_pango_df = gm_object.covariate_counts()
				cov_prevalence_df = va.analyze_dynamics(cov_counts_df, cov_region_dates_df, period, world_growth)
				if args.analyze:
					scores = vs(cov_prevalence_df, interval).composite_score(spike_sfoc = spike_sfoc, non_spike_sfoc = non_spike_sfoc)
					scores = cov_pango_df.merge(scores, on = 'Variant').sort_values(by = 'Composite Score', ascending = False).reset_index(drop=True)
				elif args.plot:
					vp(cov_prevalence_df, interval, region, period).plot_covariates(protein, pango, who)
			if args.analyze:
				if args.covariants:
					scores = scores.rename(columns = {"Variant": "Protein+Variant"})
					scores = covariants_file.merge(scores, on = "Protein+Variant", how = "outer").fillna(0).drop(columns = ["Protein+Variant"]).sort_values(by = ["Composite Score"], ascending = False).reset_index(drop=True)
				scores.to_csv("results/gisaid_composite_scores_"+period+"_"+today+".txt", sep = '\t', index = False)
				print(scores)
		elif (args.analyze == "mutations" or args.plot == "mutations"):
			aa_counts_df, aa_region_dates_df = gm_object.mutation_counts(protein)
			aa_prevalence_df = va.analyze_dynamics(aa_counts_df, aa_region_dates_df, period, world_growth)
			if args.analyze:
				if protein:
					prot_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score()
					prot_scores.to_csv("results/gisaid_"+protein+"_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
					print(prot_scores)
				else:
					spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(Spike = True)
					non_spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(nonSpike = True)
					spike_scores.to_csv("results/gisaid_spike_mutation_scores_"+period+"_"+today+".txt", sep = '\t', index = False)
					non_spike_scores.to_csv("results/gisaid_non_spike_mutation_scores_"+period+"_"+today+".txt", sep = '\t', index = False)
					print("Spike mutations ...")
					print(spike_scores)
					print('\n')
					print("Non-Spike mutations ...")
					print(non_spike_scores)
			elif args.plot:
					vp(aa_prevalence_df, interval, region, period).plot_mutations(protein, pango, who)
		elif (args.analyze == "lineages" or args.plot == "lineages"):
			if args.analyze:
				cov_counts_df, cov_region_dates_df, cov_pango_df = gm_object.covariate_counts()
				cov_prevalence_df = va.analyze_dynamics(cov_counts_df, cov_region_dates_df, period, world_growth)
				cov_prevalence_df = cov_pango_df.merge(cov_prevalence_df, on = 'Variant')
				scores = vs(cov_prevalence_df, interval).emerging_lineage_score()
				scores.to_csv("results/gisaid_lineage_scores"+period+"_"+today+".txt", sep = '\t', index = False)
				print(scores)
			elif args.plot:
				lineage_counts_df, lineage_region_dates_df = gm_object.lineage_counts()
				lineage_prevalence_df = va.analyze_dynamics(lineage_counts_df, lineage_region_dates_df, period, world_growth)
				vp(lineage_prevalence_df, interval, region, period).plot_lineages()



	# Running pipeline on GenBank/BV-BRC fasta data
	elif args.ref:
		
		start = time.time()
		clean_sequences = "data/clean_sequences.fasta"
		kept_acc = cms.clean_mol_seqs(args.query, clean_sequences, seq_length, 1 - n_content, True, 'data/SARS-CoV-2_Sequence_QC.json.gz')
		cov_counts_df, aa_counts_df, prot_counts_dict, region_dates_df = vd.spatial_temporal_prevalence(args.ref, clean_sequences, qc_storage_file = 'data/SARS-CoV-2_Sequence_QC.json.gz')

		if (args.analyze == "covariants"):
			if protein:
				cov_prevalence_df = va.analyze_protein_covariates(prot_counts_dict, region_dates_df, protein, period)
				if protein == "Spike":
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
			if (protein):
				if protein == "Spike":
					spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(Spike = True)
					spike_scores.to_csv("results/bvbrc_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
					print(spike_scores)
				else:
					non_spike_scores = vs(aa_prevalence_df, interval).mutation_prevalence_score(nonSpike = True)
					non_spike_scores = non_spike_scores[non_spike_scores['Variant'].str.contains(protein)]
					non_spike_scores.to_csv("results/bvbrc_"+protein+"_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
					print(non_spike_scores)
			else:
				spike_scores.to_csv("results/bvbrc_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
				non_spike_scores.to_csv("results/bvbrc_non_spike_mutation_scores"+period+"_"+today+".txt", sep = '\t', index = False)
				print(spike_scores)
				print('\n')
				print(non_spike_scores)

		end = time.time()
		print("Total run time:", end - start)
			







