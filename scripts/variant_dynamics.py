# Object for capturing epidemiological dynamics of SARS-CoV-2 variants based on genomic sequencing data from BV-BRC.
# The general return from the object is a dataframe of variant counts per region and date and is passed to
# the VariantAnalysis object for growth, prevalence, and jerk computations per date bin.

# NOTE: To maintin consistent FASTA format, this object only accepts a query FASTA file in the format
# downloaded from BV-BRC/ViPR. This object WILL NOT work with a GISAID FASTA file. Accepted FASTA files
# can be downloaded at https://www.viprbrc.org/brcDocs/datafiles/public_share/Corona/

import sys
import re
import pickle
import compress_json
import pandas as pd
from datetime import datetime
from Bio import SeqIO
import seq_to_covariate as stc

class VariantDynamics(object):

	# This function converts the dictionary containing the isolate counts per region and date into a dataframe
	# format that can be easily grouped into date intervals.  This can be used to get isolate counts for a
	# specific region over certain dates, or it be used to get the isolate counts for a specific date within
	# the entire world.  Primarly, this data is needed to complement the aa_prevalence data, as there's no way
	# of knowing the isolates counts from that data alone.
	"""Parameters:
		region_date_dict: Dictionary with isolate counts by region and date
	"""
	@staticmethod
	def region_date_dict_to_df(region_date_dict, db = "bvbrc"):
		region_date_df = []
		for region in region_date_dict.keys():
			for date in region_date_dict[region].keys():
				count = region_date_dict[region][date]
				dto = datetime.strptime(date, '%Y-%m-%d').date()
				data = pd.DataFrame({'Region': [region], 'Date': [dto], 'Count': [count]})
				region_date_df.append(data)

		region_date_df = pd.concat(region_date_df, ignore_index = True)

		with open('data/'+db+'_region_date_counts_df.pkl', 'wb') as handle:
			pickle.dump(region_date_df, handle)

		return region_date_df


	# This function coverts the dictionary containing all the epidemiological dynamics for each covariate
	# or AA mutation into a dataframe format that can be used in downstream data analysis.  This function
	# also saves the dataframe to a pickle.
	"""Parameters:
		prevalence_dict: Dictionary of either covariate or aa mutation region-date prevalence
		mut_type: either 'aa' for amino acid mutation or 'cov' for covariate
	"""
	@staticmethod
	def variant_dict_to_df(prevalence_dict, mut_type):
		prevalence_df = []
		for variant in prevalence_dict.keys():
			for region in prevalence_dict[variant].keys():
				for date in prevalence_dict[variant][region].keys():
					prev = prevalence_dict[variant][region][date]
					dto = datetime.strptime(date, '%Y-%m-%d').date()
					data = pd.DataFrame({'Variant': [variant], 'Region': [region], 'Date': [dto], 'Count': [prev]})
					prevalence_df.append(data)

		prevalence_df = pd.concat(prevalence_df, ignore_index = True)

		with open('data/'+mut_type+'_prevalence_df.pkl', 'wb') as handle:
			pickle.dump(prevalence_df, handle)

		return prevalence_df


	# This function relies on the covariate prevalence output from the spatial_temporal_prevalence function 
	# (a dictionary/json file of covariate date-region prevalence) and computes the prevalence of each amino 
	# acid substituion by time and region.  This function does not need to be called as it is mostly an alternative
	# approach to computing single amino acid substitution dynamics.
	"""Parameters:
		cov_prevalence: Dictionary mapping covariates to region-date prevalence
	"""
	@staticmethod
	def aa_prevalence_extra(cov_prevalence):

		aa_prevalence = {}

		for cov in cov_prevalence.keys():
			aa_muts = cov.split(",")
			for region in cov_prevalence[cov].keys():
				for date in cov_prevalence[cov][region].keys():
					for mut in aa_muts:
						if aa_prevalence.get(mut):
							if aa_prevalence[mut].get(region):
								if aa_prevalence[mut][region].get(date):
									aa_prevalence[mut][region][date] += cov_prevalence[cov][region][date]
								else:
									aa_prevalence[mut][region][date] = cov_prevalence[cov][region][date]
							else:
								aa_prevalence[mut][region] = {}
								aa_prevalence[mut][region][date] = cov_prevalence[cov][region][date]
						else:
							aa_prevalence[mut] = {}
							aa_prevalence[mut][region] = {}
							aa_prevalence[mut][region][date] = cov_prevalence[cov][region][date]

		return aa_prevalence


	# This function gets called within the spatial_temporal_prevalence function to compute the
	# region-date prevalence counts of protein specific covariates.  It returns a dictionary of 
	# dictionaries, containing the region-date covariate counts for each unique covariate of each
	# protein.
	"""Parameters:
		prot_cov_prevalence: Dictionary of region-date protein-specific covariate counts
		prot_cov: Dictionary of the covariate sequence for each protein orginating from some
			single covariate
		region: Region of isolation for sequence that mutations originate
		date: Date of isolation for sequence that mutation orginate
	"""
	@staticmethod
	def protein_covariate_prevalence(prot_cov_prevalence, prot_cov, region, date):

		for prot in prot_cov.keys():
			cov = prot_cov[prot]
			if prot_cov_prevalence.get(prot):
				if prot_cov_prevalence[prot].get(cov):
					if prot_cov_prevalence[prot][cov].get(region):
						if prot_cov_prevalence[prot][cov][region].get(date):
							prot_cov_prevalence[prot][cov][region][date] += 1
						else:
							prot_cov_prevalence[prot][cov][region][date] = 1
					else:
						prot_cov_prevalence[prot][cov][region] = {}
						prot_cov_prevalence[prot][cov][region][date] = 1
				else:
					prot_cov_prevalence[prot][cov] = {}
					prot_cov_prevalence[prot][cov][region] = {}
					prot_cov_prevalence[prot][cov][region][date] = 1
			else:
				prot_cov_prevalence[prot] = {}
				prot_cov_prevalence[prot][cov] = {}
				prot_cov_prevalence[prot][cov][region] = {}
				prot_cov_prevalence[prot][cov][region][date] = 1

		return prot_cov_prevalence



	# This function gets called from within spatial_temporal_prevalence to compute the region-date
	# prevalence of single amino acid substitions in parallel with covariates and is the main method
	# for capturing single amino acid substitution dynamics.
	"""Parameters:
		aa_prevalence: Dictionariy mapping aa mutation to region-date prevalence
		aa_muts: List of mutations from a single covariate
		region: Region of isolation for sequence that mutations originate
		date: Date of isolation for sequence that mutation orginate
	"""
	@staticmethod
	def aa_prevalence_main(aa_prevalence, aa_muts, region, date):

		for mut in aa_muts:
			if aa_prevalence.get(mut):
				if aa_prevalence[mut].get(region):
					if aa_prevalence[mut][region].get(date):
						aa_prevalence[mut][region][date] += 1
					else:
						aa_prevalence[mut][region][date] = 1
				else:
					aa_prevalence[mut][region] = {}
					aa_prevalence[mut][region][date] = 1
			else:
				aa_prevalence[mut] = {}
				aa_prevalence[mut][region] = {}
				aa_prevalence[mut][region][date] = 1

		return aa_prevalence


	# This will compute the covariates for each sequence in the seqs_fasta file and then increment the 
	# prevalence of each covariate by date and region.  It takes the Wuhan reference in fasta and some sequence
	# dump containing potential variants also in fasta.  It also takes optional files that contain previous records
	# including a dictionary of sequence to covariate mappings, a dictionary of covariate prevalence, a dictionary
	# of aa mut prevalence, and a dictionary of QC results for each ViPR/Genbank Sequence.  These four optional 
	# files must be saved as a compressed json (gzip) and will re-save the updated data back to a pickle.
	"""Parameters:
		ref_fasta: SARS-CoV-2 reference fasta
		seqs_fasta: SARS-CoV-2 query sequences fasta
		variant_mappings: Dictionary mapping all sequences to covariate, protein specific covariates, non-synonymous mutations
			synoymous mutations, and ambiguitiy location by protein
		cov_prevalence: Dictionary mapping covariates to region-date prevalence
		aa_prevalence: Dictionary mapping aa mutation to region-date prevalence
		prot_cov_prevalence: Dictionary mapping each protein to covariate region-date prevalence
			region_date_counts 
		region_date_counts: Dictionary of sequence isolate counts by region and date
		qc_storage_file: Dictionary mapping GenBank accesion to sequence to QC status
	"""
	@staticmethod
	def spatial_temporal_prevalence(ref_fasta, seqs_fasta, variant_mappings = {}, cov_prevalence = {}, 
		aa_prevalence = {}, prot_cov_prevalence = {}, region_date_counts = {}, qc_storage_file = {}):

		print("Running alignments, computing spatial-temporal counts ...")

		if variant_mappings != {}:
			variant_mappings = compress_json.load(variant_mappings)
		if cov_prevalence != {}:
			cov_prevalence = compress_json.load(cov_prevalence)
		if aa_prevalence != {}:
			aa_prevalence = compress_json.load(aa_prevalence)
		if prot_cov_prevalence != {}:
			prot_cov_prevalence = compress_json.load(prot_cov_prevalence)
		if region_date_counts != {}:
			region_date_counts = compress_json.load(region_date_counts)	
		if qc_storage_file != {}:
			qc_storage_file = compress_json.load(qc_storage_file)
			new_storage = False
		else:
			new_storage = True

		
		for seq_record in SeqIO.parse(ref_fasta, 'fasta'):
			ref = seq_record.seq

		for seq_record in SeqIO.parse(seqs_fasta, 'fasta'):
			acc = str(seq_record.id.split("|")[1])
			date = str(seq_record.id.split("|")[2])
			date = date.replace("_", "-")
			region = str(seq_record.id.split("|")[4])
			seq = str(seq_record.seq)
			if region == "USA":
				seq_id = str(seq_record.id.split("|")[0])
				loc_id = str(seq_id.split("/")[3])
				state = loc_id.split("_")[0]
				if (state.isalpha() and len(state) == 2):
					region = str("USA_"+state)
					if new_storage:
						qc_storage_file[acc] = {}
						qc_storage_file[acc]['sequence'] = seq
						qc_storage_file[acc]['QC'] = 'Pass'
					else:
						qc_storage_file[acc]['QC'] = 'Pass'
				else:
					if new_storage:
						qc_storage_file[acc] = {}
						qc_storage_file[acc]['sequence'] = seq
						qc_storage_file[acc]['QC'] = 'Fail: US Region'
					else:
						qc_storage_file[acc]['QC'] = 'Fail: US Region'
					continue
			if not(variant_mappings.get(seq)):
				all_mutations = stc.SeqToCovariate.mutations(ref, seq)
				if all_mutations == 'BAD_SEQUENCE':
					if new_storage:
						qc_storage_file[acc] = {}
						qc_storage_file[acc]['sequence'] = seq
						qc_storage_file[acc]['QC'] = 'Fail: Frameshift'
					else:
						qc_storage_file[acc]['QC'] = 'Fail: Frameshift'
					continue
				if new_storage:
					qc_storage_file[acc] = {}
					qc_storage_file[acc]['sequence'] = seq
					qc_storage_file[acc]['QC'] = 'Pass'
				else:
					qc_storage_file[acc]['QC'] = 'Pass'
				non_syn = all_mutations[0]
				syn = all_mutations[1]
				cov = all_mutations[2]
				prot_cov = all_mutations[3]
				ambig = all_mutations[4]
				aa_muts = cov.split(",")
				aa_prevalence = VariantDynamics.aa_prevalence_main(aa_prevalence, aa_muts, region, date)
				prot_cov_prevalence = VariantDynamics.protein_covariate_prevalence(prot_cov_prevalence, prot_cov, region, date)
				if cov_prevalence.get(cov):
					if cov_prevalence[cov].get(region):
						if cov_prevalence[cov][region].get(date):
							cov_prevalence[cov][region][date] += 1
						else:
							cov_prevalence[cov][region][date] = 1
					else:
						cov_prevalence[cov][region] = {}
						cov_prevalence[cov][region][date] = 1
				else:
					cov_prevalence[cov] = {}
					cov_prevalence[cov][region] = {}
					cov_prevalence[cov][region][date] = 1
					variant_mappings[seq] = {}
					variant_mappings[seq]['covariate'] = cov
					variant_mappings[seq]['protein-covariate'] = prot_cov
					variant_mappings[seq]['non-synonymous'] = non_syn
					variant_mappings[seq]['synonymous'] = syn
					variant_mappings[seq]['ambiguous'] = ambig
			else:
				cov = variant_mappings[seq]['covariate']
				prot_cov = variant_mappings[seq]['protein-covariate']
				aa_muts = cov.split(",")
				aa_prevalence = VariantDynamics.aa_prevalence_main(aa_prevalence, aa_muts, region, date)
				prot_cov_prevalence = VariantDynamics.protein_covariate_prevalence(prot_cov_prevalence, prot_cov, region, date)
				if cov_prevalence.get(cov):
					if cov_prevalence[cov].get(region):
						if cov_prevalence[cov][region].get(date):
							cov_prevalence[cov][region][date] += 1
						else:
							cov_prevalence[cov][region][date] = 1
					else:
						cov_prevalence[cov][region] = {}
						cov_prevalence[cov][region][date] = 1
				else:
					cov_prevalence[cov] = {}
					cov_prevalence[cov][region] = {}
					cov_prevalence[cov][region][date] = 1
			if region_date_counts.get(region):
				if region_date_counts[region].get(date):
					region_date_counts[region][date] += 1
				else:
					region_date_counts[region][date] = 1
			else:
				region_date_counts[region] = {}
				region_date_counts[region][date] = 1

		compress_json.dump(variant_mappings, "data/bvbrc_variant_mappings.json.gz")
		compress_json.dump(cov_prevalence, "data/bvrc_cov_prevalence.json.gz")
		compress_json.dump(aa_prevalence, "data/bvbrc_aa_prevalence.json.gz")
		compress_json.dump(prot_cov_prevalence, "data/bvbrc_prot_cov_prevalence.json.gz")
		compress_json.dump(region_date_counts, "data/bvbrc_region_date_counts.json.gz")
		compress_json.dump(qc_storage_file, "data/bvbrc_SARS-CoV-2_Sequence_QC.json.gz")

		print("Done with alignments and region-date counts")
		print('\n')
		
		return (VariantDynamics.variant_dict_to_df(cov_prevalence, 'bvbrc_cov'), 
			VariantDynamics.variant_dict_to_df(aa_prevalence, 'bvbrc_aa'),
			prot_cov_prevalence,
			VariantDynamics.region_date_dict_to_df(region_date_counts))


	# This function computes covariate prevalence for an entire pool of sequences in some fasta file. 
	# It does not compute prevalence by date and region.  It also returns accessions for each covariate
	# as a dictionary.  This function is mainly used for quality control purposes, not a step for analysis
	# in the pipeline.
	"""Parameters:
		ref_fasta: SARS-CoV-2 reference sequence as fasta
		seqs_fasta: SARS-CoV-2 query sequences as fasta
	"""
	@staticmethod
	def covariate_prevalence(ref_fasta, seqs_fasta):

		cov_mappings = {}
		cov_prevalence = {}
		cov_accessions = {}

		for seq_record in SeqIO.parse(ref_fasta, 'fasta'):
			ref = seq_record.seq

		for seq_record in SeqIO.parse(seqs_fasta, 'fasta'):
			acc = seq_record.id.split("|")[1]
			seq = seq_record.seq
			seq_file = open("last_sequence.txt", "w")
			seq_file.write(str(seq))
			seq_file.close()
			if not(cov_mappings.get(seq)):
				all_mutations = stc.SeqToCovariate.mutations(ref, seq)
				if all_mutations == 'BAD_SEQUENCE':
					continue
				cov = all_mutations[2]
				if cov_prevalence.get(cov):
					cov_prevalence[cov] += 1
					cov_accessions[acc] = cov
				else:
					cov_prevalence[cov] = 1
					cov_accessions[acc] = cov
					cov_mappings[seq] = cov
			else:
				cov = cov_mappings[seq]
				cov_prevalence[cov] += 1
				cov_accessions[acc] = cov

		for key in cov_prevalence:
			cov_prevalence[key] = cov_prevalence[key] / len(cov_accessions)

		return cov_prevalence, cov_accessions


if __name__ == "__main__":

	ref = sys.argv[1]
	seqs = sys.argv[2]
	variant_mappings = sys.argv[3]
	cov_prevalence = sys.argv[4]
	aa_prevalence = sys.argv[5]
	region_date_counts = sys.argv[6]
	qc_storage_file = sys.argv[7]

	cov_prev, aa_prev, prot_cov_prev, region_dates = VariantDynamics.spatial_temporal_prevalence(ref, seqs, variant_mappings, cov_prevalence, aa_prevalence, region_date_counts, qc_storage_file)

	for i in cov_prev.keys():
		print(i)
		for j in cov_prev[i].keys():
			print(j)
			for k in cov_prev[i][j].keys():
				print(k, ":", cov_prev[i][j][k])
		print('\n')
